# checks a split has to pass to be trusted
'''
Run this as "python tests/test_splitting.py" for the synthetic cases, which
take seconds, or add "--real shs27k ppi" to also split real benchmarks, which
takes minutes. The synthetic cases are small enough to work out by hand, so a
wrong answer is obvious rather than merely surprising: each one states the
optimum and why, and the solver has to find exactly that. They cover the three
things the method rests on — that similarity closes transitively, that the
solver returns the hand-computed optimum or proves there isn't one, and that a
datapoint survives only when everything it is about agrees on a split. The real
check builds a benchmark's actual graph, splits it, and re-reads the raw
mmseqs, foldseek and rdkit hits to confirm none of them crosses a split.
Retention is only a property of a benchmark when the solver proved it; where it
did not, the number printed is what that run reached in the time it had, and
reruns move a long way.
'''

import argparse
import sys
from collections import Counter

from plug import graphs, splitting
from plug.solver import solve_scip

PASS, FAIL = [], []

def check(name, got, want):
    '''
    record one assertion and print it either way.

    inputs:
        name: str -- what is being checked
        got: object -- what came back
        want: object -- what should have
    outputs:
        none
    '''
    (PASS if got == want else FAIL).append(name)
    print(f"  {'ok  ' if got == want else 'FAIL'}  {name}"
          + ("" if got == want else f"   got {got!r}, want {want!r}"))

def counts(assign, sigs):
    '''
    retained datapoints per split, recounted from an assignment rather than
    read off the solver, so the two can be compared.

    inputs:
        assign: list[int] -- group -> split index
        sigs: dict[frozenset[int], int] -- group set -> datapoints on it
    outputs:
        counts: list[int] -- datapoints kept in each split
    '''
    out = Counter()
    for k, n in sigs.items():
        where = {assign[cid] for cid in k}
        if len(where) == 1: out[where.pop()] += n
    return [out.get(i, 0) for i in range(3)]

def test_groups():
    '''
    similarity has to close transitively, and a datapoint has to be the set of
    groups it rides on however many proteins it is about.

    inputs:
        none
    outputs:
        none
    '''
    print("\ngroups (similarity closes transitively)")

    # a linked to b, b linked to c: all three have to end up in one group even though a
    # and c were never compared. this is the property the whole method rests on.
    got = graphs.merge(["a", "b", "c", "d"], [("a", "b"), ("b", "c")])
    check("a~b~c is one group, d is its own",
          sorted(sorted(g) for g in got), [["a", "b", "c"], ["d"]])

    # nothing linked: every entity is its own group
    check("no links, no groups", len(graphs.merge(list("abcde"))), 5)

    # a datapoint is the SET of groups it rides on — one for a lone protein, two for a
    # pair, and a pair whose ends share a group collapses to one
    class Pair:
        '''a stand-in benchmark whose items are about two proteins.'''

        key, seq, views, ligand = ("a", "b"), ("sa", "sb"), (), "smiles"
        keyfields = classmethod(lambda cls: cls.key)
        seqfields = classmethod(lambda cls: cls.seq)
    comp = {"p1": 0, "p2": 0, "p3": 1}
    sigs = graphs.signatures(Pair, [{"a": "p1", "b": "p2"},      # both in group 0
                                    {"a": "p1", "b": "p3"}], comp)
    check("a pair inside one group is a single-group datapoint",
          sorted((sorted(k), v) for k, v in sigs.items()), [([0], 1), ([0, 1], 1)])

def test_solver():
    '''
    every case here has a hand-computable optimum, including the ones where the
    right answer is that no legal split exists. the proportions are of what
    survives, so a case is infeasible when no assignment can shape the survivors
    into 80/10/10 rather than when it merely loses rows.

    inputs:
        none
    outputs:
        none
    '''
    print("\nsolver (hand-checkable optima)")

    # ten datapoints, each on its own group: 8/1/1 is exactly 80/10/10 of all ten,
    # so everything is retained and the shape is exact
    sigs = {frozenset([i]): 1 for i in range(10)}
    assign, info = solve_scip(sigs, 10, 10, tolerance=0.05, time_limit=30)
    check("ten singletons: sizes are 8/1/1", sorted(counts(assign, sigs)), [1, 1, 8])
    check("ten singletons: nothing dropped", sum(counts(assign, sigs)), 10)
    check("ten singletons: proven optimal", info["proven_optimal"], True)

    # a tolerance that IS reachable must be obeyed, and retention maximised inside it
    sigs = {frozenset([i]): 1 for i in range(100)}
    assign, info = solve_scip(sigs, 100, 100, tolerance=0.02, time_limit=30)
    check("reachable tolerance: sizes are 80/10/10", sorted(counts(assign, sigs)),
          [10, 10, 80])
    check("reachable tolerance: retention is 100", info["retained"], 100)

    # every datapoint on one indivisible group. it has to go somewhere, so the split is
    # 100/0/0 whatever the solver does, and no tolerance around 80/10/10 admits that.
    # the right answer is to prove there is none, not to return a split that misses.
    sigs = {frozenset([0]): 100}
    assign, info = solve_scip(sigs, 1, 100, tolerance=0.05, time_limit=30)
    check("one group: no assignment returned", assign, None)
    check("one group: proven infeasible", info.get("proven_infeasible"), True)

    # 80/10/10 of everything is impossible here, but 80/10/10 of what survives is exact.
    # A carries 8 datapoints alone, B and C one each, and 5 more ride on {A,B} and {A,C}.
    # parting all three drops those 10 shared rows and leaves 8/1/1 — a perfect shape on
    # half the data. keeping A with B instead retains 14 of 15 in one split, which no
    # tolerance admits. so the optimum is 10, and asked for proportions of the ORIGINAL
    # 20 the request cannot be met at all.
    sigs = {frozenset([0]): 8, frozenset([1]): 1, frozenset([2]): 1,
            frozenset([0, 1]): 5, frozenset([0, 2]): 5}
    assign, info = solve_scip(sigs, 3, 20, tolerance=0.02, time_limit=30, basis="retained")
    got = counts(assign, sigs)
    check("of-retained: keeps 10", sum(got), 10)
    check("of-retained: shape is 8/1/1", sorted(got, reverse=True), [8, 1, 1])
    check("of-retained: proven optimal", info["proven_optimal"], True)
    assign, info = solve_scip(sigs, 3, 20, tolerance=0.02, time_limit=30, basis="original")
    check("unreachable of-original: no assignment", assign, None)

    # the warm start is handed over as what it proves: a legal hint that keeps R rows
    # means retained >= R holds at the optimum too, so it is recorded as a floor.
    sigs = {frozenset([i]): 1 for i in range(10)}
    legal = [0] * 8 + [1, 2]                      # 8/1/1, which the tolerance allows
    _, info = solve_scip(sigs, 10, 10, tolerance=0.05, time_limit=30, hint=legal)
    check("legal hint becomes a retention floor", info["hint_floor"], 10)
    # ...and an illegal one bounds nothing, rather than being imposed anyway
    _, info = solve_scip(sigs, 10, 10, tolerance=0.05, time_limit=30, hint=[0] * 10)
    check("illegal hint is ignored", info["hint_floor"], None)

    # the same problem twice has to give the same answer; scip takes no seed and the
    # model is built the same way every time, which is what info["deterministic"] claims
    sigs = {frozenset([i]): i + 1 for i in range(30)}
    kw = dict(tolerance=0.05, time_limit=30)
    a1, i1 = solve_scip(sigs, 30, 465, **kw)
    a2, _ = solve_scip(sigs, 30, 465, **kw)
    check("same problem, same assignment", a1, a2)
    check("scip reports deterministic", i1["deterministic"], True)

class Fake:
    '''a stand-in benchmark whose items pair two proteins with a ligand.'''

    key, seq, views, ligand = ("a", "b"), ("sa", "sb"), ("ligand",), "smiles"
    keyfields = classmethod(lambda cls: cls.key)
    seqfields = classmethod(lambda cls: cls.seq)

def test_stamp_and_audit():
    '''
    a datapoint survives only when everything it is about agrees, and the audit
    has to see a group divided down the middle while not charging the method
    for a hit touching an entity nothing placed.

    inputs:
        none
    outputs:
        none
    '''
    print("\nstamp and audit")

    rows = [{"a": "p1", "b": "p2", "smiles": "L1"},     # all three in train
            {"a": "p1", "b": "p3", "smiles": "L1"},     # p3 in test -> dropped
            {"a": "p1", "b": "p2", "smiles": "L2"},     # ligand in val -> dropped
            {"a": "p1", "b": "p2", "smiles": ""}]       # no ligand, so only the pair counts
    sm = {"p1": "train", "p2": "train", "p3": "test", "L1": "train", "L2": "val"}
    check("a datapoint survives only when everything it is about agrees",
          splitting.stamp(Fake, sm, rows), ["train", None, None, "train"])

    # the audit counts hits whose two ends landed in different splits, and must not charge
    # a method for a hit touching an entity nothing placed — every datapoint about that
    # entity was dropped, so it cannot leak
    g = {"edges": {"protein": [("p1", "p2"), ("p1", "p3"), ("p1", "p4")]},
         "groups": {"protein": [["p1", "p2"], ["p3"], ["p4"]]}}
    a = splitting.audit({"p1": "train", "p2": "train", "p3": "test"}, g)
    check("audit: one hit crosses", a["protein_cross_split"], 1)
    check("audit: the unplaced end is counted separately",
          a["protein_edges_with_an_unplaced_end"], 1)
    check("audit: no group spans splits", a["components_spanning_splits"], 0)

    # a group split down the middle is exactly what must never happen, and the audit has
    # to see it when it does
    b = splitting.audit({"p1": "train", "p2": "test", "p3": "test", "p4": "val"}, g)
    check("audit: a divided group is caught", b["components_spanning_splits"], 1)
    check("audit: totals roll up", b["leaked_edges"], 3)

def test_real(names, rebuild=False):
    '''
    build a benchmark's real graph, split it, and prove nothing crosses. slow —
    the searches and the solve are the same ones a real split runs.

    inputs:
        names: list[str] -- benchmarks to split
        rebuild: bool -- redo the similarity searches
    outputs:
        none
    '''
    for name in names:
        print(f"\n{name} (real data)")
        cls = graphs._resolve(name)
        g = graphs.build(cls, rebuild)
        split_of, info = splitting.split(cls, g, time_limit=120)
        if split_of is None:
            check(f"{name}: a split was found", info.get("status"), "a split")
            continue
        m = splitting.report(cls, split_of, g)
        for mode in g["edges"]:
            check(f"{name}: no {mode} hit crosses a split", m[f"{mode}_cross_split"], 0)
        check(f"{name}: no group spans splits", m["components_spanning_splits"], 0)
        # every entity the graph knows about got placed, so nothing is dropped by accident
        check(f"{name}: every entity is placed",
              sum(1 for mo in g["nodes"] for k in g["nodes"][mo] if k not in split_of), 0)
        # and the retention in the report is what recounting the rows actually gives
        where = splitting.stamp(cls, split_of, g["rows"])
        check(f"{name}: retention recounts identically",
              sum(1 for w in where if w), m["retained"])
        # retention is only a property of the benchmark when the solver proved it. where
        # it did not, the number is what this run reached in the time it had and reruns
        # move a long way: shs27k lands anywhere from 57% to 91% depending on effort, so
        # printing it bare would invite reading a solver budget as a result.
        print(f"        retention {m['retention_pct']:.2f}%  "
              f"{m['counts']['train']:,}/{m['counts']['val']:,}/{m['counts']['test']:,}"
              f"  {'proven optimal' if info.get('proven_optimal') else 'ACHIEVED ONLY'}"
              f" (basis {info.get('balance_basis')},"
              f" {info.get('time_limit_per_stage'):.0f}s stages)")

def main():
    '''
    run the synthetic checks, then the real ones if --real was given, and exit
    non-zero if anything failed.

    inputs:
        none
    outputs:
        none
    '''
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", nargs="*", default=None,
                    help="benchmarks to split for real (slow); no names means all of them")
    ap.add_argument("--rebuild", action="store_true", help="redo the similarity searches")
    a = ap.parse_args()

    test_groups()
    test_solver()
    test_stamp_and_audit()
    if a.real is not None:
        from plug.benchmarks import REGISTRY
        test_real(a.real or list(REGISTRY), a.rebuild)

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    for f in FAIL: print(f"  FAILED: {f}")
    sys.exit(1 if FAIL else 0)

if __name__ == "__main__":
    main()

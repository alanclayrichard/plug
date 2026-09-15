# assign groups to train/val/test, keeping as much data as possible
'''
This is the second half of the method: where the groups go. graphs.py has
already decided what is linked to what and reduced every benchmark to the same
two things, a list of indivisible groups and one signature per datapoint saying
which groups it rides on, so everything here is one optimisation problem. There
is a boolean per group and split and every group goes to exactly one split; a
datapoint counts for a split exactly when every group in its signature is in
that split, which is the leakage rule restated as arithmetic; each split's size
is held within a tolerance of its share; and retained datapoints are what gets
maximised. No group is ever divided, so nothing above threshold can cross a
split whatever the solver decides — leakage is not penalised here, it is
unrepresentable. Because a signature is just a set, a single-protein
measurement, an interaction pair and a protein-ligand affinity are all the same
kind of object by the time they arrive. The model is built in cvxpy and solved
with SCIP, which is the stack DataSAIL uses, so a comparison between the two
methods is a comparison of formulations rather than of solvers.
'''

import sys
import time

from .graphs import SPLITS

FRAC = (0.8, 0.1, 0.1)

def peak_memory_mb():
    '''
    how much memory this process has used at its peak, for the solver report.

    inputs:
        none
    outputs:
        mb: float -- peak resident set size in megabytes
    '''
    import resource
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(kb / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)

def solve_scip(sigs, n_comp, n_items, frac=FRAC, tolerance=0.05, time_limit=600.0,
               basis="retained", verbose=False, hint=None, balance=True,
               strata=None, strat_tolerance=0.75, **_):
    '''
    assign every group, of every modality, in one model. the reified and — z is
    true when every group of a signature is in this split — linearises as z <=
    x for each member plus z >= sum(x) - (size - 1) to force it back up, and
    that second half is not optional: the size tolerance has an upper bound, so
    a solver allowed to leave a z at zero can hide datapoints out of an
    oversized split and satisfy the cap without moving a group. z is continuous
    rather than binary because its integrality is implied by x's, which shrinks
    the tree without changing the relaxation. cvxpy cannot hand SCIP a starting
    solution, so what gets handed over instead is what the hint proves: if that
    assignment is legal and keeps R rows then retained >= R holds at the
    optimum too.

    inputs:
        sigs: Counter[frozenset[int], int] -- group set -> datapoints on it
        n_comp: int -- how many groups there are
        n_items: int -- how many datapoints the benchmark has
        frac: tuple[float, float, float] -- the target train/val/test shares
        tolerance: float -- how far a split may sit from its share
        time_limit: float -- seconds per solver stage
        basis: str -- "original" or "retained"
        verbose: bool -- let SCIP talk
        hint: list[int] | None -- a starting assignment, used as a bound
        balance: bool -- run the second pass that lands the sizes on target
        strata: dict | None -- class -> its signatures, one size rule per class
        strat_tolerance: float -- relative, so 0.75 keeps every class between
                        a quarter and one and three quarters of its share
    outputs:
        assign: list[int] | None -- group -> split index, None if unsolved
        info: dict -- model size, status, retention and peak memory
    '''
    import cvxpy as cp
    import numpy as np
    import scipy.sparse as sp

    keys = list(sigs)
    multi = [k for k in keys if len(k) > 1]
    x = cp.Variable((3, n_comp), boolean=True)
    cons = [cp.sum(x, axis=0) == 1]

    # a datapoint on a single group rides along with it and needs no variable of its own
    w_single = np.zeros(n_comp)
    for k in keys:
        if len(k) == 1: w_single[next(iter(k))] += sigs[k]

    if multi:
        # built with sparse matrices: ppi has half a million of the first kind of row,
        # 123 GB written out in full and 6 MB written sparsely
        z = cp.Variable((3, len(multi)), nonneg=True)
        rows = [(i, cid) for i, k in enumerate(multi) for cid in k]
        r = np.arange(len(rows))
        A = sp.csr_matrix((np.ones(len(rows)), (r, [i for i, _ in rows])),
                          shape=(len(rows), len(multi)))
        B = sp.csr_matrix((np.ones(len(rows)), (r, [cid for _, cid in rows])),
                          shape=(len(rows), n_comp))
        for i in range(3): cons.append(A @ z[i] <= B @ x[i])
        # one row per signature per split, rather than one per member: z >= (how many of
        # its groups are here) - (all but one of them)
        C = sp.csr_matrix((np.ones(len(rows)), ([i for i, _ in rows], [c for _, c in rows])),
                          shape=(len(multi), n_comp))
        size = np.array([len(k) for k in multi], float)
        for i in range(3): cons.append(z[i] >= C @ x[i] - (size - 1))
        w_multi = np.array([sigs[k] for k in multi], float)
        kept = [w_single @ x[i] + w_multi @ z[i] for i in range(3)]
    else:
        kept = [w_single @ x[i] for i in range(3)]

    if basis == "retained":
        R = sum(kept)
        for i in range(3):
            cons += [kept[i] >= (frac[i] - tolerance) * R,
                     kept[i] <= (frac[i] + tolerance) * R]
    else:
        slack = round(tolerance * n_items)
        for i in range(3):
            cons += [kept[i] >= round(frac[i] * n_items) - slack,
                     kept[i] <= round(frac[i] * n_items) + slack]
    # val and test ask for the same share, so every answer has a mirror twin with the two
    # swapped. pinning their order halves the search and cannot cost a datapoint.
    if frac[1] == frac[2]: cons.append(kept[1] >= kept[2])

    # the same sentence, once per class. kept[class][split] is the same linear form as
    # kept[split] with a different weight vector, so this adds three rows per class and
    # not one variable, and it stays a single optimisation rather than a second pass.
    for name, w in (strata or {}).items():
        ws = np.zeros(n_comp)
        wm = np.zeros(len(multi))
        at = {k: i for i, k in enumerate(multi)}
        for k, v in w.items():
            if len(k) == 1: ws[next(iter(k))] += v
            else: wm[at[k]] += v
        kc = [ws @ x[i] + (wm @ z[i] if multi else 0) for i in range(3)]
        Rc = sum(kc)
        for i in range(3):
            cons += [kc[i] >= (1 - strat_tolerance) * frac[i] * Rc,
                     kc[i] <= (1 + strat_tolerance) * frac[i] * Rc]

    floor = None
    if hint is not None:
        k = [0, 0, 0]
        for key, w in sigs.items():
            where = {hint[cid] for cid in key}
            if len(where) == 1: k[where.pop()] += w
        R0 = sum(k)
        # only a legal hint bounds anything. the mirror twin is legal whenever this one is
        # and keeps the same rows, so the val/test order does not have to be checked.
        legal = (all((frac[i] - tolerance) * R0 <= k[i] <= (frac[i] + tolerance) * R0
                     for i in range(3)) if basis == "retained" else
                 all(abs(k[i] - round(frac[i] * n_items)) <= round(tolerance * n_items)
                     for i in range(3)))
        if R0 and legal:
            cons.append(sum(kept) >= R0)
            floor = R0

    t = time.time()
    prob = cp.Problem(cp.Maximize(sum(kept)), cons)
    built = time.time() - t

    info = {"backend": "cvxpy+scip", "components": n_comp, "signatures": len(keys),
            "binaries": 3 * n_comp, "auxiliary_continuous": 3 * len(multi),
            "and_rows": 3 * (len(rows) + len(multi)) if multi else 0, "tolerance": tolerance,
            "balance_basis": basis, "time_limit_per_stage": time_limit,
            "strata": len(strata or {}),
            "strat_tolerance": strat_tolerance if strata else None,
            "build_seconds": round(built, 1), "hint_floor": floor,
            "deterministic": True,
            "objective": "maximise retained datapoints subject to a hard size tolerance"}

    # one scip stage, returning its status and how long it took
    def run(problem):
        t = time.time()
        try:
            problem.solve(solver=cp.SCIP, verbose=verbose,
                          scip_params={"limits/time": float(time_limit)})
            return problem.status, round(time.time() - t, 2)
        except Exception as e:
            # the message, not just the class. cvxpy raises SolverError for everything from
            # "no feasible point inside the limit" to a missing solver, and reporting only
            # the type wrote "SolverError: SolverError" into every summary — true, and
            # useless for telling those apart later.
            return f"{type(e).__name__}: {e}".strip()[:300], round(time.time() - t, 2)

    status, took = run(prob)
    if prob.value is None or not np.isfinite(prob.value):
        info.update(status={"retention": status}, runtime=took,
                    peak_memory_mb=peak_memory_mb(),
                    proven_infeasible="infeasible" in str(status))
        return None, info
    best, answer = int(round(prob.value)), x.value

    # retention rarely has one maximiser: anywhere inside the tolerance box usually keeps
    # the same rows, so left alone the solver returns whichever it found first and the
    # sizes wander — on passerrank, 85/9.5/5.5 while keeping every datapoint either way.
    # this pins the retention just proved and spends a second pass landing the sizes on
    # target, which cannot cost a datapoint. measured against each split's own share, so
    # missing test by half of it costs more than missing train by a sixteenth of it.
    st2, t2 = None, 0.0
    if balance:
        # written as three scalar variables squeezed from both sides rather than with
        # cp.abs. an abs over these expressions is the same arithmetic, but it hands the
        # canonicaliser one term per component per split to rewrite — on the bigger
        # benchmarks that is where the process died, without a traceback to say so.
        R = sum(kept)
        off = cp.Variable(3, nonneg=True)
        shape = cp.Problem(cp.Minimize(cp.sum(off)), cons + [R >= best] + [
            c for i in range(3) for c in
            (off[i] * frac[i] >= kept[i] - frac[i] * R,
             off[i] * frac[i] >= frac[i] * R - kept[i])])
        st2, t2 = run(shape)
        if shape.value is not None and np.isfinite(shape.value): answer = x.value
    info.update(status={"retention": status, "balance": st2},
                runtime=round(took + t2, 2), peak_memory_mb=peak_memory_mb(),
                proven_optimal=(status == "optimal"), retained=best,
                counts=dict(zip(SPLITS, (int(round(k.value)) for k in kept))))
    return [int(np.argmax(answer[:, i])) for i in range(n_comp)], info

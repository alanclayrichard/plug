# leakage-free splits for protein benchmarks, and leakage-free pretraining sets
'''
plug answers one question twice. For a supervised model it asks where a
benchmark's own entries may go, so that no two splits hold things too alike to
be pulled apart; graphs.py decides what is linked to what, solver.py decides
where those groups go, and splitting.py runs the two together and audits the
answer against the raw search hits. For an unsupervised one it asks which
sequences of a reservoir may be trained on at all, which is the same rule
pointed the other way and lives in build_unlabeled_dataset.py. The eleven
benchmarks each say only how to read their own rows; benchmarks/ does the rest.
'''

__version__ = "0.2.0"

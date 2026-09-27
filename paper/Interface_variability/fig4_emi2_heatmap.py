#!/usr/bin/env python
"""Figure 4 — Emi2-native interface map from AF3 S1 + S2b + S3 (no Emi1/4UI9).

Delegates to build_s2b_s3_matrices_and_fig4.py so one script owns the matrices
and the figure.
"""
from build_s2b_s3_matrices_and_fig4 import build_matrices, make_fig4

if __name__ == "__main__":
    combined, n_eng = build_matrices()
    make_fig4(combined, n_eng)

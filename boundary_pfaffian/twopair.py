"""Exact dent ratios in the Aztec diamond with big integers.
S(i,j) = 2^N * R_N(i,j) = sum_{k<N} C(k,i) C(k,j) 2^(N-1-k)   (dents on SW at i, on SE at j, 0-indexed)
Multi-pair ratio M(dents)/M(AD_N) = det[R_N(i_a, j_b)]  (boundary-monomer Pfaffian -> det).
"""
from math import comb, log
import numpy as np

def binom_rows(N, idx):
    return {i: [comb(k, i) for k in range(N)] for i in idx}

def S(N, i, j, B=None):
    bi = B[i] if B else [comb(k, i) for k in range(N)]
    bj = B[j] if B else [comb(k, j) for k in range(N)]
    return sum(bi[k] * bj[k] << (N - 1 - k) for k in range(N))

def lnratio_det(N, I, J, B=None):
    """ln det[R(I,J)] using exact integers (R = S / 2^N)."""
    m = len(I)
    M = [[S(N, i, j, B) for j in J] for i in I]
    from fractions import Fraction
    # exact determinant by Bareiss
    A = [row[:] for row in M]; sign = 1; prev = 1
    for k in range(m - 1):
        if A[k][k] == 0:
            for r in range(k + 1, m):
                if A[r][k] != 0: A[k], A[r] = A[r], A[k]; sign = -sign; break
        for r in range(k + 1, m):
            for c in range(k + 1, m):
                A[r][c] = (A[r][c] * A[k][k] - A[r][k] * A[k][c]) // prev
        prev = A[k][k]
    det = sign * A[m - 1][m - 1]
    return lnbig(det) - m * N * log(2)

def lnbig(x):
    if x <= 0: return float('nan')
    b = x.bit_length(); return log(x >> max(b - 60, 0)) + max(b - 60, 0) * log(2)

"""Leading term of the edge free energy: dg(rho) = kappa rho^2 + O(rho^3), kappa = lim kappa_n = 4/3.

kappa_n = <(1/n) sum_i 1/g_i> / n in the point process mu_n(x) ~ prod_{i<j} |2 sin pi(x_i-x_j)| prod_i g_i on the unit ring
(g_i = nearest-neighbour gaps), i.e. the walls of the first layer pair of the staircase strip in the dilute limit
(see README.md, "dg의 선행 항").  Evaluated exactly for even n by Pfaffians:
    kappa_n = N1 / (2 n^2 D),
    N1 = int dc int_{(n-1)-sets y} |Delta_n(y, c)| prod gaps(y),   D = int_{n-sets x} |Delta_n(x)| prod gaps(x),
    prod gaps = Pfaffian of F(x_i - x_j) with F(u) = u(1-u) (n even) or bordered with F(u) = u(1-u)(1-2u) (n odd),
and de Bruijn's formula for the integral of a determinant times a Pfaffian.
closed(n) = 1 + 2(n^2-1)/(3n) [pi^2/8 - sum_{j<=n/2} (2j-1)^-2] agrees with the Pfaffian value to 57 digits for even n <= 40
(found by integer-relation search for n <= 14; not proved for general n).  Its limit is 4/3.
Usage: python3 dg_leading.py 2 4 6 8      (dg_leading_sampler.c is an independent Monte Carlo check)
"""
import mpmath as mp, itertools, sys
mp.mp.dps = 60
def pf(M):
    M = M.copy(); n = M.rows; r = mp.mpf(1)
    for k in range(0, n - 1, 2):
        p = max(range(k + 1, n), key=lambda i: abs(M[i, k]))
        if abs(M[p, k]) == 0: return mp.mpf(0)
        if p != k + 1:
            for j in range(n): M[k + 1, j], M[p, j] = M[p, j], M[k + 1, j]
            for j in range(n): M[j, k + 1], M[j, p] = M[j, p], M[j, k + 1]
            r = -r
        r *= M[k, k + 1]
        for i in range(k + 2, n):
            t = M[i, k] / M[k + 1, k]          # eliminate M[i,k] using row/col k+1
            if t != 0:
                for j in range(n): M[i, j] -= t * M[k + 1, j]
                for j in range(n): M[j, i] -= t * M[j, k + 1]
    return r
phi0 = lambda k: -1 / (12 * k**2) + 1 / k**4
def S(k, kp):       # sum over q in 2 pi Z \ 0 of 1/(q^3 (k+q)(kp-q)),  k, kp odd multiples of pi
    if k + kp == 0: return 1 / (6 * k**3) - 4 / k**5 + 1 / (4 * k**3)
    return (phi0(k) - phi0(kp)) / (k + kp)
def gapprod(y): 
    y = sorted(y); return mp.fprod([(y[(i + 1) % len(y)] - y[i]) % 1 if i < len(y) - 1 else 1 - y[-1] + y[0] for i in range(len(y))])
def c_odd(m):       # bordered Pfaffian with F = y(1-y)(1-2y) equals c_m * prod(gaps)
    y = [mp.mpf(i * i + 1) / (m * m + 3) for i in range(m)]; M = mp.zeros(m + 1)
    F = lambda u: u * (1 - u) * (1 - 2 * u) if u > 0 else -(-u) * (1 + u) * (1 + 2 * u)
    for i in range(m):
        for j in range(m): M[i, j] = F(y[i] - y[j]) if i != j else 0
        M[i, m] = 1; M[m, i] = -1
    return pf(M) / gapprod(y)
def c_even(n):
    y = [mp.mpf(i * i + 1) / (n * n + 3) for i in range(n)]; M = mp.zeros(n)
    F = lambda u: u * (1 - u) if u > 0 else u * (1 + u)
    for i in range(n):
        for j in range(n): M[i, j] = F(y[i] - y[j]) if i != j else 0
    return pf(M) / gapprod(y)
def kappa(n):
    m = n - 1; ks = [mp.pi * (2 * j - 1) for j in range(-n // 2 + 1, n // 2 + 1)]; tot = mp.mpc(0)
    for a in range(n):
        kb = ks[:a] + ks[a + 1:]; M = mp.zeros(m + 1)
        for i in range(m):
            for j in range(m): M[i, j] = 48j * S(kb[i], kb[j]) if i != j else 0
            M[i, m] = 2j / kb[i]; M[m, i] = -2j / kb[i]
        tot += (-1) ** a * pf(M)
    N1 = abs(tot) / abs(c_odd(m)); D = mp.fprod([4 / k**3 for k in ks if k > 0]) / abs(c_even(n))
    return N1 / (2 * n * n * D)
def closed(n):
    return 1 + mp.mpf(2) * (n * n - 1) / (3 * n) * (mp.pi ** 2 / 8 - mp.fsum(mp.mpf(1) / (2 * j - 1) ** 2 for j in range(1, n // 2 + 1)))
if __name__ == "__main__":
    for n in [int(a) for a in sys.argv[1:]]: print(n, mp.nstr(kappa(n), 20), mp.nstr(abs(kappa(n) - closed(n)), 3), flush=True)

"""Fast checks of the core modules (run: pytest -q tests)."""
import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from scipy.optimize import brentq
from scipy.special import gamma
from elliptic_sol import Sol
from exact_sol import solve
import sample_aztec as SA

def sol_at(x): return Sol(brentq(lambda r: Sol(r).x() - x, 1e-3, 0.995, xtol=1e-13))

def test_table_values():
    # values quoted in the manuscript (Table I)
    for x, sig, ds in ((0.1, 0.2357287, 0.0883), (0.5, 0.7207146, 0.4420), (0.9, 0.9693026, 0.8229)):
        s = sol_at(x); assert abs(s.sigma() - sig) < 2e-7; assert abs(s.dstar() - ds) < 1e-4

def test_series_vs_elliptic():
    # two independent evaluations of the same solution
    s = solve(1 + np.exp(brentq(lambda l: solve(1 + np.exp(l))["x"] - 0.5, -9, 11, xtol=1e-12)), m=8000)
    assert abs(s["sigma"] - sol_at(0.5).sigma()) < 1e-5 and abs(s["dstar"] - sol_at(0.5).dstar()) < 1e-5

def test_small_x_law():
    A = 8 * np.log(2) / np.pi ** 2; B = 4 * np.pi ** 1.5 * np.e / (gamma(0.25) ** 2 * np.log(2))
    for x in (0.01, 0.02): assert abs(sol_at(x).sigma() / (A * x * np.log(B / x)) - 1) < 2e-3
    assert abs(sol_at(0.02).dstar() / (4 * 0.02 * np.log(2) / np.pi) - 1) < 2e-3

def test_aztec_count():
    for N in (2, 3, 5):
        cells, _ = SA.geometry(N); assert abs(SA.lnM(cells, ()) - N * (N + 1) / 2 * np.log(2)) < 1e-9

def test_boundary_pfaffian_vs_enumeration():
    # ln(Z/Z0) from the closed-form boundary kernel against a direct sum over all vacancy sets
    from flint import arb, ctx
    import bpf_exact as BP
    ctx.prec = 200
    for N, z in ((2, 0.3), (3, 0.5)):
        cells, edge = SA.geometry(N); l0 = SA.lnM(cells, ()); Z = 0.0
        for k in range(0, 4 * N + 1, 2):
            for S in itertools.combinations(edge, k):
                l = SA.lnM(cells, S)
                if np.isfinite(l) and l > -30: Z += z ** k * np.exp(l - l0)
        ref = float(BP.lndet(N, BP.build(N, BP.Rmat(N)), arb(z)).mid())
        assert abs(np.log(Z) - ref) < 1e-8

def test_strip_isotherm():
    import tm_diag_b as T
    T.RIN = 0.0; W, t = 10, 0.8
    assert abs(T.densities(W, t ** (W / 2))[2][0] - (2 / np.pi) * np.arctan(t)) < 0.01

def test_edge_free_energy_leading_term():
    # exact Pfaffian value of kappa_n against the closed form, and the limit 4/3
    import mpmath as mp
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "full_vs_boundary"))
    import dg_leading as G
    assert abs(G.kappa(2) - mp.pi ** 2 / 8) < mp.mpf(10) ** -40
    for n in (4, 6, 8): assert abs(G.kappa(n) - G.closed(n)) < mp.mpf(10) ** -40
    assert abs(G.closed(10 ** 5) - mp.mpf(4) / 3) < 1e-9

def test_compact_form_of_the_hodograph_solution():
    # Psi = 1 - i mu - F(J)/(zeta + 1/zeta) with F = sum a_n eps^n (scripts/corner_islands/DERIVATION.md)
    from exact_sol import cheb
    JM = 5.0; c0, c1, a = cheb(JM); n = np.arange(len(a))
    for zeta in (0.3 - 0.05j, 0.7 * np.exp(-0.5j), 0.2 * np.exp(-0.7j), 0.9 * np.exp(-0.2j)):
        J = (zeta**4 + zeta**-4) / 2; B = (c0 - J) / c1; sq = np.sqrt(B * B - 1 + 0j); e = B - sq
        if abs(e) >= 1: e = B + sq
        mu = (1 - zeta**2) / (1 + zeta**2); S = (a * e**n).sum()
        old = np.exp(-1j * np.pi / 4) * np.sqrt(mu) * (1 - J) ** 0.25 * (JM - J) ** 0.5 * (-(2**0.75 / (2 * c1)) * S / (B - e)) + 1 - 1j * mu
        assert abs(old - (1 - 1j * mu - S / (zeta + 1 / zeta))) < 1e-12

def test_differential_relation_sigma_x():
    # d sigma/d x = (2 ln2/pi) K'/K(k1), k1 = 2/p_M, against finite differences of the solution
    import elliptic_sol as ES
    for r in (0.2, 0.5, 0.8):
        h = 1e-4; a, b = ES.Sol(r - h), ES.Sol(r + h)
        assert abs((b.sigma() - a.sigma()) / (b.x() - a.x()) / ES.dsigma_dx(r) - 1) < 1e-6

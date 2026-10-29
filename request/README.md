# Large-scale computation request

Full request document (Korean): https://claude.ai/artifact/7YWL447AGx4L8omreMXHEy

## Question
In the Aztec diamond with monomer fugacity z, is Q_un(N, z) - log2(1/z) a universal (z-independent) function f(N)
for N* = 2 log2(1/z) << N << xi(z)?  Here: f(32) = 1.57-1.69, f(48) = 3.26-3.35, f ~ 4 ln(N/21) (N <= 48).
And what is the screened limit c(z) of Q_un for small z?

## Task A (priority 1)
Parallel tempering in z (`aztec/aztec_pt.c`), z ladder 1e-5 ... 0.1 (must include 1e-5, 3e-5, 1e-4, 3e-4, 1e-3,
3e-3, 1e-2, 3e-2, 0.1), swap acceptance >= 0.2 everywhere, 4 independent seeds per N.

| N | cells | replicas (est.) | sweeps therm + meas | seeds | core-hours |
|---|---|---|---|---|---|
| 32 | 2,112 | 30 | 2e6 + 8e6 | 4 | ~14 |
| 48 | 4,704 | 36 | 4e6 + 1.6e7 | 4 | ~75 |
| 64 | 8,320 | 42 | 6e6 + 2.4e7 | 4 | ~230 |
| 96 | 18,624 | 50 | 1e7 + 4e7 | 4 | ~1,030 |
| 128 (optional) | 33,024 | 60 | 1.6e7 + 6.4e7 | 4 | ~3,500 |

Do N = 32, 48, 64 first (~320 core-hours) and report back.

## Validation (mandatory, before production)
`aztec/validate.sh` must reproduce the exact N = 10 values within 2 sigma.  N = 20 exact references:

| | z=0.001 | 0.003 | 0.01 | 0.03 | 0.1 |
|---|---|---|---|---|---|
| Q_un | 9.2822 | 8.4134 | 7.0521 | 5.4232 | 3.3854 |
| Psi | 0.3994 | 0.3393 | 0.2602 | 0.1836 | 0.1063 |

Errors must come from seed-to-seed scatter: at small z single-run block errors underestimate by 3-5x.

## Deliverables
Output directories + `summary.csv` from `aztec/pt_analyze.py`.

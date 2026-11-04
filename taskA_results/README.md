# Task A results (parallel tempering, N = 32-64)

Imported from branch `claude/new-session-xc9yyg` (commit 6972a1f), produced by a separate session running the
large-scale request in `request/`.  Only reports, scripts and per-run summaries are kept here; the raw per-replica
time series (~100 MB) remain on that branch.

| folder | contents |
|---|---|
| `reports/` | interim reports 1-5 (N = 32, 48 hi/lo z, 64 plateau/hi z), diagnostic N = 32 report, figures |
| `scripts/` | task-A versions of `aztec_pt.c`, `run_pt.sh` (differ from `aztec/`), analysis and queue scripts |
| `summaries/<run>/` | `summary.csv`, z ladders (`*.ladder`, `*.zl`) per run |
| `REQUEST_as_run.md`, `README_branch.md` | the request and README as they stood on the branch |

Main conclusions (see reports): f(N) = Q_un - log2(1/z) data up to N = 64; the 1/2 log2 N hypothesis is rejected by
N = 64; c(0.01) ~ 8.75; Psi(s) collapse.  These are superseded for sigma_inf(x) by the exact boundary Pfaffian
(`boundary_pfaffian/`) and the variational theory, and serve as independent full-monomer checks.

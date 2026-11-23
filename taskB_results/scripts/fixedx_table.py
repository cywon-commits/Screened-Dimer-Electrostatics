#!/usr/bin/env python3
"""Fixed-x table for request B: Q_un - (N/2) sigma_inf(x) at z = 2^(-x N / 2).

usage (from taskB_results/):  python3 scripts/fixedx_table.py RUNDIR [RUNDIR ...] > fixedx_table.csv
  RUNDIR holds N<N>_s<seed>_r<k>.txt from aztec_pt.  Runs ../../scripts/aztec/pt_analyze.py on the whole
  directory and on cold (brick-wall start) / hot (AZTEC_INIT=empty start) subsets.  Seeds whose number
  mod 4 is 1 or 2 are cold starts, 3 or 0 are hot starts (queue_B.sh convention).
Error = seed-to-seed scatter (pt_analyze.py: std of per-seed means / sqrt(#seeds)).
"""
import sys, os, re, csv, glob, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PTA = os.path.normpath(os.path.join(HERE, '..', '..', 'scripts', 'aztec', 'pt_analyze.py'))
SIG = {0.3: 0.52089, 0.4152: 0.64349, 0.5: 0.72071, 0.7: 0.86753}

def analyze(d):
    subprocess.run([sys.executable, PTA, d], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return list(csv.DictReader(open(os.path.join(d, 'summary.csv'))))

def subset(d, seeds, tmp):
    os.makedirs(tmp)
    for f in glob.glob(os.path.join(d, 'N*_s*_r*.txt')):
        if int(re.search(r'_s(\d+)_', os.path.basename(f))[1]) in seeds:
            os.symlink(os.path.abspath(f), os.path.join(tmp, os.path.basename(f)))
    return analyze(tmp)

def pick(rows, z):
    r = min(rows, key=lambda r: abs(float(r['z']) / z - 1))
    return r if abs(float(r['z']) / z - 1) < 1e-4 else None

def min_acc(d):
    m = 1.0
    for f in glob.glob(os.path.join(d, '*.ladder')):
        for l in open(f):
            if l.startswith('swap acceptance'):
                m = min([m] + [float(v) for v in l.split(':')[1].split()])
    return m

w = csv.writer(sys.stdout)
w.writerow(['x', 'N', 'z', 'Q_un', 'err', 'seeds', 'Psi', 'sigma_inf', 'Q_un-(N/2)sigma_inf',
            'Q_un_cold', 'err_cold', 'Q_un_hot', 'err_hot', 'min_swap_acc', 'run'])
for d in sys.argv[1:]:
    d = d.rstrip('/')
    rows = analyze(d)
    seeds = sorted({int(re.search(r'_s(\d+)_', os.path.basename(f))[1]) for f in glob.glob(os.path.join(d, 'N*_s*_r0.txt'))})
    N = int(rows[0]['N'])
    with tempfile.TemporaryDirectory() as t:
        cold = subset(d, [s for s in seeds if s % 4 in (1, 2)], os.path.join(t, 'c'))
        hot = subset(d, [s for s in seeds if s % 4 in (3, 0)], os.path.join(t, 'h'))
    for x, sig in SIG.items():
        z = 2.0 ** (-x * N / 2)
        r = pick(rows, z)
        if r is None: continue
        c, h = pick(cold, z), pick(hot, z)
        q = float(r['Q_un'])
        w.writerow([x, N, f'{z:.10g}', f'{q:.4f}', f"{float(r['err']):.4f}", r['seeds'], f"{float(r['Psi']):.4f}", sig,
                    f'{q - N / 2 * sig:.4f}', f"{float(c['Q_un']):.4f}", f"{float(c['err']):.4f}",
                    f"{float(h['Q_un']):.4f}", f"{float(h['err']):.4f}", f'{min_acc(d):.3f}', os.path.basename(d)])

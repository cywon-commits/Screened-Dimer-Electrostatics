import pickle, glob, re, numpy as np
rows = []
for fn in glob.glob('az/N*_z*.pkl'):
    N = int(re.search(r'N(\d+)', fn)[1]); z = float(re.search(r'_z([\d.e-]+)\.pkl', fn)[1])
    d = pickle.load(open(fn, 'rb'))
    Q = {'top': 0, 'bot': 0, 'left': 0, 'right': 0}; cen = []; fro = []
    for (x, y), (q, n, f) in d.items():
        k = x + y - (N - 1); b = x - y; u = k - N
        if abs(u) > abs(b): w = {('top' if u < 0 else 'bot'): 1}
        elif abs(b) > abs(u): w = {('left' if b < 0 else 'right'): 1}
        else: w = {('top' if u < 0 else 'bot'): 0.5, ('left' if b < 0 else 'right'): 0.5}
        for kk, vv in w.items(): Q[kk] += vv * q
        if abs(u) <= 2 and abs(b) <= 2: cen.append(n)
        fro.append(f)
    fro = np.array(fro)
    rows.append((N, z, Q['top'], Q['bot'], Q['left'], Q['right'], np.mean(cen), np.mean(fro)))
rows.sort()
if __name__ == '__main__':
    print(" N      z      Q_top    Q_bot    Q_left   Q_right   q_e=Q_top/N   n_center  <phi>")
    for r in rows:
        print(f"{r[0]:2d} {r[1]:7.4f} {r[2]:8.4f} {r[3]:8.4f} {r[4]:8.4f} {r[5]:8.4f}   {r[2]/r[0]:.4f}   {r[6]:.5f}  {r[7]:.4f}")
    pickle.dump(rows, open('az/summary.pkl', 'wb'))

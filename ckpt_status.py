#!/usr/bin/env python3
"""Progress and neighbour swap acceptance of running aztec_pt jobs, read from their checkpoints.
usage: ckpt_status.py OUTDIR/*.ckpt"""
import sys, struct, numpy as np
for fn in sys.argv[1:]:
    b = open(fn, 'rb').read()
    N, NR = struct.unpack_from('<ii', b, 0); therm, meas, per = struct.unpack_from('<qqq', b, 8)
    o = 32 + 8 * NR; sw = struct.unpack_from('<q', b, o)[0]; M = 2 * N
    o += 8 + 32 + NR * (4 * M * M + 4 + 8 + 40 * M * M)
    tr = np.frombuffer(b, '<i8', NR, o); ac = np.frombuffer(b, '<i8', NR, o + 8 * NR)
    a = ac[:NR - 1] / np.maximum(tr[:NR - 1], 1)
    print(f"{fn}: sweep {sw:,} / {therm + meas:,} ({100 * sw / (therm + meas):.1f}%), min acc {a.min():.3f}, acc {' '.join(f'{x:.2f}' for x in a)}")

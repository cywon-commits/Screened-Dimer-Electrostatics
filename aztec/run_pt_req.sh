#!/bin/bash
# usage: REQ="z1 z2 ..." run_pt_req.sh N SEED R ZMIN ZMAX THERM MEAS OUTDIR [NBLK] [EVERY]
# Same as run_pt.sh, but the z values listed in $REQ are inserted into the geometric ladder (fixed-x runs).
set -e
N=$1; SEED=$2; R=$3; ZMIN=$4; ZMAX=$5; THERM=$6; MEAS=$7; OUT=$8; NBLK=${9:-100}; EVERY=${10:-100}
mkdir -p "$OUT"
LADDER=$(python3 - "$R" "$ZMIN" "$ZMAX" $REQ <<'PY'
import sys, numpy as np
R, a, b = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]); req = [float(v) for v in sys.argv[4:]]
z = set(np.geomspace(a, b, R).round(12)) | {x for x in req if a <= x <= b}
print(' '.join(f'{x:.8g}' for x in sorted(z)))
PY
)
echo "N=$N seed=$SEED ladder: $LADDER" | tee "$OUT/N${N}_s${SEED}.ladder"
./aztec_pt "$N" "$THERM" "$MEAS" "$NBLK" "$SEED" "$EVERY" "$OUT/N${N}_s${SEED}" $LADDER | tee -a "$OUT/N${N}_s${SEED}.ladder"

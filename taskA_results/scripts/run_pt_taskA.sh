#!/bin/bash
# usage: run_pt.sh N SEED R ZMIN ZMAX THERM MEAS OUTDIR [NBLK] [EVERY]
# R = number of replicas: geometric z ladder of R replicas between ZMIN and ZMAX, plus the required points.
# R = a<T> (e.g. a0.35): ladder from ladder.py with predicted neighbour swap acceptance >= T.
# The ladder is saved in OUTDIR/N<N>_s<SEED>.zl on first start and reused when a run is resumed.
set -e
N=$1; SEED=$2; R=$3; ZMIN=$4; ZMAX=$5; THERM=$6; MEAS=$7; OUT=$8; NBLK=${9:-100}; EVERY=${10:-100}
mkdir -p "$OUT"
ZL="$OUT/N${N}_s${SEED}.zl"
if [ -s "$ZL" ]; then LADDER=$(cat "$ZL")
elif [[ $R == a* ]]; then LADDER=$(python3 "$(dirname "$0")/ladder.py" "$N" "$ZMIN" "$ZMAX" "${R#a}")
else LADDER=$(python3 - "$R" "$ZMIN" "$ZMAX" <<'PY'
import sys, numpy as np
R, a, b = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
req = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1]
z = set(np.geomspace(a, b, R).round(12)) | {x for x in req if a <= x <= b}
print(' '.join(f'{x:.6g}' for x in sorted(z)))
PY
)
fi
echo "$LADDER" > "$ZL"
echo "N=$N seed=$SEED ladder: $LADDER" | tee "$OUT/N${N}_s${SEED}.ladder"
./aztec_pt "$N" "$THERM" "$MEAS" "$NBLK" "$SEED" "$EVERY" "$OUT/N${N}_s${SEED}" $LADDER | tee -a "$OUT/N${N}_s${SEED}.ladder"

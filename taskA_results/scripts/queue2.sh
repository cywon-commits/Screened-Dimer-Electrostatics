#!/bin/bash
# Revised plan after the N=32 diagnostic (reports/diag_N32.md): split ladders.
#   N=48 large z [1e-3, 0.1]: short (tau ~ 3e4 sweeps), 4 seeds            -> out48hi/
#   N=48 small z [1e-5, 1e-3]: long, 2 cold (brick-wall) + 2 hot (all-monomer) starts,
#                              recorded from sweep 0 (discard the start in the analysis) -> out48lo/
# At most NPAR jobs at once. Safe to rerun after a restart: finished jobs are skipped, interrupted ones resume.
# usage: setsid nohup ./queue2.sh >> queue2.log 2>&1 < /dev/null &
cd "$(dirname "$0")"
NPAR=${NPAR:-$(nproc)}
JOBS=(
  "cold|48 1 a0.35 1e-3 0.1 1000000 3000000 out48hi 150 200"
  "cold|48 2 a0.35 1e-3 0.1 1000000 3000000 out48hi 150 200"
  "cold|48 3 a0.35 1e-3 0.1 1000000 3000000 out48hi 150 200"
  "cold|48 4 a0.35 1e-3 0.1 1000000 3000000 out48hi 150 200"
  "cold|48 1 a0.45 1e-5 1e-3 0 60000000 out48lo 1500 400"
  "cold|48 2 a0.45 1e-5 1e-3 0 60000000 out48lo 1500 400"
  "empty|48 3 a0.45 1e-5 1e-3 0 60000000 out48lo 1500 400"
  "empty|48 4 a0.45 1e-5 1e-3 0 60000000 out48lo 1500 400"
)
done_job(){ grep -q '^swap acceptance' "$3/N$1_s$2.ladder" 2>/dev/null && [ ! -e "$3/N$1_s$2.ckpt" ]; }
running_job(){ pgrep -f "aztec_pt $1 .* $3/N$1_s$2 " >/dev/null; }
for j in "${JOBS[@]}"; do
  init=${j%%|*}; args=${j#*|}
  set -- $args
  if done_job $1 $2 $8; then echo "$(date +%F\ %T) skip $8 N=$1 s=$2 (done)"; continue; fi
  if running_job $1 $2 $8; then echo "$(date +%F\ %T) $8 N=$1 s=$2 already running"; continue; fi
  while [ "$(pgrep -c -x aztec_pt)" -ge "$NPAR" ]; do sleep 60; done
  echo "$(date +%F\ %T) start $8 N=$1 s=$2 init=$init"
  if [ "$init" = empty ]; then AZTEC_INIT=empty ./run_pt.sh $args >> "${8}_N$1_s$2.log" 2>&1 < /dev/null &
  else ./run_pt.sh $args >> "${8}_N$1_s$2.log" 2>&1 < /dev/null & fi
  sleep 5
done
wait
while pgrep -x aztec_pt >/dev/null; do sleep 60; done
echo "$(date +%F\ %T) queue2 finished"

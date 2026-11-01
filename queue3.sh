#!/bin/bash
# N=64 after the N=48 report: plateau ladder for H1 (f(64)), then large z for c(z).
# At most NPAR jobs at once. Safe to rerun after a restart: finished jobs are skipped, interrupted ones resume.
#   out64pl: z in [5e-5, 1e-3], 2 cold + 2 hot starts, recorded from sweep 0 (discard by cold/hot convergence)
#   out64hi: z in [1e-3, 0.1], 4 seeds, short
# usage: setsid nohup ./queue3.sh >> queue3.log 2>&1 < /dev/null &
cd "$(dirname "$0")"
NPAR=${NPAR:-$(nproc)}
JOBS=(
  "cold|64 1 a0.45 5e-5 1e-3 0 30000000 out64pl 750 400"
  "cold|64 2 a0.45 5e-5 1e-3 0 30000000 out64pl 750 400"
  "empty|64 3 a0.45 5e-5 1e-3 0 30000000 out64pl 750 400"
  "empty|64 4 a0.45 5e-5 1e-3 0 30000000 out64pl 750 400"
  "cold|64 1 a0.35 1e-3 0.1 1000000 3000000 out64hi 150 200"
  "cold|64 2 a0.35 1e-3 0.1 1000000 3000000 out64hi 150 200"
  "cold|64 3 a0.35 1e-3 0.1 1000000 3000000 out64hi 150 200"
  "cold|64 4 a0.35 1e-3 0.1 1000000 3000000 out64hi 150 200"
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
echo "$(date +%F\ %T) queue3 finished"

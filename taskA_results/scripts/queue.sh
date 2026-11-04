#!/bin/bash
# Production queue for task A (priority N=32, 48, 64; 4 seeds each), at most NPAR jobs at once.
# Safe to rerun after a restart: finished jobs are skipped, interrupted ones resume from their checkpoint.
# Ladders (R=a0.35) come from ladder.py; N=64 waits until all N=48 runs are done and uses their n(z) as reference.
# Last: N=32 again with the acceptance-designed ladder (out32b), since the first N=32 ladder had acceptance < 0.2 at z > 0.03.
# usage: setsid nohup ./queue.sh > queue.log 2>&1 &
cd "$(dirname "$0")"
NPAR=${NPAR:-$(nproc)}
JOBS=(
  "32 1 30 1e-5 0.1 2000000 8000000 out 200 200"
  "32 2 30 1e-5 0.1 2000000 8000000 out 200 200"
  "32 3 30 1e-5 0.1 2000000 8000000 out 200 200"
  "32 4 30 1e-5 0.1 2000000 8000000 out 200 200"
  "48 1 a0.35 1e-5 0.1 4000000 16000000 out 200 200"
  "48 2 a0.35 1e-5 0.1 4000000 16000000 out 200 200"
  "48 3 a0.35 1e-5 0.1 4000000 16000000 out 200 200"
  "48 4 a0.35 1e-5 0.1 4000000 16000000 out 200 200"
  "64 1 a0.35 1e-5 0.1 6000000 24000000 out 200 200"
  "64 2 a0.35 1e-5 0.1 6000000 24000000 out 200 200"
  "64 3 a0.35 1e-5 0.1 6000000 24000000 out 200 200"
  "64 4 a0.35 1e-5 0.1 6000000 24000000 out 200 200"
  "32 1 a0.35 1e-5 0.1 2000000 8000000 out32b 200 200"
  "32 2 a0.35 1e-5 0.1 2000000 8000000 out32b 200 200"
  "32 3 a0.35 1e-5 0.1 2000000 8000000 out32b 200 200"
  "32 4 a0.35 1e-5 0.1 2000000 8000000 out32b 200 200"
)
done_job(){ grep -q '^swap acceptance' "$3/N$1_s$2.ladder" 2>/dev/null && [ ! -e "$3/N$1_s$2.ckpt" ]; }
running_job(){ pgrep -f "aztec_pt $1 .* $3/N$1_s$2 " >/dev/null; }
for j in "${JOBS[@]}"; do
  set -- $j
  if done_job $1 $2 $8; then echo "$(date +%F\ %T) skip N=$1 s=$2 (done)"; continue; fi
  if running_job $1 $2 $8; then echo "$(date +%F\ %T) N=$1 s=$2 already running"; continue; fi
  if [ "$1" = 64 ] && [ ! -s "$8/N64_s$2.zl" ] && [ ! -e ladder_ref_N48.txt ]; then
    for s in 1 2 3 4; do while ! done_job 48 $s out; do sleep 60; done; done
    python3 ladder.py ref out 48 && echo "$(date +%F\ %T) wrote ladder_ref_N48.txt"
  fi
  while [ "$(pgrep -c -x aztec_pt)" -ge "$NPAR" ]; do sleep 60; done
  echo "$(date +%F\ %T) start N=$1 s=$2"
  ./run_pt.sh $j > "$8_N$1_s$2.log" 2>&1 < /dev/null &
  sleep 5
done
wait
while pgrep -x aztec_pt >/dev/null; do sleep 60; done
echo "$(date +%F\ %T) queue finished"

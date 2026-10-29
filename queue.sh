#!/bin/bash
# Production queue for task A (priority N=32, 48, 64; 4 seeds each), at most NPAR jobs at once.
# Safe to rerun after a restart: finished jobs are skipped, interrupted ones resume from their checkpoint.
# usage: setsid nohup ./queue.sh > queue.log 2>&1 &
cd "$(dirname "$0")"
NPAR=${NPAR:-$(nproc)}
JOBS=(
  "32 1 30 1e-5 0.1 2000000 8000000 out 200 200"
  "32 2 30 1e-5 0.1 2000000 8000000 out 200 200"
  "32 3 30 1e-5 0.1 2000000 8000000 out 200 200"
  "32 4 30 1e-5 0.1 2000000 8000000 out 200 200"
  "48 1 36 1e-5 0.1 4000000 16000000 out 200 200"
  "48 2 36 1e-5 0.1 4000000 16000000 out 200 200"
  "48 3 36 1e-5 0.1 4000000 16000000 out 200 200"
  "48 4 36 1e-5 0.1 4000000 16000000 out 200 200"
  "64 1 42 1e-5 0.1 6000000 24000000 out 200 200"
  "64 2 42 1e-5 0.1 6000000 24000000 out 200 200"
  "64 3 42 1e-5 0.1 6000000 24000000 out 200 200"
  "64 4 42 1e-5 0.1 6000000 24000000 out 200 200"
)
done_job(){ grep -q '^swap acceptance' "out/N$1_s$2.ladder" 2>/dev/null && [ ! -e "out/N$1_s$2.ckpt" ]; }
running_job(){ pgrep -f "aztec_pt $1 .* out/N$1_s$2 " >/dev/null; }
for j in "${JOBS[@]}"; do
  set -- $j
  if done_job $1 $2; then echo "$(date +%F\ %T) skip N=$1 s=$2 (done)"; continue; fi
  if running_job $1 $2; then echo "$(date +%F\ %T) N=$1 s=$2 already running"; continue; fi
  while [ "$(pgrep -c -x aztec_pt)" -ge "$NPAR" ]; do sleep 60; done
  echo "$(date +%F\ %T) start N=$1 s=$2"
  ./run_pt.sh $j > "out_N$1_s$2.log" 2>&1 < /dev/null &
  sleep 5
done
wait
while pgrep -x aztec_pt >/dev/null; do sleep 60; done
echo "$(date +%F\ %T) queue finished"

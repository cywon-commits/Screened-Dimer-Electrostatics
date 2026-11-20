#!/bin/bash
# Request B (request/REQUEST_B_ko.md on main): fixed-x series, z = 2^(-x N / 2).
# Job = "init|extra z on the ladder|run_pt.sh arguments". At most NPAR jobs at once; rerun after a restart to resume.
# usage: setsid nohup ./queue4.sh >> queue4.log 2>&1 < /dev/null &
cd "$(dirname "$0")"
NPAR=${NPAR:-$(nproc)}
VAL="0.0625 0.020617311 0.03125 0.0078125"
JOBS=(
  "cold|$VAL|16 1 a0.35 0.005 0.125 500000 4000000 valB 100 100"
  "cold|$VAL|16 2 a0.35 0.005 0.125 500000 4000000 valB 100 100"
  "cold|$VAL|16 3 a0.35 0.005 0.125 500000 4000000 valB 100 100"
  "cold|$VAL|16 4 a0.35 0.005 0.125 500000 4000000 valB 100 100"
  "cold|$VAL|20 1 a0.35 0.005 0.125 500000 4000000 valB 100 100"
  "cold|$VAL|20 2 a0.35 0.005 0.125 500000 4000000 valB 100 100"
  "cold|$VAL|20 3 a0.35 0.005 0.125 500000 4000000 valB 100 100"
  "cold|$VAL|20 4 a0.35 0.005 0.125 500000 4000000 valB 100 100"
)
done_job(){ grep -q '^swap acceptance' "$3/N$1_s$2.ladder" 2>/dev/null && [ ! -e "$3/N$1_s$2.ckpt" ]; }
running_job(){ pgrep -f "aztec_pt $1 .* $3/N$1_s$2 " >/dev/null; }
for j in "${JOBS[@]}"; do
  IFS='|' read -r init extra args <<< "$j"
  set -- $args
  if done_job $1 $2 $8; then echo "$(date +%F\ %T) skip $8 N=$1 s=$2 (done)"; continue; fi
  if running_job $1 $2 $8; then echo "$(date +%F\ %T) $8 N=$1 s=$2 already running"; continue; fi
  while [ "$(pgrep -c -x aztec_pt)" -ge "$NPAR" ]; do sleep 30; done
  echo "$(date +%F\ %T) start $8 N=$1 s=$2 init=$init"
  if [ "$init" = empty ]; then AZTEC_INIT=empty LADDER_EXTRA="$extra" ./run_pt.sh $args >> "${8}_N$1_s$2.log" 2>&1 < /dev/null &
  else LADDER_EXTRA="$extra" ./run_pt.sh $args >> "${8}_N$1_s$2.log" 2>&1 < /dev/null & fi
  sleep 3
done
wait
while pgrep -x aztec_pt >/dev/null; do sleep 30; done
echo "$(date +%F\ %T) queue4 finished"

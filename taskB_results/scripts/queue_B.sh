#!/bin/bash
# Request B (request/REQUEST_B_ko.md): fixed x = 2 log2(1/z)/N, z = 2^(-x N/2), several x on one ladder per N.
# Binary: taskA_results/scripts/aztec_pt_taskA.c (checkpoint/resume, AZTEC_INIT=empty for the all-monomer start).
# Job = "init|z values that must sit on the ladder|run_pt_B.sh arguments". At most NPAR jobs at once.
# Run from taskB_results/run/ (holds ./aztec_pt); rerun after a restart to resume from checkpoints.
# usage: cd taskB_results/run && setsid nohup ../scripts/queue_B.sh >> queue_B.log 2>&1 < /dev/null &
S="$(cd "$(dirname "$0")" && pwd)"
NPAR=${NPAR:-$(nproc)}
Z32="0.0004250735172 0.00390625 0.03589682359"   # x = 0.7, 0.5, 0.3
Z48="0.000244140625 0.006801176276"              # x = 0.5, 0.3
Z24="0.002960383919"                             # x = 0.7
Z64a="0.001288581944"                            # x = 0.3
Z64b="1.525878906e-05"                           # x = 0.5 (slow region: plateau-style ladder, low precision)
JOBS=()
for s in 1 2 3 4; do i=cold; [ $s -ge 3 ] && i=empty; JOBS+=("$i|$Z32|32 $s a0.4 0.0004250735172 0.1 2000000 8000000 B32 200 200"); done
for s in 1 2 3 4; do i=cold; [ $s -ge 3 ] && i=empty; JOBS+=("$i|$Z48|48 $s a0.4 0.000244140625 0.1 4000000 16000000 B48 400 200"); done
for s in 1 2 3 4; do i=cold; [ $s -ge 3 ] && i=empty; JOBS+=("$i|$Z24|24 $s a0.4 0.002960383919 0.1 1000000 4000000 B24 100 200"); done
for s in 1 2 3 4; do i=cold; [ $s -ge 3 ] && i=empty; JOBS+=("$i|$Z64a|64 $s a0.4 0.001288581944 0.1 1000000 5000000 B64a 125 200"); done
# B48 seeds 5-8 (added after B48 1-4 gave +-0.026 at x=0.5; target +-0.02)
for s in 5 6 7 8; do i=cold; [ $s -ge 7 ] && i=empty; JOBS+=("$i|$Z48|48 $s a0.4 0.000244140625 0.1 4000000 16000000 B48 400 200"); done
for s in 1 2 3 4; do i=cold; [ $s -ge 3 ] && i=empty; JOBS+=("$i|$Z64b|64 $s a0.45 1.525878906e-05 0.001 0 40000000 B64b 1000 400"); done
done_job(){ grep -q '^swap acceptance' "$3/N$1_s$2.ladder" 2>/dev/null && [ ! -e "$3/N$1_s$2.ckpt" ]; }
running_job(){ pgrep -f "aztec_pt $1 .* $3/N$1_s$2 " >/dev/null; }
for j in "${JOBS[@]}"; do
  IFS='|' read -r init extra args <<< "$j"
  set -- $args
  if done_job $1 $2 $8; then echo "$(date +%F\ %T) skip $8 N=$1 s=$2 (done)"; continue; fi
  if running_job $1 $2 $8; then echo "$(date +%F\ %T) $8 N=$1 s=$2 already running"; continue; fi
  while [ "$(pgrep -c -x aztec_pt)" -ge "$NPAR" ]; do sleep 30; done
  echo "$(date +%F\ %T) start $8 N=$1 s=$2 init=$init"
  if [ "$init" = empty ]; then AZTEC_INIT=empty LADDER_EXTRA="$extra" "$S/run_pt_B.sh" $args >> "${8}_N$1_s$2.log" 2>&1 < /dev/null &
  else LADDER_EXTRA="$extra" "$S/run_pt_B.sh" $args >> "${8}_N$1_s$2.log" 2>&1 < /dev/null & fi
  sleep 3
done
wait
while pgrep -x aztec_pt >/dev/null; do sleep 30; done
echo "$(date +%F\ %T) queue_B finished"

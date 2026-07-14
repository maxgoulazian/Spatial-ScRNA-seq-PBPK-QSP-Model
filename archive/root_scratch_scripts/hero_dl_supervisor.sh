#!/bin/bash
# Accelerated resumable hero D2-D4 download supervisor (16-conn aria2c, relaunch-until-complete).
BASE="https://genome-scale-tcell-perturb-seq.s3.us-east-1.amazonaws.com/marson2025_data"
HERO="/media/balthasar-lab/RAID4/costim_engager_counterscreen/data/hero_cd4_perturbseq_zhu2025"
LOG="/tmp/hero_dl_supervisor.log"
cd "$HERO" || exit 1
STEMS="D2_Stim48hr D3_Rest D3_Stim8hr D3_Stim48hr D4_Rest D4_Stim8hr D4_Stim48hr"
> /tmp/hero_urls.txt
declare -A TGT
for stem in $STEMS; do
  F="${stem}.assigned_guide.h5ad"
  T=$(curl -sI "$BASE/$F" | grep -i content-length | tr -d '\r' | awk '{print $2}')
  TGT[$F]=$T
  printf '%s\n  out=%s\n' "$BASE/$F" "$F" >> /tmp/hero_urls.txt
done
echo "[$(date +%H:%M:%S)] SUPERVISOR START — 7 files queued, targets:" >> $LOG
for F in "${!TGT[@]}"; do echo "    $F = ${TGT[$F]}" >> $LOG; done
for round in $(seq 1 40); do
  ALL=1; PENDING=""
  for F in "${!TGT[@]}"; do
    D=$(stat -c%s "$F" 2>/dev/null || echo 0)
    if [ "$D" != "${TGT[$F]}" ]; then ALL=0; PENDING="$PENDING $F"; fi
  done
  if [ "$ALL" = "1" ]; then echo "[$(date +%H:%M:%S)] ALL 7 COMPLETE — supervisor exiting" >> $LOG; break; fi
  echo "[$(date +%H:%M:%S)] round $round | pending:$PENDING | launching aria2c -j2 -x16" >> $LOG
  aria2c -i /tmp/hero_urls.txt -j2 -x16 -s16 -k1M -c \
    --file-allocation=none --max-tries=8 --retry-wait=15 --timeout=60 \
    --auto-file-renaming=false --console-log-level=warn --summary-interval=60 >> $LOG 2>&1
  sleep 8
done
echo "[$(date +%H:%M:%S)] SUPERVISOR DONE" >> $LOG

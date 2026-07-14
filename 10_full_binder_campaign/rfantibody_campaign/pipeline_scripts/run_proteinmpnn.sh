#!/bin/bash
set -e
mpnn_one(){ # $1=target $2=loop_string
  local t=$1 loops=$2
  local indir=/home/campaign_out/$t
  local outdir=/home/campaign_out/${t}_mpnn
  local n=$(ls $indir/*.pdb 2>/dev/null | wc -l)
  if [ "$n" -eq 0 ]; then echo "!! $t: 0 backbones, SKIP"; return; fi
  mkdir -p $outdir
  echo "=== [$(date +%H:%M:%S)] ProteinMPNN $t: $n backbones, loops=$loops, 8 seqs/struct ==="
  poetry run python /home/scripts/proteinmpnn_interface_design.py \
    -pdbdir $indir -outpdbdir $outdir \
    -seqs_per_struct 8 -loop_string "$loops" \
    -temperature 0.000001 -omit_AAs CX -num_connections 48 \
    -checkpoint_name /home/campaign_out/${t}_mpnn.checkpoint 2>&1 | tail -4
  echo "=== [$(date +%H:%M:%S)] DONE $t: $(ls $outdir/*.pdb 2>/dev/null | wc -l) seq-designed ==="
}
# VHH arms: H-only loops.  VH/VL arms (CD3, CEA5): all 6 CDRs.
mpnn_one 41BB "H1,H2,H3"
mpnn_one CD27 "H1,H2,H3"
mpnn_one CD3  "H1,H2,H3,L1,L2,L3"
mpnn_one CEA5 "H1,H2,H3,L1,L2,L3"
echo "=== ProteinMPNN COMPLETE $(date +%H:%M:%S) ==="

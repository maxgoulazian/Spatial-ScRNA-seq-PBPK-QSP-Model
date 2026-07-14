#!/bin/bash
# $1=target  $2=loop_string
t=$1; loops=$2
indir=/home/campaign_out/$t
outdir=/home/campaign_out/${t}_mpnn
n=$(ls $indir/${t}_*.pdb 2>/dev/null | grep -vE 'traj|Xt|pX0' | wc -l)
mkdir -p $outdir
echo "=== [$(date +%H:%M:%S)] ProteinMPNN(CPU) $t: $n backbones, loops=$loops, 10 seqs/struct ==="
poetry run python /home/scripts/proteinmpnn_interface_design.py \
  -pdbdir $indir -outpdbdir $outdir \
  -seqs_per_struct 10 -loop_string "$loops" \
  -temperature 0.000001 -omit_AAs CX -num_connections 48 \
  -checkpoint_name /home/campaign_out/${t}_mpnn_cpu.checkpoint
echo "=== [$(date +%H:%M:%S)] DONE $t: $(ls $outdir/*.pdb 2>/dev/null | wc -l) seq-designed ==="

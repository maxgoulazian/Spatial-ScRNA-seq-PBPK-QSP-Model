#!/bin/bash
set -e
run_one(){ # $1=name $2=target $3=hotspots
  echo "=== [$(date +%H:%M:%S)] START $1 (100 designs) ==="
  poetry run python /home/scripts/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb=/home/input/$2 \
    antibody.framework_pdb=/home/input/framework_VHH.pdb \
    inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt \
    "ppi.hotspot_res=$3" \
    "antibody.design_loops=[H1:7,H2:6,H3:5-13]" \
    inference.num_designs=100 \
    inference.deterministic=True \
    inference.output_prefix=/home/campaign_out/$1/$1 2>&1 | grep -E "Finished design|START|Timestep 2" | tail -5
  echo "=== [$(date +%H:%M:%S)] DONE $1: $(ls /home/campaign_out/$1/*.pdb 2>/dev/null | wc -l) backbones ==="
}
run_one 41BB target_41BB.pdb "[T26,T31,T40,T41,T44]"
run_one CD27 target_CD27.pdb "[T72,T80,T83,T85,T113,T115]"
echo "=== CAMPAIGN COMPLETE $(date +%H:%M:%S) ==="

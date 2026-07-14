#!/bin/bash
set -e
echo "=== [$(date +%H:%M:%S)] START CD3 VH/VL (100 designs) ==="
poetry run python /home/scripts/rfdiffusion_inference.py \
  --config-name antibody \
  antibody.target_pdb=/home/input/target_CD3.pdb \
  antibody.framework_pdb=/home/input/framework_VHVL.pdb \
  inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt \
  "ppi.hotspot_res=[T186,T141,T189,T188,T155,T187]" \
  "antibody.design_loops=[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]" \
  inference.num_designs=100 \
  inference.deterministic=True \
  inference.output_prefix=/home/campaign_out/CD3/CD3 2>&1 | grep -E "Finished design|START|Timestep 2" | tail -5
echo "=== [$(date +%H:%M:%S)] DONE CD3: $(ls /home/campaign_out/CD3/*.pdb 2>/dev/null | wc -l) backbones ==="

#!/bin/bash
set -e
echo "=== [$(date +%H:%M:%S)] START CEA5 B3 single-domain VH/VL (100 designs) ==="
poetry run python /home/scripts/rfdiffusion_inference.py \
  --config-name antibody \
  antibody.target_pdb=/home/input/target_CEA5_B3.pdb \
  antibody.framework_pdb=/home/input/framework_VHVL.pdb \
  inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt \
  "ppi.hotspot_res=[T607,T630,T642,T644,T648,T675]" \
  "antibody.design_loops=[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]" \
  inference.num_designs=100 \
  inference.deterministic=True \
  inference.output_prefix=/home/campaign_out/CEA5/CEA5 2>&1 | grep -E "Finished design|START|Timestep 2" | tail -5
echo "=== [$(date +%H:%M:%S)] DONE CEA5: $(ls /home/campaign_out/CEA5/*.pdb 2>/dev/null | wc -l) backbones ==="

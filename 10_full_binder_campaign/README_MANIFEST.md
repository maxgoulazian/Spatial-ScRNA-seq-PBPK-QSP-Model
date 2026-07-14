# Binder Design + Animation — Complete Local Archive

Everything from the de novo binder-design campaign (RFdiffusion → ProteinMPNN → AlphaFold3)
and the 4-arm diffusion timelapse, consolidated into one folder.

Location: /media/balthasar-lab/RAID4/costim_engager_counterscreen/BINDER_AND_ANIMATION_COMPLETE/

## Folder layout

### rfantibody_campaign/   (the full RFantibody pipeline: input → output)
- `campaign_out/`        — 4 targets (41BB, CD27, CD3, CEA5): 100 RFdiffusion backbones each +
                           200 denoising-trajectory PDBs each (traj/). 1,200 PDBs total.
                           Docked winners: 41BB/41BB_62.pdb, CD27/CD27_24.pdb, CD3/CD3_23.pdb, CEA5/CEA5_52.pdb
- `input/`               — target structures + VHH/VH-VL framework inputs (14 PDBs):
                           target_41BB/CD27/CD3/CEA5_A3B3/CEA5_B3.pdb, framework_VHH.pdb, framework_VHVL.pdb, etc.
- `pipeline_scripts/`    — every run script + inference code: run_campaign.sh, run_cd3_campaign.sh,
                           run_cea5_campaign.sh, run_proteinmpnn.sh, chain_*.sh, rfdiffusion_inference.py,
                           proteinmpnn_interface_design.py, rf2_predict.py, Dockerfile, rfantibody.def, README.md, CLAUDE.md
- `rfantibody_config/`   — RFdiffusion config + util
- `rfdiffusion_outputs/` — raw RFdiffusion outputs (dated run dirs + 8bw0 test set)
- `proteinmpnn_cloud/`   — ProteinMPNN cloud setup: backbone tarballs (bb_41BB/CD27/CD3), env, weights, source
- `model_weights/`       — RFdiffusion_Ab.pt, RF2_ab.pt, ProteinMPNN_v48_noise_0.2.pt (the weights actually used)
- `smoke_in/`, `smoke_out/` — pipeline smoke tests

### binder_design_deliverable/   (deposited AF3 campaign results)
- `af3_campaign/structures/` — 39 AF3-folded finalist CIFs incl. co-leads:
                               41bb_r0867_41bb_62, cd27_r0011_cd27_24, cd3_r0478_cd3_65, cea5_r0420_cea5_52
- `af3_campaign/sequences/`  — designed binder sequences (FASTA + panels)
- `af3_campaign/scores/`     — ipSAE, fold-scRMSD, AF3 confidence scores
- `af3_campaign/analysis/`   — funnel, retained-vs-filtered, CDR-property tables + figures
- `af3_campaign/docs/`       — AF3_TESTABLE_PANELS.md, AF3_FINAL_FINALISTS.md, methods, README
- `mpnn_sequences/`          — ProteinMPNN unique-CDR sequence tables (temp 0.2)
- `rfab_inputs/`, `rfab_outputs/`, `cea5/`, `templates/`

### animation_4arm/   (current deck timelapse: source + video)
- `S_*.pdb`               — the 4 oriented receptor structures (41bb trimer, cd27 full, cd3, cea5 full)
- `F_41bb_*.pdb`, `F_cd27/`, `F_cd3_*.pdb`, `F_cea5/` — 50 binder diffusion frames per target
- `frames_png/`           — 50 rendered PNG frames (noise→docked)
- `fourarm_diffusion_final.mp4` — final single-pass video (noise→docked)
- `fourarm_diffusion_loop.mp4`  — forward-only 3× loop for slides
- `membrane_long.pdb`, `domain_colors.json`, `cdr_windows.json`, `hotspots.json`, `make_frames.cxc` — scene/render config

### animation_2arm_prior/   (earlier 2-arm 41BB+CD3 scene, incl. 6MGP multimer assembly.pdb)

## Referenced, NOT copied (recoverable, avoids multi-GB duplication)
- RFantibody repo source (`src/`, `.venv/`, `RFantibody-main.zip`) — pip/git installable
- `scripts/examples/` (4.3 GB) — RFantibody repo DEMO data, not our campaign
- AlphaFold3 databases (627 GB) + AF3 weights — external; see af3_campaign/docs for the split-DB/cloud-GPU method
- Full AF3 per-fold raw outputs — the scored summary (scores/) + finalist CIFs (structures/) are the distilled record;
  full reproducibility tarball is also saved as artifact df9a7067-76dc-48eb-8722-845a8ab69f93

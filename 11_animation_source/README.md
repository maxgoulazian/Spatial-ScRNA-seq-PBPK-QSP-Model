# 4-Arm Binder Diffusion Timelapse — Complete Asset Folder

The RFdiffusion de novo binder-design timelapse: four designed binders (against 4-1BB, CD27,
CD3, CEACAM5) diffusing from noise onto their epitopes on a shared lipid membrane.

Location: /media/balthasar-lab/RAID4/costim_engager_counterscreen/TIMELAPSE_4ARM/

## Folders

- **01_videos/** — THE DELIVERABLES
  - `fourarm_diffusion_final.mp4`  — single pass, noise→docked (7.5 s)
  - `fourarm_diffusion_loop.mp4`   — forward-only 3× loop for slides (22.5 s)
  - `fourarm_docked_still.png`     — final docked key frame

- **02_rendered_frames/** — the 50 rendered PNG frames (f_000 = noise → f_049 = docked)
  that ffmpeg encodes into the videos.

- **03_receptor_structures/** — the 4 oriented receptor structures placed in the scene
  (S_41bb_receptor.pdb = 6MGP trimer, S_cd27_fullreceptor.pdb, S_cd3_receptor.pdb,
  S_cea5_fullreceptor.pdb) + ligand/multimer overlays.

- **04_binder_diffusion_frames/** — 50 binder poses per target (41bb, cd27, cd3, cea5),
  transformed into scene coordinates. These are the moving objects, one file per frame.

- **05_scene_config/** — membrane_long.pdb (seamless lipid bilayer), scene_membrane.pdb,
  domain_colors.json (per-protein by-domain gradients), cdr_windows.json (CDR loops → gold),
  hotspots.json / paratopes.json (epitopes), AF-P26842_CD27_full.pdb (full CD27 overlay).

- **06_render_scripts/** — ChimeraX .cxc scripts. `make_frames.cxc` is the final render
  recipe (static scene built once, 50 frames looped, balls + gold CDRs, supersample).

- **07_source_trajectories/** — the raw RFdiffusion denoising trajectories (*_Xt-1_traj.pdb,
  50 frames each) the diffusion motion was derived from, per target. This is the ground-truth
  physics: Rg goes 13→7 Å (noise→docked).

## Re-render
    cd <this folder>
    /usr/bin/chimerax --offscreen --nogui 06_render_scripts/make_frames.cxc   # writes 02_rendered_frames/
    ffmpeg -framerate 10 -i 02_rendered_frames/f_%03d.png \
      -vf "scale=1920:-2,tpad=stop_mode=clone:stop_duration=2.5" \
      -c:v libx264 -pix_fmt yuv420p -crf 17 01_videos/fourarm_diffusion_final.mp4
(paths in make_frames.cxc point at the anim4 working dir; adjust to this folder's structure if re-running standalone)

## Style (locked)
By-domain gradient receptor surfaces (41BB reds / CD27 oranges / CD3 greys / CEA5 blues, 30% transparent),
binder as balls (grey framework + bright-gold CDR loops), house-style lipid membrane, zoomed, binder-faces to camera.

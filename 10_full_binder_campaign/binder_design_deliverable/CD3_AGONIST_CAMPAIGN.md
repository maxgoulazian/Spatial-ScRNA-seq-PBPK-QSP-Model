# CD3 Agonist Arm — VH/VL de novo Design Campaign

## Why VH/VL (not VHH) for CD3
- Every clinically validated anti-CD3 (OKT3, SP34, UCHT1 families -> blinatumomab,
  the SP34-derived class, tebentafusp) is a conventional VH/VL. No validated
  anti-CD3 VHH in the clinic. The paratope is a two-domain interface.
- CD3 is the Signal-1 REDIRECTOR arm: its job is TCR engagement; affinity is a
  TUNABLE QSP knob (KD_CD3), not a maximize target (native pHLA:TCR is ~1000x
  weaker than these mAbs — over-affinity drives CRS).

## Agonism method — SAME approach as 4-1BB (literature-defined epitope), DIFFERENT mechanism
- 4-1BB agonism = clustering geometry (6MGP 3:3 assembly; urelumab CRD1 filter).
- CD3 agonism = CROSSLINKING-dependent: the TCE bridges CD3<->TAA to cluster TCRs.
  There is NO clustering-clash filter for CD3 — the agonist determinant is
  BINDING THE VALIDATED IMMUNODOMINANT CD3-epsilon EPITOPE (OKT3/UCHT1 site),
  which is what confers TCE T-cell activation.
- Structural analog of 6MGP/6MHR = PDB 1SY6: OKT3 Fab + CD3-gamma-epsilon
  heterodimer (Kjer-Nielsen et al., PNAS 2004, 2.1 A). This is the therapeutic-
  mAb-bound structure — the literature-defined agonism epitope.

## Epitope (verified: crystal contacts in 1SY6 + HDX in patent literature)
- OKT3 literature epitope on CD3-epsilon: residues 29-37, 79-84, 87-89 (mature
  numbering); conformational, requires gamma-epsilon association (OKT3 does not
  bind epsilon alone). UCHT1 overlaps this immunodominant site.
- In 1SY6 chain-A (single-chain gamma-epsilon) numbering, the crystallographic
  OKT3 footprint maps to residues 139-192 (epsilon portion).
- HOTSPOTS selected by #interface contacts (>=3, top 6):
  T186(R,9) T141(E,7) T189(K,7) T188(S,6) T155(E,4) T187(G,4).

## Campaign (identical engine/params to 4-1BB & CD27, VH/VL variant)
- Framework: hu-4D5-8_Fv (trastuzumab Fv, chains H+L) — the RFantibody VH/VL scaffold.
- design_loops = [L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]  (all 6 CDRs; VHH used H-only).
- target = CD3 gamma-epsilon heterodimer (chain T), hotspots above.
- RFdiffusion_Ab.pt, 100 backbones, deterministic. Then ProteinMPNN 8 seqs/backbone,
  then AF3/Boltz-2 iPTM screen (per AF3_SCREENING_PROTOCOL.md).
- Agonist filter = epitope-fidelity (does the design engage the OKT3/UCHT1 site) +
  affinity-in-tunable-range, NOT clustering-clash.

## Queue
- run_cd3_campaign.sh written; inputs staged (target_CD3.pdb, framework_VHVL.pdb).
- Launch AFTER 41BB+CD27 finish (sequential on the single A2000), OR concurrently
  (RFdiffusion uses 3.8/12 GB; 2 procs ~7.6 GB, fits — but time-slices, no speedup).

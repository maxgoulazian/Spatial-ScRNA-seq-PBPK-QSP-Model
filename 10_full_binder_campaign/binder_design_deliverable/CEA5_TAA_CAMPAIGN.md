# CEACAM5 (CEA) TAA Arm — VH/VL de novo Design Campaign

## Rationale
- CEACAM5 = top-2 balanced TAA finalist, highest tumor-restriction with CEACAM6
  (z_restriction 2.39), replicates across both CRC cohorts (Lee et al. 2020
  SMC GSE132465 + KUL3 GSE144735). The likely winning TAA arm.
- CLINICALLY PRECEDENTED (cibisatamab CEA-TCB, ADCs M9140/tusamitamab) -> serves
  as a de novo POSITIVE CONTROL: pipeline recovering a CEA binder against a target
  with clinical ground truth validates the whole binder-design workflow.
- UNLIKE the costim (4-1BB/CD27) and CD3 arms, this is a TAA anchor: NO agonist
  filter. The design goal is high-affinity membrane-proximal binding, full stop.
  Gate = plain binder gate (AF3 iPTM + Rosetta interaction), no clustering/epitope
  agonism constraint.

## Format: VH/VL (matches clinical anti-CEA format; framework hu-4D5-8_Fv trastuzumab Fv)
- design_loops = [L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13] (all 6 CDRs)

## Target epitope: A3-B3 membrane-proximal domain
- CEACAM5 (UniProt P06731) domain architecture (UniProt-verified):
  N/IgV 35-144 (distal) | A1-B1 145-315 | A2-B2 323-495 | A3 501-588 | B3 593-675 (proximal, nearest GPI).
- A3-B3 = residues 501-675. Membrane-PROXIMAL by design: a proximal epitope forces
  a tight immunological synapse -> maximal T-cell cytotoxicity (same principle as
  cibisatamab/CEA-TCB, which deliberately targets the proximal A3-B3 / IgC domains
  rather than the distal N domain).
- Structure: AlphaFold AF-P06731-F1 v6 (mean pLDDT 91 across A3-B3).
- HOTSPOTS: spatially-clustered exposed patch on B3 (18.5 A diameter, one paratope),
  side-chain-bearing residues only (no P/G anchors):
  T607(S) T630(N) T642(I) T644(K) T648(N) T675(S).

## Queue
- run_cea5_campaign.sh written; inputs staged (target_CEA5_A3B3.pdb, framework_VHVL.pdb).
- chain_cea5.sh watcher (PID logged) auto-launches CEA5 after CD3 finishes.
- Order: 41BB -> CD27 -> CD3 -> CEA5. Then ProteinMPNN (run_proteinmpnn.sh,
  8 seqs/backbone, all-6-loops for VH/VL arms) -> AF3/Boltz-2 iPTM screen.

## ProteinMPNN — smoke-tested OK
- 41BB backbone through proteinmpnn_interface_design.py: 8 seqs in 15 s, GPU,
  loops correctly detected. ~25 min per 100-backbone target.
- VHH arms (41BB,CD27): loop_string H1,H2,H3. VH/VL arms (CD3,CEA5): H1,H2,H3,L1,L2,L3.
- settings: temperature 1e-6, omit_AAs CX, num_connections 48 (Baker/Bennett defaults).

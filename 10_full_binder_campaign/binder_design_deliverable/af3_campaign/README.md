# De Novo AF3 Binder Campaign — Model Deposit
_Deposited 2026-07-13 10:43. Complete outputs of the RFdiffusion→ProteinMPNN→AlphaFold3 binder-design campaign for the costim counter-screen nominated targets._

## What this is
A genome-scale CD4 counter-screen nominated **4-1BB + CD27** as liability-clean costim co-leads. This directory holds the de novo antibody binders designed against those co-leads plus the CD3 redirector arm and the CEA5 tumor-antigen anchor — the buildable-molecule layer of the project.

## Campaign funnel
RFdiffusion backbones (100/target) → ProteinMPNN sequences (temp 0.2, ~1000 unique-CDR/target) → **AF3 screen 3,964 designs (1-seed)** → **150 refold (3-seed)** → **39 finalists (10-seed)**.
Ran on 3 cloud A100 pods. Screen completeness 3,964/3,974 = 99.75%.

## Scoring (3 orthogonal gates)
- **ipSAE** (MIN-of-directions, Dunbrack 10/15) — interface quality, primary ranker
- **fold-scRMSD** (binder-on-binder Kabsch, single-sqrt RMSD, <2Å = hard gate) — self-consistency, ChimeraX-validated
- **10-seed consistency** (ranking-score mean/std across 10 independent AF3 seeds) — annotated

## Co-lead binders (top LEAD per target)
| Target | Design | Format | ipSAE | fold-scRMSD |
|---|---|---|---|---|
| 4-1BB | 41bb_r0867_41bb_62 | VHH | 0.492 | 0.86Å |
| CD27 | cd27_r0011_cd27_24 | VHH | 0.552 | 0.60Å |
| CD3 | cd3_r0478_cd3_65 | VH/VL | 0.491 | 1.00Å |
| CEA5 | cea5_r0420_cea5_52 | VH/VL | 0.652 | 1.04Å |

Testable panel (LEAD+CONFIRM, fold<2Å) = **27 designs**: 4-1BB 2 / CD27 5 / CD3 8 / CEA5 12.

## Directory layout
```
af3_campaign/
├── sequences/                    ← BINDER SEQUENCES (the deliverable)
│   ├── finalists_39.fasta         all 39, rich headers (target/role/class/ipSAE/fold)
│   ├── panel_27_testable.fasta    LEAD+CONFIRM order-and-test set
│   ├── coleads.fasta              4 co-lead binders (1 per target)
│   ├── cdr_sequences.csv          full CDR + per-CDR H1/H2/H3 + charge/GRAVY
│   └── per_target/                {41bb,cd27,cd3,cea5}_finalists.fasta
├── structures/                   ← 39 AF3 best-model complexes (.cif; chain T=target, H/L=binder)
├── scores/                       ← screen_3964, refold, panel39, scrmsd, combined_verdict
├── analysis/                     ← T1-T4 tables, 2 figures, af3_analysis_universe.parquet (3964×48)
└── docs/                         ← testable-panels doc, state-recovery index, reproducibility.tgz, backbone map
```

## Sequence format notes
- **VHH (4-1BB, CD27):** single chain H (humanized camelid framework, ~116 aa). CDRs on H1/H2/H3.
- **VH/VL (CD3, CEA5):** two chains H (~117–121 aa) + L (~108 aa), trastuzumab-derived scaffold.
- Sequences extracted directly from the AF3 folded structures (authoritative, per-chain).
- Design naming: `{target}_r{mpnn_index}_{target}_{backbone}` — e.g. 41bb_r0867_41bb_62 = MPNN unique-CDR record 867, RFdiffusion backbone 62.

## Provenance / reproducibility
- `docs/af3_design_backbone_map.txt` — 3,974-row design→pod→backbone map
- `docs/af3_campaign_reproducibility.tgz` — all runner scripts, ProteinMPNN params, parse code, panel defs
- `docs/AF3_STATE_RECOVERY.md` — master artifact index (durable artifact-store IDs)
- Key artifact IDs: panel39_final 326c8d3f · analysis_universe a70c2829 · combined_verdict(v2) 03829d25

## Honest limits (for downstream users)
- **In-silico only** — no wet-lab binding yet. SPR/BLI KD + agonism cell assay = bench next step.
- **4-1BB panel is thin** (2 designs pass both gates) — hard VHH-vs-cysteine-rich-TNFR target. CEA5 is deep (12/12).
- Animation/epitope-overlay showcase design **41BB_23** is a DIFFERENT design from AF3 co-lead **r0867** — epitope recapitulation verified for 41BB_23, AF3 interface confidence for r0867.

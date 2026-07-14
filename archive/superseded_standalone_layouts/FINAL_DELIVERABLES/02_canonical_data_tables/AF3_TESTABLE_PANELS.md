# AF3 Binder Design — FINAL TESTABLE PANELS

_Generated 2026-07-13 10:19. 39 designs 10-seed-confirmed; per-target testable groups, not single leads._

## Selection logic
Each design carries three orthogonal metrics: **ipSAE**(MIN-of-directions, interface quality), **fold-scRMSD** (binder-on-binder Kabsch vs RFdiffusion design — does the ProteinMPNN sequence fold as designed; ChimeraX-validated), and **10-seed consistency** (ranking-score mean/std across 10 independent AF3 seeds).
- **fold-scRMSD < 2Å = hard gate** (a design whose sequence doesn't fold into its designed backbone is not real, regardless of interface score).
- **ipSAE grades within the fold-passing set**: LEAD (>0.35), CONFIRM (0.25–0.35), EXPLORATORY (<0.25).
- **Seed-stability annotated, not gated** — a hard std cutoff was found to flip leads (demoting the strongest binder over a hair of variance), so it is reported per-design instead.

## Recommended panel per target (LEAD + CONFIRM = order-and-test set)

| Target | Role | Format | Panel size | Lead | Lead ipSAE | Lead fold-scRMSD |
|---|---|---|---|---|---|---|
| 4-1BB (TNFRSF9) | costim co-lead | VHH | 2 | `41bb_r0867_41bb_62` | 0.492 | 0.86Å |
| CD27 (TNFRSF7) | costim co-lead | VHH | 5 | `cd27_r0011_cd27_24` | 0.552 | 0.60Å |
| CD3ε | redirector arm | VH/VL | 8 | `cd3_r0478_cd3_65` | 0.491 | 1.00Å |
| CEA5 (CEACAM5) | TAA tumor anchor | VH/VL | 12 | `cea5_r0420_cea5_52` | 0.652 | 1.04Å |

**Total: 27 designs across 4 targets** (2 / 5 / 8 / 12 for 4-1BB / CD27 / CD3 / CEA5).

## Full graded panels

### 4-1BB (TNFRSF9) — VHH, costim co-lead
| design | class | ipSAE | pDockQ | fold-scRMSD | 10-seed mean | seed std | stable |
|---|---|---|---|---|---|---|---|
| 41bb_r0025_41bb_14 | CONFIRM | 0.303 | 0.370 | 0.92 | 0.758 | 0.008 | ✓ |
| 41bb_r0658_41bb_59 | DROP_foldfail | 0.366 | 0.316 | 5.40 | 0.625 | 0.036 | ✓ |
| 41bb_r0022_41bb_66 | DROP_foldfail | 0.014 | 0.199 | 7.11 | 0.401 | 0.123 | · |
| 41bb_r0750_41bb_3 | EXPLORATORY | 0.224 | 0.300 | 1.02 | 0.592 | 0.160 | · |
| 41bb_r0492_41bb_3 | EXPLORATORY | 0.206 | 0.297 | 0.52 | 0.689 | 0.012 | ✓ |
| 41bb_r0147_41bb_94 | EXPLORATORY | 0.174 | 0.311 | 0.80 | 0.593 | 0.145 | · |
| 41bb_r0030_41bb_14 | EXPLORATORY | 0.118 | 0.312 | 0.88 | 0.483 | 0.172 | · |
| 41bb_r0867_41bb_62 | LEAD | 0.492 | 0.426 | 0.86 | 0.797 | 0.023 | ✓ |

### CD27 (TNFRSF7) — VHH, costim co-lead
| design | class | ipSAE | pDockQ | fold-scRMSD | 10-seed mean | seed std | stable |
|---|---|---|---|---|---|---|---|
| cd27_r0513_cd27_51 | CONFIRM | 0.347 | 0.312 | 0.77 | 0.735 | 0.013 | ✓ |
| cd27_r0240_cd27_80 | CONFIRM | 0.292 | 0.165 | 0.91 | 0.680 | 0.026 | ✓ |
| cd27_r0204_cd27_55 | CONFIRM | 0.281 | 0.302 | 0.72 | 0.648 | 0.036 | ✓ |
| cd27_r0254_cd27_54 | DROP_foldfail | 0.485 | 0.400 | 7.76 | 0.723 | 0.020 | ✓ |
| cd27_r0235_cd27_36 | DROP_foldfail | 0.255 | 0.491 | 3.70 | 0.635 | 0.016 | ✓ |
| cd27_r0011_cd27_24 | LEAD | 0.552 | 0.335 | 0.60 | 0.792 | 0.103 | · |
| cd27_r0846_cd27_13 | LEAD | 0.451 | 0.416 | 0.44 | 0.770 | 0.015 | ✓ |

### CD3ε — VH/VL, redirector arm
| design | class | ipSAE | pDockQ | fold-scRMSD | 10-seed mean | seed std | stable |
|---|---|---|---|---|---|---|---|
| cd3_r0664_cd3_65 | CONFIRM | 0.337 | 0.187 | 1.57 | 0.660 | 0.036 | ✓ |
| cd3_r0834_cd3_14 | CONFIRM | 0.334 | 0.165 | 1.15 | 0.708 | 0.007 | ✓ |
| cd3_r0618_cd3_30 | CONFIRM | 0.327 | 0.210 | 0.79 | 0.705 | 0.012 | ✓ |
| cd3_r0868_cd3_81 | CONFIRM | 0.301 | 0.214 | 1.12 | 0.734 | 0.014 | ✓ |
| cd3_r0069_cd3_52 | CONFIRM | 0.279 | 0.196 | 1.14 | 0.675 | 0.012 | ✓ |
| cd3_r0293_cd3_30 | EXPLORATORY | 0.243 | 0.212 | 0.90 | 0.689 | 0.017 | ✓ |
| cd3_r0973_cd3_29 | EXPLORATORY | 0.237 | 0.072 | 1.38 | 0.750 | 0.012 | ✓ |
| cd3_r0470_cd3_5 | EXPLORATORY | 0.232 | 0.167 | 1.04 | 0.668 | 0.012 | ✓ |
| cd3_r0147_cd3_30 | EXPLORATORY | 0.230 | 0.207 | 0.87 | 0.661 | 0.021 | ✓ |
| cd3_r0478_cd3_65 | LEAD | 0.491 | 0.180 | 1.00 | 0.724 | 0.136 | · |
| cd3_r0925_cd3_65 | LEAD | 0.426 | 0.154 | 0.95 | 0.726 | 0.021 | ✓ |
| cd3_r0303_cd3_73 | LEAD | 0.404 | 0.312 | 0.95 | 0.711 | 0.070 | ✓ |

### CEA5 (CEACAM5) — VH/VL, TAA tumor anchor
| design | class | ipSAE | pDockQ | fold-scRMSD | 10-seed mean | seed std | stable |
|---|---|---|---|---|---|---|---|
| cea5_r0420_cea5_52 | LEAD | 0.652 | 0.215 | 1.04 | 0.839 | 0.114 | · |
| cea5_r0075_cea5_10 | LEAD | 0.587 | 0.300 | 1.03 | 0.880 | 0.012 | ✓ |
| cea5_r0543_cea5_67 | LEAD | 0.552 | 0.313 | 0.96 | 0.847 | 0.073 | ✓ |
| cea5_r0160_cea5_93 | LEAD | 0.544 | 0.281 | 1.20 | 0.827 | 0.096 | ✓ |
| cea5_r0538_cea5_16 | LEAD | 0.538 | 0.276 | 0.93 | 0.742 | 0.116 | · |
| cea5_r0473_cea5_65 | LEAD | 0.527 | 0.241 | 1.12 | 0.865 | 0.009 | ✓ |
| cea5_r0753_cea5_37 | LEAD | 0.514 | 0.270 | 1.09 | 0.865 | 0.008 | ✓ |
| cea5_r0094_cea5_61 | LEAD | 0.483 | 0.350 | 1.17 | 0.847 | 0.015 | ✓ |
| cea5_r0269_cea5_43 | LEAD | 0.467 | 0.268 | 1.11 | 0.849 | 0.031 | ✓ |
| cea5_r0179_cea5_72 | LEAD | 0.462 | 0.244 | 0.81 | 0.702 | 0.210 | · |
| cea5_r0231_cea5_35 | LEAD | 0.458 | 0.281 | 1.14 | 0.802 | 0.013 | ✓ |
| cea5_r0299_cea5_61 | LEAD | 0.453 | 0.318 | 1.32 | 0.832 | 0.016 | ✓ |

## Key reads
- **4-1BB (hard target)**: thinnest bench — only 23 designs in the entire 3,964 screen cleared btiPTM>0.45, and just 2 pass fold+interface (`r0867` ipSAE 0.49 fold 0.86Å is the lead; `r0025` fold 0.92Å is the seed-stable backup, std 0.008). Expected for de novo VHH against a small cysteine-rich TNFR domain.
- **CD27**: 5-design panel; `r0011` strongest (ipSAE 0.553) though seed-variable; `r0846` (ipSAE 0.451, fold 0.44Å, seed-stable) is the safest single pick.
- **CD3 (redirector)**: 8-design panel; deep, well-behaved — `r0478` top ipSAE, several sub-0.02-std alternatives.
- **CEA5 (TAA)**: 12/12 all LEAD-class — easiest target; `r0420` highest ipSAE (0.652) but `r0075` (0.587, std 0.012) is the cleaner reproducible pick.
- **4 costim fold-failures dropped** (2 per co-lead): 41bb_r0658 (5.4Å), 41bb_r0022 (7.1Å), cd27_r0254 (7.8Å), cd27_r0235 (3.7Å) — high ipSAE but sequence doesn't fold into designed backbone; correctly deprioritized by the fold gate.

## Data
- `panel39_final.csv` — all 39 designs, all metrics, panel_class
- `AF3_STATE_RECOVERY.md` — full artifact index
- Raw structures: pod{1,2,3}_tenfull.tgz (39 designs, full per-seed .cif)
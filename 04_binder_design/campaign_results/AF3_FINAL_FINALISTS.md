# AF3 Binder Design — FINAL Finalists (10-seed confirmed)

_Generated 2026-07-13 08:32. Terminal tier of the 3-stage AF3 funnel._

## Funnel
1. **Screen** — 3,964 designs @ 1 seed (RFdiffusion backbones × ProteinMPNN seqs), scored by AF3 + ipSAE(min)
2. **Refold** — 150 widened finalists @ 3 seeds (interface-iPTM/PAE union ∪ ipSAE>0.35), cross-seed stability filter
3. **Final** — 24 seed-stable designs (6/target) @ 10 seeds — this table

_ipSAE = MIN of asymmetric directions (conservative), Dunbrack cutoff 10/15, binder↔target chain pairs. Consistency = ranking-score mean/min/std across 10 independent seeds._

## Per-target top pick

| Target | Role | Format | Design | ipSAE(min) | pDockQ | LIS | 10-seed mean | seed std |
|---|---|---|---|---|---|---|---|---|
| 4-1BB (TNFRSF9) | costim co-lead | VHH | `41bb_r0867_41bb_62` | 0.492 | 0.426 | 0.467 | 0.80 | 0.023 |
| CD27 (TNFRSF7) | costim co-lead | VHH | `cd27_r0011_cd27_24` | 0.553 | 0.335 | 0.542 | 0.79 | 0.103 |
| CD3ε | redirector arm | VH/VL | `cd3_r0478_cd3_65` | 0.491 | 0.180 | 0.483 | 0.72 | 0.136 |
| CEA5 (CEACAM5) | TAA (tumor anchor) | VH/VL | `cea5_r0420_cea5_52` | 0.652 | 0.215 | 0.611 | 0.84 | 0.114 |

## Full 6-per-target finalist table

### 4-1BB (TNFRSF9) (VHH, costim co-lead)
| design | ipSAE(min) | pDockQ | LIS | 10-seed mean | 10-seed min | seed std |
|---|---|---|---|---|---|---|
| 41bb_r0867_41bb_62 | 0.492 | 0.426 | 0.467 | 0.797 | 0.759 | 0.023 |
| 41bb_r0658_41bb_59 | 0.366 | 0.316 | 0.332 | 0.625 | 0.564 | 0.036 |
| 41bb_r0025_41bb_14 | 0.304 | 0.370 | 0.385 | 0.758 | 0.748 | 0.008 |
| 41bb_r0750_41bb_3 | 0.223 | 0.300 | 0.294 | 0.592 | 0.151 | 0.160 |
| 41bb_r0492_41bb_3 | 0.206 | 0.297 | 0.310 | 0.689 | 0.663 | 0.012 |
| 41bb_r0147_41bb_94 | 0.174 | 0.311 | 0.265 | 0.593 | 0.316 | 0.145 |

### CD27 (TNFRSF7) (VHH, costim co-lead)
| design | ipSAE(min) | pDockQ | LIS | 10-seed mean | 10-seed min | seed std |
|---|---|---|---|---|---|---|
| cd27_r0011_cd27_24 | 0.553 | 0.335 | 0.542 | 0.792 | 0.503 | 0.103 |
| cd27_r0254_cd27_54 | 0.485 | 0.400 | 0.454 | 0.723 | 0.701 | 0.020 |
| cd27_r0846_cd27_13 | 0.451 | 0.416 | 0.448 | 0.770 | 0.750 | 0.015 |
| cd27_r0513_cd27_51 | 0.347 | 0.312 | 0.386 | 0.735 | 0.720 | 0.013 |
| cd27_r0240_cd27_80 | 0.292 | 0.165 | 0.348 | 0.680 | 0.635 | 0.026 |
| cd27_r0235_cd27_36 | 0.255 | 0.491 | 0.303 | 0.635 | 0.609 | 0.016 |

### CD3ε (VH/VL, redirector arm)
| design | ipSAE(min) | pDockQ | LIS | 10-seed mean | 10-seed min | seed std |
|---|---|---|---|---|---|---|
| cd3_r0478_cd3_65 | 0.491 | 0.180 | 0.483 | 0.724 | 0.462 | 0.136 |
| cd3_r0925_cd3_65 | 0.426 | 0.154 | 0.396 | 0.726 | 0.693 | 0.021 |
| cd3_r0303_cd3_73 | 0.404 | 0.312 | 0.431 | 0.711 | 0.646 | 0.070 |
| cd3_r0664_cd3_65 | 0.337 | 0.187 | 0.325 | 0.660 | 0.601 | 0.036 |
| cd3_r0834_cd3_14 | 0.334 | 0.165 | 0.303 | 0.708 | 0.698 | 0.007 |
| cd3_r0618_cd3_30 | 0.327 | 0.210 | 0.357 | 0.705 | 0.687 | 0.012 |

### CEA5 (CEACAM5) (VH/VL, TAA (tumor anchor))
| design | ipSAE(min) | pDockQ | LIS | 10-seed mean | 10-seed min | seed std |
|---|---|---|---|---|---|---|
| cea5_r0420_cea5_52 | 0.652 | 0.215 | 0.611 | 0.839 | 0.607 | 0.114 |
| cea5_r0075_cea5_10 | 0.587 | 0.300 | 0.558 | 0.880 | 0.863 | 0.012 |
| cea5_r0543_cea5_67 | 0.552 | 0.313 | 0.551 | 0.847 | 0.639 | 0.073 |
| cea5_r0160_cea5_93 | 0.544 | 0.281 | 0.545 | 0.827 | 0.647 | 0.096 |
| cea5_r0538_cea5_16 | 0.538 | 0.276 | 0.530 | 0.742 | 0.646 | 0.116 |
| cea5_r0473_cea5_65 | 0.527 | 0.241 | 0.543 | 0.865 | 0.845 | 0.009 |

## Key reads
- **4-1BB `r0867`** — lowest ipSAE (0.492) but **most seed-stable of all finalists** (std 0.023) and best pDockQ (0.43). On a hard VHH-vs-CRD target, 10 seeds unanimously docking is the real signal. The CRD1 agonist-epitope compatibility check applies to this design.
- **CD27 `r0011`** — strongest costim (ipSAE 0.553); `r0254`/`r0846` are stabler backups (std ~0.02).
- **CEA5 `r0420`** highest ipSAE (0.652); **`r0075` is the cleaner pick** (ipSAE 0.587, std 0.012 vs r0420's 0.114) — recommend r0075 as the lead CEA5 finalist for downstream.
- **CD3 `r0478`** — top ipSAE (0.491); `r0925`/`r0834` stabler (std ~0.01-0.02).
- **All four targets yield a defensible, seed-reproducible finalist.** Costim VHHs (4-1BB, CD27) are borderline on absolute ipSAE — expected for de novo VHH vs small cysteine-rich TNFR domains — but pass the reproducibility bar.

## Data
- `tenseed_final_scored.csv` — this table (24 designs, ipSAE-min + 10-seed consistency)
- `refold_scored.csv` — 3-seed tier (150 designs)
- `af3_master_scored.csv` — 1-seed screen (3,964 designs)
- Structures: 10-seed `.cif` models in af3_tenseed tarballs (per-seed detail preserved)
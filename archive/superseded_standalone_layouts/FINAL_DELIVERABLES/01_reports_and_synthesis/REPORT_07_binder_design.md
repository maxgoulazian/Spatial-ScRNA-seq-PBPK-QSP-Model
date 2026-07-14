# 7. De novo binder design against the nominated targets

The three-axis screen nominates **4-1BB (TNFRSF9)** and **CD27 (TNFRSF7)** as the liability-clean costimulatory arms. This section converts that nomination into buildable molecules by generating de novo antibody binders against the two nominated costim receptors and against the two remaining arms a tumor-conditional costimulatory engager requires — the **CD3ε** redirector and the **CEACAM5 (CEA5)** tumor anchor — so that the receptor-selection result is delivered as a bench-ready starting panel rather than a target list. All results in this section are in silico; they are prioritization metrics that define an ordered wet-lab test set, not measured affinities or activities (Section 7.7).

## 7.1 Design pipeline

Binders were generated with a three-stage generative funnel. Backbones were sampled with **RFdiffusion** [1], sequences were assigned to each backbone with **ProteinMPNN** at sampling temperature 0.2 [2], and every binder–target complex was refolded and scored with **AlphaFold3** [3]. The costimulatory targets (4-1BB, CD27) were built as single-domain **VHH** binders and the redirector/anchor arms (CD3ε, CEACAM5) as paired **VH/VL** fragments, matching the formats used for each arm in a bispecific redirector. This RFdiffusion → ProteinMPNN → structure-prediction-filter workflow is the de novo antibody design protocol validated experimentally for single-domain and scFv binders against defined epitopes [4].

For the two costim receptors, diffusion was seeded on the ligand-competitive cysteine-rich-domain (CRD) hotspots so that designs are biased toward agonism-relevant epitopes: the 4-1BB campaign templated on the 4-1BB/4-1BBL complex (PDB 6MGP) and the CD27 campaign on the CD27/CD70 complex (PDB 7KX0).

Each design carries three orthogonal quality metrics:

1. **ipSAE** — the interface predicted Score from Aligned Errors, taken as the minimum of the two asymmetric chain directions (the conservative choice) with the Dunbrack 10/15 cutoff. ipSAE restricts the interface confidence estimate to residue pairs with low predicted aligned error, and separates true from false complexes more reliably than AlphaFold's global ipTM [5].
2. **fold-scRMSD** — the self-consistency check: the binder from the AlphaFold3 refold is superimposed (Kabsch) on the RFdiffusion design backbone, testing whether the ProteinMPNN sequence actually folds into the backbone it was designed for (ChimeraX-validated).
3. **10-seed consistency** — mean, minimum and standard deviation of the AlphaFold3 ranking score across ten independent seeds.

pDockQ [6] and LIS are reported alongside ipSAE as secondary interface descriptors.

## 7.2 Funnel attrition

The screen began with **3,964 designs** scored at a single seed (981 / 986 / 997 / 1,000 for 4-1BB / CD27 / CD3 / CEA5). A binder–target interface-iPTM gate at 0.45 (unioned with ipSAE > 0.35) retained **152 designs** for a 3-seed refold, and cross-seed-stable survivors were promoted to a **39-design** 10-seed finalist tier (**Table 7.1**). Overall pass rates were 3.8% to refold and 0.98% to the finalist tier.

**Table 7.1 — Funnel attrition per target** (`T1_funnel_attrition.csv`)

| Target | Screened (1-seed) | Refold (3-seed) | Finalist (10-seed) | % to refold | % to finalist |
|---|---|---|---|---|---|
| 4-1BB | 981 | 23 | 8 | 2.3 | 0.82 |
| CD27 | 986 | 51 | 7 | 5.2 | 0.71 |
| CD3ε | 997 | 36 | 12 | 3.6 | 1.20 |
| CEA5 | 1,000 | 42 | 12 | 4.2 | 1.20 |
| **All** | **3,964** | **152** | **39** | **3.8** | **0.98** |

The interface gate is what separates retained from filtered designs, and it does so cleanly: binder–target interface iPTM is bimodal, with the retained set (n = 152) sitting entirely above the 0.45 gate and the filtered set (n = 3,812) below it (**Figure 7.1b**). 4-1BB is the hardest target — only 23 of its 981 designs cleared the interface gate, expected for de novo VHH binding against a small, disulfide-constrained TNFR CRD.

![Figure 7.1 — AF3 screening funnel and interface-based separation. (a) Per-target attrition across the 1-seed screen, 3-seed refold and 10-seed finalist tiers (log scale). (b) Binder–target interface iPTM for retained (n=152) vs filtered-out (n=3,812) designs; the retained set separates cleanly above the 0.45 refold gate.]({{artifact:b0472532-f11c-4725-a646-b975051adf6d}})

## 7.3 Fold gate and interface grading

Because a sequence that does not fold into its designed backbone is not a real binder regardless of its interface score, **fold-scRMSD < 2 Å is applied as a hard gate**; ipSAE then grades the fold-passing set into **LEAD** (ipSAE > 0.35), **CONFIRM** (0.25–0.35) and **EXPLORATORY** (< 0.25). Seed stability is annotated per design rather than gated, because a hard standard-deviation cutoff was found to demote the strongest binders over marginal variance.

Of the 39 finalists, **35 pass the fold gate and 4 fail it**. The four fold-failures are informative: each carries an interface score that would otherwise place it in the panel but a scRMSD that shows the sequence does not adopt the designed fold — `cd27_r0254_cd27_54` (ipSAE 0.485, the second-highest CD27 interface score, scRMSD 7.8 Å), `41bb_r0658_41bb_59` (ipSAE 0.366, scRMSD 5.4 Å), `cd27_r0235_cd27_36` (ipSAE 0.255, scRMSD 3.7 Å) and `41bb_r0022_41bb_66` (scRMSD 7.1 Å). Grading on interface score alone would have advanced `cd27_r0254` as a co-lead; the fold gate correctly removes it. This is the single most consequential filter in the funnel and the reason absolute ipSAE is not used as a standalone ranking.

## 7.4 Testable panel

The order-and-test set is the union of LEAD and CONFIRM designs: **27 designs across the four targets** (2 / 5 / 8 / 12 for 4-1BB / CD27 / CD3 / CEA5). Per-target co-leads and their metrics are given in **Table 7.2**.

**Table 7.2 — Per-target co-leads** (`AF3_TESTABLE_PANELS.md`, `T3_finalist_headtohead.csv`, `af3_master_scored.csv`)

| Target | Role · format | Panel (LEAD+CONFIRM) | Co-lead design | ipSAE | fold-scRMSD (Å) | 10-seed mean | seed SD |
|---|---|---|---|---|---|---|---|
| 4-1BB | costim · VHH | 2 | `41bb_r0867_41bb_62` (LEAD) | 0.49 | 0.86 | 0.80 | 0.023 |
| | | | `41bb_r0025_41bb_14` (CONFIRM, stable backup) | 0.30 | 0.92 | 0.76 | 0.008 |
| CD27 | costim · VHH | 5 | `cd27_r0011_cd27_24` (LEAD) | 0.55 | 0.60 | 0.79 | 0.103 |
| | | | `cd27_r0846_cd27_13` (LEAD, seed-stable pick) | 0.45 | 0.44 | 0.77 | 0.015 |
| CD3ε | redirector · VH/VL | 8 | `cd3_r0478_cd3_65` (LEAD) | 0.49 | 1.00 | 0.72 | 0.136 |
| CEA5 | tumor anchor · VH/VL | 12 | `cea5_r0420_cea5_52` (LEAD) | 0.65 | 1.04 | 0.84 | 0.114 |
| | | | `cea5_r0075_cea5_10` (LEAD, reproducible pick) | 0.59 | 1.03 | 0.88 | 0.012 |

The costim VHHs sit at lower absolute ipSAE than the VH/VL arms — expected for single-domain binders against small cysteine-rich TNFR domains — but pass the fold and reproducibility bars. For the nominated 4-1BB arm the bench is thin by design: the fold-passing interface set contains only two VHHs, `r0867` (the lead: highest interface score, best pDockQ 0.43, and the most seed-stable finalist at SD 0.023) and `r0025` (the seed-stable backup, SD 0.008). CD27 offers a deeper five-design panel; where `r0011` is the strongest interface but seed-variable, `r0846` (ipSAE 0.45, fold-scRMSD 0.44 Å, SD 0.015) is the safer reproducible pick. CD3 gives a well-behaved eight-design panel and CEA5 is the easiest target, returning twelve LEAD-class designs.

## 7.5 Selection is on interface metrics, not CDR sequence bias

To confirm that finalist selection reflects predicted interface quality rather than an incidental sequence signature, retained and filtered designs were compared across interface metrics, CDR physicochemical properties and CDR sequence-liability motif counts (**`T2_retained_vs_filtered.csv`**). The five interface descriptors separate the two groups with overwhelming significance (binder–target chain iPTM Mann–Whitney p = 1×10⁻⁹⁷, Δmedian +0.49; interface PAE p = 2×10⁻⁹⁴, Δmedian −11.0 Å; screen iPTM p = 3×10⁻⁶³; screen ipSAE p = 2×10⁻²³; screen pDockQ p = 5×10⁻⁹). In contrast, **none of the CDR sequence properties discriminate** — CDR length (p = 0.18), net charge (p = 0.34), GRAVY hydrophobicity (p = 0.19), aromatic fraction (p = 0.54) and ProteinMPNN score (p = 0.78) are all non-significant, as are all four developability-liability motif counts (deamidation, isomerization, N-glycosylation, oxidation; p = 0.20–0.98). The funnel is therefore selecting on the structural interface, not filtering to a particular paratope chemistry.

## 7.6 CDR property profiles and developability

Across all four targets the fold-passing paratopes are net-negative at pH 7.4 (mean CDR net charge −1.1 to −2.1) with negative GRAVY (−0.32 to −0.72), i.e. hydrophilic and acidic — a low-nonspecificity profile (**Figure 7.2a–b**; `T4_pertarget_cdr_profile.csv`). The only fold-failure with a net-positive paratope (`41bb_r0022`, net charge +1) is removed by the fold gate, consistent with positive-charge paratopes being a known nonspecificity risk. CDR-H3 lengths are moderate (mean 8.7–10.2 residues; **Figure 7.2c**). The finalist panel's sequence-liability motif burden is comparable to that of the full 3,964-design screen and is not elevated (**Figure 7.2d**): relative to the screened population the 39 finalists carry a marginally higher isomerization-motif count (DG/S/T mean 0.67 vs 0.62) and a marginally lower N-glycosylation-motif count (NxS/T mean 0.21 vs 0.24; `T3_finalist_headtohead.csv`). This is a distinct comparison from the interface-gated split of Section 7.5: in that retained (n=152) vs filtered (n=3,812) comparison none of the four liability-motif counts differs significantly (deamidation p=0.98, isomerization p=0.74, N-glycosylation p=0.50, oxidation p=0.20; `T2_retained_vs_filtered.csv`), consistent with selection acting on the interface rather than on CDR sequence chemistry. Any residual motifs are addressable by point substitution during optimization.

![Figure 7.2 — CDR property profiles of the finalist panel. (a) Fold-passing paratopes are net-negative at pH 7.4 across all targets; the one net-positive design is a dropped fold-failure. (b) Finalists (colored) sit within the charge/hydrophobicity envelope of the full screen (grey). (c) CDR-H3 length distributions. (d) CDR sequence-liability motif burden, finalists vs all screened — not elevated.]({{artifact:38627667-8533-40b5-a44a-f2ba45b8e71f}})

## 7.7 The 4-1BB binder targets the non-blocking, agonism-permissive epitope

The 4-1BB arm carries a specific mechanistic requirement that distinguishes a costimulatory agonist from a ligand-blocking antagonist. The two clinical anti-4-1BB antibodies bind in structurally distinct modes: **urelumab** binds the membrane-distal CRD1, away from the 4-1BBL binding site (non-ligand-blocking), whereas **utomilumab** binds at the CRD3/CRD4 junction and competes with the natural ligand [7]. A costimulatory arm intended to amplify — not block — 4-1BB signaling should recapitulate the urelumab-class CRD1 epitope.

Mapping the designed 4-1BB VHH footprints onto the receptor confirms this. The diffusion-seeded CRD1 design shares **5 of its 7** interface residues with the urelumab non-blocking epitope and **zero** with the utomilumab ligand-blocking epitope (`41bb_23_clustering_compat.json`); a companion CRD1 design reaches **7 of 8** urelumab residues (Jaccard 0.44), again with zero utomilumab overlap (`41bb_epitope_comparison.json`). Both designs sit on the ligand-independent agonist face, not the blocking face. This is the desired configuration for a costim agonist arm and is consistent with the mechanistic rationale connecting CRD1 (urelumab-class) epitopes to receptor-clustering agonism [8]. It also reframes the known urelumab liability: urelumab's dose-limiting hepatotoxicity is an Fc/FcγR-crosslinking and systemic-agonism problem [8], which a tumor-conditional bispecific format addresses at the delivery layer while retaining the favorable non-blocking epitope geometry — orthogonal to, and compatible with, the receptor-selection result of this report.

## 7.8 Limitations and wet-lab next step

These are computational designs prioritized by structure-prediction confidence, and structure-prediction interface scores are ranking metrics, not affinities; the de novo antibody literature is explicit that current in-silico filters are weak binder/non-binder discriminators and that experimental screening remains necessary [4]. The 27-design LEAD+CONFIRM panel is accordingly framed as an ordered test set. The recommended validation sequence is: (i) express the panel as VHH (4-1BB, CD27) and scFv/Fab (CD3, CEA5); (ii) measure binding by SPR/BLI, resolving apparent affinity and, for the costim arms, epitope binning against urelumab (non-blocking) and utomilumab (blocking) reference antibodies to confirm the predicted CRD1 footprint; (iii) test agonism functionally (e.g. NF-κB reporter and primary-T-cell costimulation) for the 4-1BB and CD27 arms, and redirected-lysis activity for the assembled CD3×CEA5 bispecific; and (iv) affinity-mature the seed-variable leads. The fold-scRMSD gate is a self-consistency proxy for expressibility and should be confirmed by expression yield. Designs, per-seed structures and all scoring tables are deposited (`af3_master_scored.csv`, `T1`–`T4` tables, and the per-target `.cif` models) for reproducibility.

---

### References

1. Watson JL, Juergens D, Bennett NR, et al. De novo design of protein structure and function with RFdiffusion. *Nature* 620, 1089–1100 (2023). doi:10.1038/s41586-023-06415-8
2. Dauparas J, Anishchenko I, Bennett N, et al. Robust deep learning–based protein sequence design using ProteinMPNN. *Science* 378, 49–56 (2022). doi:10.1126/science.add2187
3. Abramson J, Adler J, Dunger J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630, 493–500 (2024). doi:10.1038/s41586-024-07487-w
4. Bennett NR, Watson JL, Ragotte RJ, et al. Atomically accurate de novo design of antibodies with RFdiffusion. *Nature* 649, 183–193 (2026); epub 5 Nov 2025; bioRxiv 2024.03.14.585103. doi:10.1038/s41586-025-09721-5
5. Dunbrack RL. Rēs ipSAE loquuntur: what's wrong with AlphaFold's ipTM score and how to fix it. *bioRxiv* 2025.02.10.637595 (2025). doi:10.1101/2025.02.10.637595
6. Bryant P, Pozzati G, Elofsson A. Improved prediction of protein–protein interactions using AlphaFold2. *Nat Commun* 13, 1265 (2022). doi:10.1038/s41467-022-28865-w
7. Chin SM, Kimberlin CR, Roe-Zurz Z, et al. Structure of the 4-1BB/4-1BBL complex and distinct binding and functional properties of utomilumab and urelumab. *Nat Commun* 9, 4679 (2018). doi:10.1038/s41467-018-07136-7
8. Qi X, Li F, Wu Y, et al. Optimization of 4-1BB antibody for cancer immunotherapy by balancing agonistic strength with FcγR affinity. *Nat Commun* 10, 2141 (2019). doi:10.1038/s41467-019-10088-1

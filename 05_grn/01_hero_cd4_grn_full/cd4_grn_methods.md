# CD4 CRISPRi Perturb-seq — genome-scale gene-regulatory network (GRN)

**Tier label:** *genome-scale GRNBoost2 co-expression GRN, validated against CRISPRi perturbations.*

This is the real Tier-2 GRN. It is **not** the prior lane's 12-node Modular-Response-Analysis
(MRA): that computed a 12×12 local interaction matrix among the costim receptors only. Here we
infer a genome-scale weighted directed regulatory network over **1,598 regulators → 18,127
target genes** (1,653,594 edges) with GRNBoost2, the gradient-boosting core of the SCENIC
algorithm, then validate inferred edges against the independent CRISPRi differential-expression
matrix and propagate each costim receptor's signal onto the effector / CRS / suppression / help
programs.

## 1. Substrate
- **Source:** `GWCD4i.pseudobulk_merged.h5ad`, from the provided dataset directory
  `hero_cd4_perturbseq_zhu2025` — the Marson/Pritchard genome-scale CD4⁺ T-cell CRISPRi Perturb-seq
  (the `zhu2025` directory tag is the dataset provenance as provided; the underlying publication /
  accession is not independently verified here). 278,684 pseudobulk profiles
  (guide×donor×condition) × 18,129 genes,
  raw summed UMI counts (verified integer-valued; library sizes 6.0k–5.2M).
- **State:** we use the **Stim48hr (sustained-activation)** condition — the activated state in which
  costim wiring and the CRS/suppressive programs are engaged — restricted to `keep_for_DE=True`
  QC-pass profiles (78,280 profiles).
- **Normalization:** CP10K + log1p per profile.
- **Collapse:** n_cells-weighted mean per **perturbed_gene** → one expression profile per CRISPRi
  target (11,332 perturbation profiles × 18,127 genes after dropping the non-ENSG spike-in feature
  and 1 zero-variance gene). This is the GRNBoost2 substrate: the perturbation-driven variance across
  knockdowns is exactly the signal GRN inference exploits. Genes are handled in Ensembl-ID space
  (unambiguous); symbols are mapped only for reporting/validation.

## 2. Regulators
- **Canonical TFs:** the human TF catalogue supplied in the project reference directory as
  `Lambert_2018_HumanTF.csv` (filtered to is_TF=Yes; the filename attributes it to the Lambert
  et al. 2018 human-TF compendium), intersected with the matrix → 1,569 TFs.
- **UNION** the 12 costim receptors (TNFRSF9, CD27, TNFRSF4, TNFRSF18, TNFRSF25, TNFRSF14, TNFRSF8,
  CD28, ICOS, CD226, CD2, CD40) and the program genes (CRS ∪ SUPP ∪ HELP).
- **Total regulators = 1,598.** Targets = all 18,127 genes (genome-scale).

## 3. GRNBoost2 (the GRN)
- **Algorithm:** `arboreto` GRNBoost2 — per-target stochastic gradient-boosting regression
  (`SGBM_KWARGS`: learning_rate 0.01, n_estimators 5000, max_features 0.1, subsample 0.9) with
  early stopping (window 25), feature importances = directed regulator→target edge weights. Identical
  regressor/importance extraction to the reference `grnboost2()` entrypoint.
- **Driver:** because arboreto 0.1.6's Dask-graph assembly is incompatible with the installed Dask
  2026.7 query-planning backend (`from_delayed` "Must supply at least one delayed object"), we drive
  the **same** `arboreto.core.infer_partial_network` regressor per target across a 20-process pool
  (capped at 20 workers — a second GRN job shares this 80-core box). Seed 777. Targets ordered
  regulators+program-genes first so the costim→program subnetwork is complete even under a timebox.
- **Signing:** GRNBoost2 importances are unsigned. Each edge is signed by the Pearson correlation of
  regulator and target across the 11,332 profiles (the CellOracle convention), giving
  `signed_importance` for directional propagation.

## 4. Costim → program DRIVE (network propagation)
For each costim receptor we compute its downstream drive onto four gene programs by signed
network propagation from the receptor node through a bounded multi-hop kernel
K = I + α·T + α²·T² (α=0.5) over the regulator subspace, then a final regulator→gene hop
(row-normalized signed weights). Program gene sets:
- **effector** = IFNG, GZMB, PRF1, TNF, IL2, TBX21, EOMES, NKG7, FASLG, KLRG1 (CD8-effector readout)
- **CRS** = TNF, IL2, IFNG
- **SUPP** = CTLA4, EBI3, ENTPD1, FOXP3, HAVCR2, IKZF2, IL10, IL12A, IL2RA, LAG3, LRRC32, NT5E, TGFB1, TIGIT, TNFRSF18
- **HELP** = BCL6, MAF, TOX2, CXCR5, TCF7

## 5. LOF → GOF (agonism) sign flip
The CRISPRi screen is **loss-of-function** (knockdown); an engager costim arm is **gain-of-function**
(agonism). We apply the sign flip explicitly: **agonism = −1 × knockdown-direction.** Operationally,
the network-propagation drive is computed as a **+Δ (agonism) shift** on the costim node — raising the
receptor and propagating the signed response onto each program. The knockdown (LOF) direction is the
negative of this. All reported drive/z values are therefore on the **agonism axis** (positive drive
= agonism raises that program), which is the axis the engager arm acts on and the axis the QSP model
consumes.

## 6. CRISPRi validation (Perturb-seq gold standard)
Inferred edges are validated against each receptor's **own CRISPRi knockdown signature** in the
independent DE matrix (`GWCD4i.DE_stats.h5ad`, `layers['log_fc']`, Stim48hr). Three complementary
metrics per receptor:
1. **Magnitude** — fraction of DE-testable out-edges whose target is significantly moved by the
   receptor's knockdown (adj_p < 0.1).
2. **Enrichment** — are |log_fc| of edge-targets systematically larger than non-edge genes?
   (one-sided Mann–Whitney; rank-biserial AUC effect). Continuous, so receptors with narrow KD
   signatures remain testable.
3. **Direction** — does edge sign predict KD direction? A positive (activating) co-expression edge
   predicts that knockdown lowers the target; concordance = sign(log_fc)==−edge_sign among moved
   targets (one-sided binomial).
Plus two aggregate cross-checks on the whole network: **importance→magnitude** (does a stronger
GRNBoost2 weight predict a larger KD effect, per broad-KD TF) and **top-25 enrichment** (the
enrichment metric restricted to each regulator's 25 strongest edges — the regime a thresholded
GRNBoost2 network is actually used in).

**Genome-scale validation results (full 1,653,594-edge network):**
- **Top-25 edge enrichment:** mean AUC **0.583** across 420 CRISPRi-targeted regulators, **69.8 %**
  with AUC > 0.5 (p = 1.7×10⁻¹⁶). A regulator's strongest GRN targets are systematically more
  perturbed by its own knockdown than random genes.
- **Importance → KD magnitude:** mean within-TF Pearson r = **+0.097**, positive for **78.7 %** of
  broad-KD regulators (p = 7.6×10⁻¹⁰). Higher edge weight ⇒ larger knockdown effect.
- **Edge sign → KD direction:** **61.7 %** concordance among moved edges (p = 3.3×10⁻¹⁹⁸).
- **All-edge enrichment AUC = 0.547** (reported transparently): the binary "edge target moved at
  adj_p<0.1" fraction (~1 %) is deflated by strict 10k-gene multiple-testing plus the low-importance
  edge tail, which is exactly why the thresholded top-K metric — not the all-edge fraction — is the
  headline. The `validated` boolean in the edge list marks the strict adj_p<0.1 hits and is a
  conservative floor, not the primary evidence.

## 7. In-silico perturbation (proxy)
**CellOracle is NOT available in this environment** (import fails: missing igraph/louvain/gimmemotifs,
and celloracle 0.20.0 pins numpy==1.26.4 / pandas<=1.5.3 / scanpy / matplotlib<3.7, all incompatible
with the arboreto/GRNBoost2 stack; force-installing them would break the running GRN). As explicitly
permitted, we use **signed linear network propagation as the in-silico perturbation proxy**: a +1
(agonism) shift on each costim node propagated through the signed GRN kernel to the program gene sets;
the KO proxy is the negation. This is a linear-propagation stand-in for CellOracle's simulate_shift,
**not** the trained CellOracle simulation — stated plainly here and in the outputs.

## 8. QSP alignment
`cd4_grn_qsp_drive.csv` carries `receptor, effector_z, crs_z, supp_z, help_z` as **empirical-null
z** where the null is the same-propagation drive of 600 random non-costim regulators, so
**0 = backbone-typical** on the QSP's empirical-null scale (matching `qsp_input_matrix_Xv2.csv`
consumed by the frozen QSP model, artifact 81b3006b). Higher crs_z/supp_z = more liability; higher
help_z = more help preserved. Verified: `arm_from_row()` in the frozen QSP model loads all 12 rows.

## 9. Three-axis nomination (CD8 efficacy × CD4 toxicity)
**Critical attribution.** This CD4 CRISPRi GRN resolves the *toxicity* axes (suppression, CRS) and a
CD4-intrinsic effector-*program* drive. The CD4-intrinsic effector drive is **not** the efficacy axis
of the nomination and must not be read as CD8 killing. The **efficacy axis is the CD8 effector benefit**
carried in the prior lane's `qsp_input_matrix_Xv2.csv` (`effector_z`) — a gain-of-function,
IFN-γ-sorted CD8 CRISPRa/CRISPRi readout (the project attributes this anchor to the Schmidt et al.
primary-human-T-cell screen; used here as the prior lane produced it, not re-derived). It supplies the
agonism-direction efficacy evidence a CD4 knockdown screen structurally cannot. The nomination pairs:
- **Efficacy:** CD8 `effector_z` from the prior lane's Schmidt-anchored matrix.
- **Toxicity:** this GRN's CD4 `supp_z` + `crs_z`.
- **Modifier:** this GRN's CD4 `help_z` (help preserved is favorable).

Because the two sources are on different z-scales (Schmidt effector_z spans up to +12; the GRN z is
clipped to ±3.5), we combine them by **within-panel rank-percentile** per axis and score the window as
the balanced mean of (CD8-efficacy↑, CD4-suppression↓, CD4-CRS↓) with a small help bonus. The rank
form is deliberately **robust to CD28's effector outlier** — a naive z-sum lets CD28's +12.11 dominate
and masks its CD4 suppression liability.

**Result (worked example — a Treg-aware arm beats the pan-costim arm on window).**
Two axes must be kept distinct: the balanced **window-rank score**, and the raw **CD4 toxicity** (supp+CRS).
- **CD4 toxicity, ascending among the CD8-effector hits:** TNFRSF8 (supp −1.65, CRS −2.24; sum −3.90)
  < CD27 (sum −0.11) < CD28 (sum +0.85) < TNFRSF9 (+1.26) < CD40 < TNFRSF4 < CD2. Both **TNFRSF8 (CD30)**
  and **CD27** are real CD8-effector hits (Schmidt +3.2 / +4.3) with **lower CD4 toxicity than CD28** —
  the Treg-aware profile: agonism lands on CD8 effectors while sparing the CD4/Treg/CRS axis.
- **Window-rank score:** TNFRSF8 (0.746) > CD28 (0.724) > CD27 (0.686). **Only TNFRSF8 outranks CD28
  on the balanced window.** CD27 is Treg-aware on the toxicity axis but does *not* outrank CD28 on the
  window score, because CD28's extreme CD8 benefit (+12.1) still lifts its balanced rank above CD27's
  (Schmidt +4.3). TNFRSF8 is therefore the single receptor that beats the pan-costim arm on window;
  CD27 is the secondary Treg-aware candidate (favorable toxicity, strong CD8 benefit, window ≈ CD28).
- **CD28** is the **pan-costim reference**: maximal CD8 benefit (+12.1) with the project's canonical
  liability profile — high CD4 **suppression** drive (supp_z +2.36; its CRS drive is low, −1.51),
  the Treg-expanding TGN1412 axis. (Its supp_z is high but not the panel maximum — TNFRSF9 and
  TNFRSF4 score +2.64; CD28 is called out because it is the clinically archetypal pan-costim arm,
  not because it is the single most suppressive node.)
- **CD226** is effector-*program*-high and CD4-toxicity-low in CD4, **but is not a CD8-effector hit**
  in Schmidt (+0.64) — the clearest illustration of why the CD4-intrinsic effector drive cannot stand
  in for the CD8 efficacy axis.
Output: `cd4grn_schmidt_integrated_nomination.csv`, `cd4_grn_nomination.png`.

## 10. Outputs
- `cd4_grn_edges.parquet` — genome-scale GRN: source_TF, target_gene, importance, pearson_r, sign, signed_importance.
- `cd4_costim_program_drive.csv` — per-receptor effector/CRS/SUPP/HELP drive + edge-corroboration metrics.
- `cd4_grn_qsp_drive.csv` — QSP-aligned empirical-null z (backbone=0).
- `cd4grn_schmidt_integrated_nomination.csv` — CD8 efficacy (Schmidt) × CD4 toxicity (this GRN), window rank.
- `cd4_grn_figure.png` — drive heatmap + costim-neighborhood network view.
- `cd4_grn_nomination.png` — efficacy-vs-toxicity nomination map + window-rank ranking.
- `cd4_grn_edges_TFtarget.csv` / `costim_edges.csv` — GRNBoost2-native `TF,target,importance,validated` edge lists.

# A26 — Tier-2 costim signaling network from CD4 Perturb-seq (STRETCH LANE)

**Status: CONVERGED (reduced scope, honestly bounded).** A fitted causal network was obtained
on the costim-receptor subnetwork; full dynamical GRN simulation (CellOracle) is scoped but not
run — see "What remains". This note states exactly what tier was achieved and why.

---

## TL;DR — tier achieved (precise)

Three linked objects, in decreasing strength of causal claim:

1. **Fitted causal network — receptor subnetwork (Modular Response Analysis, Kholodenko 2002 PNAS).**
   Solved the 11-node costim-receptor connectivity matrix `r` by deconvolving indirect effects from
   the ground-truth CRISPRi perturbation-response matrix. This IS a fitted causal-inference method,
   not a correlation graph. **Scope boundary: it is a TRANSCRIPTIONAL cross-regulation network among
   the receptor genes** (readout = receptor mRNA), NOT a signaling-flux network. Most receptors
   signal post-transcriptionally, so this layer is genuinely sparse — that sparsity is a result,
   not a failure.

2. **Directed perturbation→program map (causal-by-construction).** Each receptor KD → response of
   the CRS / SUPP / HELP program genes. Causal because the perturbed node is known ground truth;
   this is a direct readout, not a fit.

3. **Linear-regime agonism & co-agonism simulation.** Agonism = sign-flipped (−1×) KD response,
   genome-wide-null-normalized to z. Co-agonism of a pair = MRA linear superposition of the two
   single-agonism vectors. These are LINEAR-REGIME PREDICTIONS, **not** a dynamical simulation
   (no iterative signal propagation, no saturation).
   Ranking pairs by (CRS<0, SUPP<0, HELP>0): 2 of 8 qualify — 4-1BB+CD27 (best-balanced) and
   4-1BB+CD28 (higher HELP). See the co-agonism headline under Key results — the earlier "only pair"
   framing was corrected against the CSV (4 pairs lower both liabilities; 2 also preserve HELP).

**What this is NOT:** it is not a genome-scale dynamical GRN with iterative signal propagation
(the CellOracle deliverable). It is not a signaling-flux model (caveat #2 below). It does not
measure CD8 killing (that is the Schmidt CD8 axis, a separate lane).

---

## Why MRA on the perturbations, not CellOracle/SCENIC on co-expression

The task offered pip-installing pyscenic/celloracle into a forked env. I did **not** install them,
by design: SCENIC and CellOracle infer edges from **co-expression + TF motifs** and use
perturbations only as downstream validation. In *this* dataset the perturbations are the crown
jewel — 11,526 ground-truth CRISPRi knockdowns. MRA fits the network **directly from those
perturbations**, which is the more principled use of a genome-scale Perturb-seq and avoids
discarding its defining feature. MRA is also exactly right-sized: 11 clean single-node
perturbations → 11-node network is the canonical MRA setup (one perturbation per module).

Nothing was installed. Only pre-existing sc-analysis packages (scanpy/anndata/numpy/networkx 3.6.1)
were used. The pinned path-venv was not modified.

---

## Method

- **Substrate:** hero DE matrix `derived_DE/GWCD4i.DE_stats.h5ad` (33,983 regulator×timepoint rows ×
  10,282 readout genes), `log_fc` layer, opened backed='r'. All 12 costim receptors are present as
  BOTH perturbed regulators AND measured readouts; CRS(3/3)/SUPP(15/15)/HELP(5/5) all present.
- **On-target QC gate:** 11/12 receptor KDs are strongly on-target (effect −7 to −50, significant).
  **CD40 fails QC** (on-target −0.5 to −2.2, never significant, ~0 DE genes) → dropped as an MRA
  *perturbation* (kept as a downstream node only). Consistent with CD40 being APC-side, not
  cis-costim in CD4 T cells (gene-set metadata flag).
- **MRA:** per condition, R_p[i,j] = log_fc of node i under KD of node j; per-row linear solve for
  local connectivity r (r_ii=−1), which is invariant to perturbation strength (so unequal KD
  magnitudes do not bias it). Conditioning was excellent: max condition number 4.1–10.9 across all
  rows/conditions (well-conditioned <100) → stable deconvolution.
- **Programs:** CRS=[TNF,IL2,IFNG]; SUPP=15-gene Treg/IL-10; HELP=5-gene Tfh. Agonism z =
  −1×(receptor module-mean log_fc − genome-wide-null mean)/null-sd, per condition.

## Key results (re-derived from primary data)

- **Receptor transcriptional network is sparse & ICOS-hubbed.** Significant directed
  receptor→receptor edges: Rest 13 / Stim8hr 26 / Stim48hr 16 (of 132 possible). ICOS is the
  dominant source hub (5 out-edges @48hr); strongest MRA-direct edges are ICOS⊣OX40 (r=−0.47) and
  DNAM1→OX40 (r=+0.47). Most receptors act transcriptionally independently → post-transcriptional
  signaling dominates (caveat #2).
- **Agonism→program drive (48hr, z):** CD30 (−2.5) and ICOS (−2.8) most strongly *suppress* the
  SUPP/Treg program; GITR (+2.0) and OX40 (+1.5) *feed* it. CD2 (+1.9) and CD28 (+1.5) top CRS.
  4-1BB is CRS-negative (−1.6) — consistent with the prior A17 low-liability picture.
- **Co-agonism headline:** of 8 candidate pairs, 4 lower BOTH liability axes (CRS<0 & SUPP<0):
  4-1BB+CD27 (−1.38,−1.38), 4-1BB+DR3 (−2.37,−0.60), 4-1BB+CD28 (−0.14,−0.75), CD27+DR3 (−0.51,−0.38).
  Adding the benefit criterion (HELP>0) narrows this to **2 qualifying pairs: 4-1BB+CD27 (CRS −1.4,
  SUPP −1.4, HELP +1.1) and 4-1BB+CD28 (CRS −0.1, SUPP −0.8, HELP +1.7)** (the two DR3 pairs are
  dropped — DR3 co-agonism drives HELP negative). Among the two qualifiers, **4-1BB+CD27 gives the
  strongest and most balanced dual-liability reduction** (combined CRS+SUPP = −2.77 vs −0.89 for
  4-1BB+CD28), while 4-1BB+CD28 trades a larger HELP boost for weaker liability reduction. CD27+CD2
  maximizes HELP (+3.4) but crosses into CRS liability (+2.1). The network thus endorses 4-1BB+CD27
  as the best-balanced low-liability costim pair and 4-1BB+CD28 as the higher-HELP alternative —
  a ranked prediction only a network (not per-receptor scores) can make.

## Standing caveats (all apply)

1. **LOF→GOF sign-flip:** a KD that lowers a program is a *prior* that agonism raises it, not proof.
2. **Transcript ≠ signaling flux:** this is a TRANSCRIPTIONAL surrogate network. Sparse receptor
   cross-talk reflects the transcriptional layer only.
3. **RTCC/engager context:** costim receptors are activation-induced; Stim48hr is the relevant panel.
4. **CD4-side only:** effector benefit is CD8-side (Schmidt) — a complete drive vector integrates both.
5. **Linear regime:** co-agonism = additive superposition; real receptor synergy/saturation is not modeled.
6. **MRA assumes quasi-steady-state** and one perturbation per node; 48hr snapshot approximates this.

## What remains — exact next steps for the FULL dynamical GRN (CellOracle)

To upgrade from linear-regime prediction to iterative dynamical simulation (signal propagation to
convergence), the substrate is the QSP-lane steer's pseudobulk (NOT per-cell):
`data/hero_cd4_perturbseq_zhu2025/pseudobulk/GWCD4i.pseudobulk_merged.h5ad` (~44GB, tractable).

    # FORKED env — do NOT touch sc-analysis (pinned path-venv):
    conda create -n grn-tier2 python=3.10 -y
    conda activate grn-tier2
    pip install celloracle scanpy  # celloracle pulls velocyto/pyarrow/gimmemotifs
    # base GRN: CellOracle ships a human promoter-TF base (hg38, 10x scATAC or built-in).
    # 1. Load pseudobulk AnnData, standard scanpy preprocess (HVG, PCA, neighbors, Leiden).
    # 2. oracle.import_TF_data(base_GRN); oracle.get_links() -> cluster-wise GRN.
    # 3. oracle.simulate_shift(perturb_condition={receptor: high}) for each costim receptor.
    # 4. read propagated program shift (CRS/SUPP/HELP) as the dynamical agonism vector.
    # VALIDATION anchor: CellOracle edges among the 11 receptors should recover the MRA-direct
    #   edges here (ICOS-OX40, DNAM1-OX40); the ground-truth KD map is the held-out test set.

Est. effort: day-scale (base-GRN setup + get_links on 44GB pseudobulk is the bottleneck).

## Artifacts
- A26_grn_network_48hr.png — directed receptor perturbation-response network + MRA direct/indirect
- A26_grn_agonism_simulation.png — (A) receptor×program agonism heatmap, (B) co-agonism pairs
- A26_grn_receptor_receptor_edges.parquet — 396 directed edges (3 conds × 132), z + MRA-direct r
- A26_grn_receptor_program_edges.parquet — 828 receptor→program edges, KD + agonism-signed
- A26_grn_agonism_simulation.csv — per-receptor CRS/SUPP/HELP agonism z, 3 conditions
- A26_grn_coagonism_prediction.csv — 8 pair predictions (MRA superposition, 48hr)
- A26_grn_qsp_drive_tier2.csv — QSP from_row-compatible Tier-2 network drive vector
- A26_grn_mra_connectivity_{Rest,Stim8hr,Stim48hr}.csv — fitted MRA connectivity matrices

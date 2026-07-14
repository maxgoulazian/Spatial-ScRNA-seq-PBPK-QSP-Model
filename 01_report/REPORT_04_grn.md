# 4. Gene-regulatory-network mechanism & QSP drive operators

The three-axis nomination (Section 3) scores each costimulatory arm as a set of scalar
liability and benefit values. This section derives the mechanism *underneath* those scalars —
a genome-scale gene-regulatory network (GRN) inferred from the hero CD4⁺ CRISPRi Perturb-seq —
and converts it into the **drive operator** the QSP lane consumes: a signed, receptor-anchored
propagation kernel that maps each candidate agonism into a magnitude, a program shape, and a
timescale for the downstream dynamical model. The section closes on the single mechanistic
result that carries the nomination story — the contrast between the CD8-effector and CD4
regulatory rings.

The QSP model outputs themselves (therapeutic-window shifts, dose rationale) are owned by a
separate lane and are **not** reported here; this section defines the operator and hands it off.

---

## 4.1 The genome-scale CD4 regulatory backbone

The substrate is the Marson/Pritchard genome-scale CD4⁺ T-cell CRISPRi Perturb-seq
(Zhu et al. 2025, bioRxiv; dataset used as provided — the underlying accession is a dataset
tag, not independently re-verified here). From 278,684 pseudobulk guide×donor×condition
profiles we restricted to the **Stim48hr** sustained-activation state (the state in which
costimulatory wiring and the CRS/suppressive programs are engaged) at `keep_for_DE=True` QC-pass,
and collapsed to one expression profile per CRISPRi target — **11,332 per-perturbation profiles
× 18,127 genes** (CP10K + log1p, Ensembl-ID space). The perturbation-driven variance across
knockdowns is precisely the signal GRN inference exploits.

Regulators were the human TF catalogue (Lambert et al. 2018, filtered to `is_TF=Yes`;
1,569 TFs) unioned with the 12 costimulatory receptors and the program genes, for
**1,598 regulators**. Edges were inferred with **GRNBoost2** (`arboreto`; per-target
stochastic gradient-boosting regression, feature importances = directed regulator→target
weights), the gradient-boosting core of the SCENIC pipeline (Moerman et al. 2019,
doi:10.1093/bioinformatics/bty916). GRNBoost2 importances are unsigned; each edge was signed by
the Pearson correlation of regulator and target across the 11,332 profiles — the CellOracle
convention, `signed_importance = sign(r) × importance` (Kamimoto et al. 2023,
doi:10.1038/s41586-022-05688-9).

**The backbone is a directed, signed, genome-scale network of 1,653,594 edges over
1,598 regulators × 18,127 target genes** (`grn_operator_shared_backbone.parquet`). Edge
correlations span r ∈ [−0.781, +0.944]; unsigned importances reach 44.67; **59.2 % of edges
are activating and 40.8 % repressing.** This shared backbone is the *single* regulatory
substrate over which every per-arm drive is computed — it is, in the QSP framing, the
"CD3-only" wiring of the resting-to-activated CD4 compartment before any costimulatory arm is
applied.

### Perturb-seq validation of the inferred edges

Because GRNBoost2 edges are co-expression-derived, we validated them against the independent
CRISPRi differential-expression matrix (each receptor's *own* knockdown signature, Stim48hr).
The network passes on all three complementary criteria reported in `cd4_grn_methods.md`:

- **Top-25 edge enrichment** — mean AUC **0.583** across 420 CRISPRi-targeted regulators,
  **69.8 %** with AUC > 0.5 (p = 1.7×10⁻¹⁶): a regulator's strongest GRN targets are
  systematically more perturbed by its own knockdown than random genes.
- **Importance → knockdown magnitude** — mean within-TF Pearson r = **+0.097**, positive for
  **78.7 %** of broad-knockdown regulators (p = 7.6×10⁻¹⁰).
- **Edge sign → knockdown direction** — **61.7 %** concordance among moved edges
  (p = 3.3×10⁻¹⁹⁸).

The strict per-edge `validated` boolean (adj_p < 0.1 in the DE matrix) marks **8,457 edges
(0.51 % of the backbone)** and is retained as a conservative floor; the thresholded top-K
enrichment — not the all-edge fraction (AUC 0.547, deflated by 18k-gene multiple testing) — is
the headline evidence that the inferred edges reflect real causal knockdown effects.

---

## 4.2 Receptor → program drive: the propagation operator

The operator converts each receptor into a signed drive onto four gene programs. Program sets
(all present in the backbone) are the effector/cytotoxic set (GZMB, PRF1, GZMK, NKG7, IFNG,
GZMA, GNLY, FASLG), the CRS cytokines (TNF, IL2, IFNG), the suppressive/Treg program (FOXP3,
IL10, CTLA4, IKZF2, ENTPD1, IL2RA, TGFB1, LRRC32) and the help program (BCL6, MAF, TOX2, CXCR5,
TCF7).

For each receptor the drive is computed by signed network propagation from the receptor node
through a bounded multi-hop kernel **K = I + α·T + α²·T²** (α = 0.5) over the regulator subspace,
followed by a row-normalized regulator→gene hop. Directionality matters: the propagation
represents a **+Δ agonism shift** on the receptor node — the gain-of-function direction an
engager arm acts on. Because the CRISPRi screen is loss-of-function, we apply the LOF→GOF sign
flip explicitly (**agonism = −1 × knockdown-direction**); all reported drive values are on the
agonism axis, which is the axis the QSP model consumes. This linear signed-propagation scheme is
used as the in-silico-perturbation proxy in place of a trained CellOracle simulation (CellOracle's
compiled dependency stack is incompatible with the GRNBoost2 environment here); it is a
propagation stand-in for `simulate_shift`, stated as such, not the trained simulator.

Drive is reported as **empirical-null z**, where the null is the same-propagation drive of 600
random non-costimulatory regulators, so **0 = backbone-typical** (`cd4_grn_qsp_drive_full.csv`;
this is the matrix visualized in the heatmap of Figure 4.1a and consumed by the QSP lane). On
this scale, higher CRS-z / suppression-z is more liability and higher help-z is more help
preserved. The co-lead values are:

| Arm | effector-drive z | CRS-drive z | suppression-drive z | help-drive z |
|---|---|---|---|---|
| 4-1BB (TNFRSF9) | +1.33 | **−1.38** | +2.64 | +1.27 |
| CD27 | +1.98 | **−1.38** | +1.27 | +1.27 |
| CD28 (pan-costim ref.) | −2.43 | −1.51 | +2.36 | +2.52 |
| CD30 (TNFRSF8) | −1.82 | −2.24 | −1.65 | −1.69 |

Both co-leads carry a **favorable (negative) CRS-drive** on this propagation axis. Their
suppression-drive is positive (the receptors are wired, through their TRAF/NF-κB targets, toward
some suppression-program genes more than a random regulator) — but the network-propagation drive
is a *mechanism-magnitude for the dynamical model*, **not** the liability gate. The CLEAN/GATED
liability calls in the nomination (Section 3) are set by the direct CRISPRi DE screen in the T1
lane, on which both 4-1BB and CD27 are CLEAN; the GRN operator supplies the shape and magnitude
of the drive, not the veto.

![Figure 4.1 — GRN-derived costim→program drive (a: empirical-null z heatmap, 0 = CD3-only backbone) and the signed costim-receptor→program-gene neighborhood (b). Genome-scale GRN: 1,653,594 edges, 1,598 regulators × 18,127 targets, GRNBoost2 on CD4 CRISPRi Perturb-seq Stim48hr.]({{artifact:ff249e44-d7d8-42da-a341-deb178c96210}})

---

## 4.3 Receptor input-node map (43/43)

The propagation cannot start on the receptor transcript alone — a surface receptor acts through
its proximal signaling adapters, and those are the nodes that carry the signal into the TF
backbone. We therefore mapped each of the 12 receptors to its canonical proximal-signaling input
nodes and confirmed every node is present in the backbone (`receptor_input_node_map.csv`):
**43 of 43 input nodes are in the backbone (100 %; every arm complete).** For the TNF-receptor
superfamily co-leads the input nodes are the TRAF→NF-κB signalosome:

- **4-1BB (TNFRSF9)** → TRAF1, TRAF2, NFKB1, RELA. CD137 has no intrinsic enzymatic activity and
  relies on TRAF adaptors to build the CD137 signalosome: TRAF1, TRAF2 and TRAF3 are recruited to
  its cytoplasmic domain, and TRAF-mediated signaling activates NF-κB (Arch & Thompson 1998,
  Mol. Cell Biol. 18:558–565, doi:10.1128/MCB.18.1.558; Zapata et al. 2018, Front. Immunol.
  9:2618, doi:10.3389/fimmu.2018.02618, PMID 30524423). The operator anchors 4-1BB on the
  NF-κB-activating TRAF1/TRAF2 pair plus the NF-κB subunits NFKB1/RELA.
- **CD27** → TRAF2, TRAF3, NFKB1, MAP3K14 (NIK). CD27 engages TRAF adapters and the
  NF-κB-inducing kinase NIK (MAP3K14) to activate NF-κB: the CD27 cytoplasmic PIQEDYR motif is
  required for interaction with TRAF2 and TRAF5, dominant-negative TRAF2 or TRAF5 blocks
  CD27-induced NF-κB activation, and NIK is the common downstream kinase (Akiba et al. 1998,
  J. Biol. Chem. 273(21):13353–13358, doi:10.1074/jbc.273.21.13353, PMID 9582383). Note the
  operator anchored CD27 on TRAF2/TRAF3 as the TRAF-adapter pair rather than the canonical
  TRAF2/TRAF5 — both are TNFRSF TRAF paralogs, but this is a node-choice detail worth flagging
  (see Unresolved).

Anchoring on the validated adapter nodes rather than the receptor transcript is what makes the
per-arm drive a *signaling* drive rather than a co-expression artifact, and it is why the operator
is receptor-specific even where two receptors share downstream TFs.

---

## 4.4 Per-arm drive magnitude, uncertainty, and timescale (the QSP handoff)

The QSP model needs each arm's drive as three separable quantities — **magnitude** (how large a
program shift), **shape** (which program), and **timescale** (how fast). The
`per_arm_drive_magnitude_uncertainty.csv` operator provides all three with bootstrap uncertainty
(`per_arm_prolif_exh_drive.csv` adds the proliferation and exhaustion programs). Two features are
load-bearing for the handoff:

**Bootstrap uncertainty on the effector benefit is not uniform across arms.** On the Schmidt
CD8-effector axis (the efficacy anchor; Section 3), the two co-leads have nearly identical
point estimates — **CD27 effector_z = 4.05 (SD 0.74; tier `T1_favorable_robust`)** and
**4-1BB effector_z = 4.08 (SD 2.38; tier `T3_neutral_or_insufficient`)** — but very different
robustness. CD27's effector benefit is tight and robust; **4-1BB's propagation-drive effector
estimate is uncertain** (wide bootstrap). This is exactly the "complementary evidence" structure:
4-1BB's CD8-selectivity is not carried by its propagation-drive magnitude (which is uncertain)
but by its **wiring topology** (Section 4.5), while CD27's is carried by a robust drive. The
pan-costim reference CD28 is high-magnitude (effector_z = 11.94, SD 0.84) but `T4_liability`.

**Timescale is resolved.** The operator carries 8-hour and 48-hour drive z for the CRS and
suppression programs (the headline CRS/suppression columns equal the 48-hour values). This
8 h→48 h split is the timescale term the QSP consumes to separate fast cytokine kinetics from
slower suppressive-program consolidation — for example CD2's CRS drive falls from +4.43 (8 h) to
+1.85 (48 h), and CD30's suppression drive deepens from −1.14 (8 h) to −2.50 (48 h). The QSP lane
scales the effector-drive magnitude into the effector-activation arm through a single
handoff coupling scalar (**kE = 0.11**, the GRN→QSP effector gain specified for this integration;
see §4.6 and Unresolved).

*Estimator note.* This table's CRS/suppression z use a bootstrap + timepoint estimator that
differs in scale from the empirical-null propagation z of §4.2 (`cd4_grn_qsp_drive_full.csv`); the
two are not directly comparable in absolute value (see Unresolved). Absolute per-program drive
magnitude for the QSP should be read from the operator the QSP lane is wired to; this table is
used for the **uncertainty and timescale** terms.

---

## 4.5 The two effector-lineage rings: the contrast is the story

The operator was instantiated on two matched effector-lineage subnetworks — a CD8-effector ring
and a CD4 ring — each holding the 12 costimulatory receptors as sources and the effector program
as the target set (`grn_operator_cd8_effector_lineage.parquet`, 4,881 edges;
`grn_operator_cd4_effector_lineage.parquet`, 4,934 edges; 86 sources × 98 targets each). Node
fill = each receptor's net signed drive onto the effector program; the layout is identical
between rings, so any difference is wiring, not geometry.

**In the CD8-effector ring, 4-1BB is the dominant effector hub; in the CD4 ring, that same wiring
collapses.** 4-1BB's single strongest receptor→effector edge is TNFRSF9→GZMB (coef +24.24 — the
largest of any receptor→effector edge in the CD8 ring), with TNFRSF9→IFNG (+8.95) and
TNFRSF9→NKG7 (+8.02) behind it; in the CD4 ring its top effector edge is only TNFRSF9→FASLG
(+1.42). Aggregated to net signed effector drive, the network-selectivity operator
(`A29_network_effector_selectivity.csv`) scores **4-1BB at net +54.36 — CD8 drive 57.4 vs CD4
drive 3.04, a ~19× CD8-vs-CD4 bias, rank #1 of 12** and the single dominant CD8-selective
effector hub. An independent recomputation of the receptor→effector drive directly from the two
ring operators reproduces the direction and magnitude ratio (4-1BB: CD8 ≈ 43.7 vs CD4 ≈ 2.1).
CD27, by contrast, is lineage-balanced by wiring (CD8 3.97 vs CD4 3.48, net +0.49); its
CD8-bias is carried by *where it is expressed* rather than *how it is wired* (Section 3), which
is the complementary half of the co-lead pair. Note also that the two rings are drawn on
different drive scales — the CD8 ring spans ≈ ±9 effector-drive units while the CD4 ring spans
only ≈ ±1.2 — so the effector amplitude itself collapses roughly an order of magnitude in the
CD4 lineage even before per-receptor selectivity is considered.

This is the mechanistic statement behind "the boost lands on the killers": agonizing 4-1BB drives
the cytotoxic-effector program through a network that is wired for it in CD8 cells and largely
inert in CD4 cells, so a 4-1BB costimulatory arm amplifies effector output where killing happens
while contributing little drive to the CD4 compartment where the suppressive and cytokine-release
liabilities live.

![Figure 4.2 — Costim-receptor GRN, CD8-effector ring. Directed edges = GRNBoost2 regulatory links among costim arms; node fill = effector drive score (scale ≈ ±9). 4-1BB is the darkest-red node (maximal effector drive).]({{artifact:b684c6ce-c5bb-40a2-a8bc-4fda4ce4e798}})

![Figure 4.3 — Costim-receptor GRN, CD4 ring, identical layout (node fill scale ≈ ±1.2). The 4-1BB effector-hub wiring seen in the CD8 ring is absent; effector amplitude collapses ~7–8× across the lineage.]({{artifact:22f9edd2-f72f-4317-9402-f2e77fe2dc85}})

---

## 4.6 Integration handoff to the QSP lane

The GRN delivers to the QSP a **drive operator**, not a set of endpoints. Concretely, the handoff
is the receptor-anchored propagation operator over the 1,653,594-edge shared backbone (§4.1–4.2),
instantiated per arm through the 43/43 validated input-node map (§4.3), and decomposed into the
three terms the dynamical model consumes:

1. **Magnitude** — per-arm effector/CRS/suppression/help drive, coupled into the effector-activation
   arm by the handoff scalar kE = 0.11;
2. **Shape** — the four-program signed vector (which program each arm feeds), the axis that
   distinguishes a Treg-aware arm from a pan-costim arm;
3. **Timescale** — the 8 h vs 48 h drive split (§4.4), separating fast cytokine kinetics from
   slower suppressive-program consolidation;

each carried with its bootstrap uncertainty so the QSP can propagate confidence. The QSP lane
parameterizes the T-cell activation and PD arms from these terms and reports the therapeutic-window
consequences; those results are out of scope for this section by design.

---

### Methods artifacts (this section)

- `grn_operator_shared_backbone.parquet` — signed genome-scale backbone (1,653,594 edges).
- `cd4_grn_qsp_drive_full.csv` — per-arm empirical-null drive z (heatmap / QSP source).
- `per_arm_drive_magnitude_uncertainty.csv`, `per_arm_prolif_exh_drive.csv` — magnitude + bootstrap
  SD + 8 h/48 h timescale; proliferation/exhaustion drive.
- `receptor_input_node_map.csv` — 43/43 proximal-signaling input nodes.
- `grn_operator_cd8_effector_lineage.parquet`, `grn_operator_cd4_effector_lineage.parquet` — the two
  effector-ring operators.
- `cd4_grn_methods.md` — full inference, signing, propagation, and CRISPRi-validation methods.
- Figures: `cd4_grn_figure.png`, `grn_ring_cd8_effector.png`, `grn_ring_cd4_effector.png`.

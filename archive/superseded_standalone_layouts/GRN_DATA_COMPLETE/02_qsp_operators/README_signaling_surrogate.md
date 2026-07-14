# COSTIM GRN SIGNALING SURROGATE — propagation operator for the QSP
_generated 2026-07-12T04:03:25Z · consumer: per-cell costim signaling integration in whole-body QSP_

CellOracle-style: inject occupancy-scaled perturbation at each receptor's input node,
propagate the signed shift through the GRN coefficient matrix, score per-cell program
shifts (prolif/exh/cytokine/supp/effector) as continuous rate modifiers. Veto gate
(PROLIF/EXH hard cut) stays intact — these shifts complement it, don't override it.

## ARCHITECTURE (honest — read before wiring)
The hero is a CD4 Perturb-seq, so there is NOT a separately-fit genome-scale GRN per
cell type. Cell-type identity enters two documented ways:
  - SHARED genome-scale signed backbone (hero CD4) = the propagation operator. It is the
    ONLY network reaching all 5 program modules AND all 12 receptor input nodes.
  - LINEAGE effector-subnetwork refinements (CD8, CD4-GOF from Schmidt CRISPRa lineage
    split) = signed, but effector-restricted (98 nodes: effector+cytokine ONLY; ZERO
    prolif/exh/supp nodes). Use to refine the EFFECTOR arm of propagation per lineage;
    do NOT expect toxicity-axis coverage from them.
  - Cell-type SPECIALIZATION = baseline state (deliv 5, per cell type) + per-arm drive
    amplitude (deliv 4). Propagation computes shift-from-baseline; the baseline is what
    makes CD8 vs CD4-conv vs Treg differ under the shared topology.

## FILES
(1) grn_operator_shared_backbone.parquet  — source_TF,target_gene,coef(signed),coef_corr,importance,validated
      1,653,594 edges; coef range [-22.9,+44.7]; 41% negative. THE propagation operator.
    grn_operator_cd8_effector_lineage.parquet / grn_operator_cd4_effector_lineage.parquet
      — signed effector-subnetwork refinement (~4.9k edges each, 98 nodes, effector+cytokine only).
(2) receptor_input_node_map.csv  — 12 arms → immediate downstream signaling nodes; 43/43 present in backbone.
(3) program_gene_sets.json  — prolif/exh/cytokine(split IL6,IFNG,TNF,IL2)/suppression/effector; symbols match operator.
(4) per_arm_drive_magnitude_uncertainty.csv  — *_z + *_z_sd (concordance CIs) per arm/axis; kinetic 8hr/48hr cols. kE=0.11 LOCKED.
(5) baseline_CD4conv.parquet (hero baseMean, GRN-consistent) / baseline_CD8.parquet / baseline_Treg.parquet
      (RTCC CITE-seq gated: CD8 13374, Treg 423, CD4conv from hero). gene,baseline_mean_expr,cell_type,source.

## CONTRACT CHECKS
  - gene symbols consistent across (1)/(3)/(5): operator symbols ⊇ all program genes (verified).
  - coefficients signed: yes (deliv 1 coef; sign from pearson_r on the fitting data).
  - one baseline file per cell type: yes (CD8/CD4conv/Treg).
  - veto gate untouched: propagated shifts are continuous rate modifiers only.


## ADDENDUM 2026-07-12T04:52:53Z — answers to QSP lane's 3 questions
(1) PROPAGATION READOUT — CONFIRMED, your split is exactly the intent.
    magnitude = deliverable-4 z-scores (drive table); per-cell shape = normalized
    GRN-propagated shift × baseline; occupancy scales over PK time. Do NOT read absolute
    program magnitude from genome-scale propagation (it dilutes to ~0, as you found — the
    backbone is a shape/heterogeneity operator, not a magnitude operator). Drive-table-as-
    magnitude is correct and is what the calibration (Schmidt kE=0.11 locked) is anchored to.

(2) CD4conv BASELINE SCALE — FIXED, redelivered on CITE-seq scale.
    baseline_CD4conv.parquet is now RTCC-gated (2,057 CD4conv cells, same log-CP10K
    normalization as CD8/Treg). Σbaseline: CD8 2851 / Treg 2741 / CD4conv 2635 (comparable).
    Cross-type suppression baseline now biologically correct: Treg 0.489 >> CD8 0.090 >
    CD4conv 0.054. Use these three as the source of truth — no CP10K re-fix needed on your side.

(3) PROLIF/EXH MAGNITUDE — NEW deliverable so you don't have to infer it from flat kinetics.
    per_arm_prolif_exh_drive.csv: per-arm prolif_agon/exh_agon (agonism effect size) + sd
    (from BH 95% CIs) + q + call. This is ARM-SPECIFIC continuous magnitude the flat kinetics
    miss. NO double-count: the binary PROLIF/EXH veto gate still owns the hard cut; this table
    is the continuous per-arm RATE modifier among arms (use it to scale, not to gate).
    Key values: CD28 prolif +0.078 (q=0.009, up*conc liability); 4-1BB prolif -0.040/exh -0.041
    (both ns → ~0 modifier); CD27 prolif +0.004 ns, exh -0.066 (q=0.001, exhaustion-NEGATING
    = FAVORABLE, CD27 protects against exhaustion). Adopt this as authoritative prolif/exh
    magnitude; keep your hero k_on/k_off for the TIMESCALE. That removes the double-count risk.

# T2b — Per-cell receptor layer: panel-native extraction, scVI overlay imputation, RNA→copies conversion

**Status: LOAD-BEARING METHOD. This is the step that makes the whole-body per-cell model runnable.**
Every organ/tumor agent (a real Xenium cell) must carry a per-cell copies/cell value for *every* antibody
target arm (TAA, CD3, costim) — otherwise the per-cell Rhoden binder (T3) and synapse kill (T5) have no
receptor density to solve against. Xenium panels are small (~250–541 genes) and rarely contain all targets,
so the missing targets are **imputed onto every cell** from a matched single-cell reference. This document is
the definitive record of how that is done, because the organ-overlay *builder script* was not retained (see
§5 Retention gap) — this file + the surviving tumor builder `tumor_builds/build_cancer_abm.py` ARE the recipe.

---

## 1. Why this step exists (the problem it solves)

A Xenium slide measures only the genes on its probe panel. A target that is **off-panel** (e.g. EPCAM on the
CRC panel, or any costim receptor not designed into the probe set) has **no transcript for any cell** — the
column is a hard zero, not low expression. Running the model against that zero would silently report "no
receptor / no binding / no toxicity" for that arm in that tissue, which is a data gap masquerading as biology.

The overlay fixes this: it borrows the off-panel gene's expression from a **real reference cell that is
transcriptionally most similar** to each query cell, so every cell ends up with a biologically plausible value
for every target. This is what the user means by *"the tabula sapiens overlay makes it so all receptors are
present in every cell — that's the whole point of training."*

---

## 2. The three-tier receptor pipeline (per cell, per gene)

```
                 is the gene on the Xenium panel for this tissue?
                        │                        │
                       YES                       NO
                        │                        │
          ┌─────────────┘                        └──────────────┐
   TIER 1 native extraction                 TIER 2 scVI overlay imputation
   (CP10K-log1p)                            (topk_spread, §3)
          │                                        │
          └──────────────┬─────────────────────────┘
                         ▼
             per-cell EXPRESSION value (log-normalized)
                         │
                 TIER 3 RNA→copies conversion  (rna_to_receptor / convert_copies_ALL, §4)
                         │
                         ▼
              per-cell RECEPTOR COPIES/CELL  ('<GENE>_copies' column)
```

- **Tier 1 — native panel extraction (CP10K-log1p).** For a gene ON the panel: normalize each cell to
  1e4 total counts, `log1p`. `Mc = Mc·(1e4/lib); Mc.data = log1p(Mc.data)` (`build_cancer_abm.py:50`).
  This is exact and deterministic — verified byte-reproducible (CEACAM5/ERBB2 in the ovary build reproduce
  at corr = 1.0000).
- **Tier 2 — scVI overlay (§3).** For a gene OFF the panel: impute from the matched reference. This is the
  "training" step.
- **Tier 3 — copies conversion (§4).** Both tiers yield an *expression* value; the binder needs *copies/cell*.
  Path A (HPA-IHC-anchored genes) and Path B (HPA single-cell nCPM→Glassman ladder) convert
  expression → copies. Native "ORIG7" genes (MS4A1, CD19, TNFRSF17, GPRC5D, EPCAM, DLL3, CD3E) are stored
  directly as copies by an earlier native pipeline.

---

## 3. The overlay method — `topk_spread` (scVI co-embed + softmax-weighted reference transfer)

Canonical source: `tumor_builds/build_cancer_abm.py` lines 60–92 (tumor; method tag `<cancer>ref_topk_spread_K15`).
Organs used the identical algorithm (method tag `topk_spread_K15_softmax`).

**Inputs.** Query = the Xenium ABM cells for the tissue. Reference = a matched single-cell dataset pulled from
the CZ CELLxGENE Census (`census_version = "2025-11-08"`, `dataset_id = cfg["ref_did"]`) — Tabula Sapiens for
normal organs, a cancer-matched atlas for each tumor.

**Step 1 — scVI co-embedding.** Concatenate reference + query on shared genes; select 2000 HVGs
(`seurat_v3`, batched by source); train an scVI model (`n_latent=20, n_layers=2`) on raw counts for
`max_epochs=60` on CPU with `early_stopping`, `seed=0`. Extract the latent representation Z; split into
`Zref` (reference cells) and `Zq` (query cells). scVI's batch correction is what lets a query cell and a
reference cell be compared in a shared, technology-corrected latent space.

**Step 2 — K-nearest-neighbor assignment (K=15) with softmax weighting.**
```
nn   = NearestNeighbors(n_neighbors=15).fit(Zref)
dist,ind = nn.kneighbors(Zq)                 # 15 nearest reference cells per query cell
tau  = median(dist[:,1:]) + 1e-9             # adaptive bandwidth
w    = exp(-dist/tau); w /= w.sum(axis=1)    # softmax over the 15 neighbors
pick = (rand < cumsum(w)).argmax(axis=1)     # ONE neighbor drawn ∝ w (stochastic, seed=0)
assign = ind[arange(n), pick]                # the chosen reference cell index, PER query cell
```
The stochastic draw (rather than deterministic argmax) *spreads* assignments across many reference cells,
preserving expression heterogeneity instead of collapsing every query cell onto a single nearest neighbor.
`ref_cell_assign = assign` is stored per query cell — this is the retained provenance of the overlay.

**Step 3 — transfer off-panel expression.** For each off-panel target present in the reference,
`d[target] = refX[assign]` — each query cell copies its assigned reference cell's (normalized) expression of
that gene. `antigen_overlay_added` records what was transferred; `antigen_still_missing` records targets
absent from the reference too.

**Step 4 — cell-type labels (same latent).** `ref_celltype = ref.cell_type[assign]`; the malignant flag
`is_tum` is set when the assigned reference label contains malign/tumor/carcinoma/epithel. The overlay thus
also supplies the per-cell cell-type used everywhere downstream (kill accounting, cytokine sources, figures).

**Determinism.** `numpy` RNG seeded 0 and `scvi.settings.seed=0`, so a re-run reproduces the same assignment
*given the same reference cell ordering* — which is exactly the caveat in §5.

---

## 4. RNA→copies conversion (Tier 3, `rna_to_receptor.py` + `convert_copies_ALL.py`)

- **Path A — HPA-IHC-anchored (24 genes):** `copies = HPA_IHC_anchor(gene, organ, lineage) × per-cell
  within-lineage rank modifier (0.1–3×)`. Anchor map: `model/params/hpa_ihc_anchor_map.json`.
- **Path B — no-IHC genes (HPA single-cell nCPM → Glassman ladder):**
  `copies = anchor_nCPM(gene, lineage) × same rank modifier`, where the per-lineage baseline comes from HPA
  single-cell nCPM converted via `log10(copies)=0.7768·IHC+3.9723`. Anchor: `noihc_ncpm_anchor.json`.
  DLL3, PMEL, EPCAM, TNFRSF9, FCRL5, FLT3, … live here (DLL3/PMEL/EPCAM added 2026-07-13 from HPA release 25.1;
  see `dll3/pmel/epcam_ncpm_anchor_addition.json`, each with a `sources` block of verbatim HPA nCPM values).
- **ORIG7 native (MS4A1, CD19, TNFRSF17, GPRC5D, EPCAM, DLL3, CD3E):** already copies/cell in the organ ABMs
  from an earlier native pipeline; not re-converted.
- **Organ pool (nM) consumed by the PBPK sink:**
  `pool_nM = (Σ_cells copies · scale) / (Vis_organ · N_Avogadro) · 1e9`, where `scale = organ_true_cells /
  section_cells`. Verified: kidney EPCAM reproduces 247.637 nM exactly through this formula. **Basis-check reconciliation:** the raw-formula path is exercised by the kidney EPCAM check (247.637 = 247.637, exact). The spleen CD20 check returns **196.6 nM**, which is the *committed corrected* pool (grid-derived, decision `spleen_CD20_corrected`) — it takes the `corrected_pools` override branch (returns the stored value directly, NOT the raw density formula), so it is a pass against the corrected number, not an 18% miss against the stale ~239-241 builder-docstring figure.

**[PROVENANCE FLAG]** These are scRNA-seq transcripts converted to protein copies — *derived*, not measured
surface densities. The conversion is HPA/Glassman-anchored but not independently QIFIKIT-validated per gene
(see T7 §3.3). Do not label figure legends "measured receptor copies."

---

## 5. Retention gap — READ THIS (the "retain everything" lesson)

**What was retained:** `ref_cell_assign` (per-cell reference index) and `overlay_method` in every organ/tumor
ABM npz; the tumor builder `build_cancer_abm.py`; all `<gene>_copies` columns that were converted.

**What was NOT retained:** (a) the **organ** overlay builder script (only the tumor one survives — same
algorithm), and (b) the **reference AnnData** that `ref_cell_assign` indexes into. Because the reference was
pulled live from the Census and not saved, the integer indices in `ref_cell_assign` are only meaningful
against a reference materialized in the *identical row order* — which a fresh Census pull does not guarantee.

**Consequence:** off-panel targets that were NEVER added to the organ overlay set (EPCAM, DLL3, PMEL) cannot be
imputed onto the organ cells by a fast index-into-cached-reference; doing it correctly requires re-running the
scVI overlay (hours, 12 organs). Their **tumor** copies ARE complete (imputed/converted from the tumor builds,
which retain their columns). The consequence is bounded to **off-tumor organ toxicity scoring for three
solid-only arms** (catumaxomab/EPCAM, tarlatamab/DLL3, tebentafusp/PMEL); the validated heme panel and the
CRC solid pair (CEACAM5, ERBB2) have full organ coverage and are unaffected.

**Fix for next time (recorded):** persist the reference AnnData (`ref.write_h5ad`) alongside every ABM, and
snapshot the organ overlay builder into the repo, so `ref_cell_assign` remains resolvable and any new target
can be overlaid by a pure lookup with no scVI re-train.

---

## 6. Where it runs / files

| Concern | File:location |
|---|---|
| Native extraction + scVI overlay + labels + BEC/LEC | `tumor_builds/build_cancer_abm.py` (steps 1–6) — canonical surviving recipe |
| Copies conversion driver | `model/engine/convert_copies_ALL.py` |
| Anchor logic (Path A/B, rank modifier) | `model/engine/rna_to_receptor.py` |
| IHC anchor map (24 genes) | `model/params/hpa_ihc_anchor_map.json` |
| no-IHC nCPM anchor (Path B) | `model/rundir/handoff/noihc_ncpm_anchor.json` (+ dll3/pmel/epcam additions) |
| Per-cell copies consumed by PK/PD | `coupled_percell_pk.py:53-59`, `coupled_percell_pd.py:26-43` |
| Overlay provenance (per ABM) | npz keys `ref_cell_assign`, `overlay_method`, `antigen_overlay_added`, `antigen_still_missing` |

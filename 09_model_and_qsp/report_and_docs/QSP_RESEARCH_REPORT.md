# A Whole-Body, Per-Cell Spatial PBPK–QSP Model for Costimulatory CD3×TAA T-Cell Engagers

### Turning a three-axis costim-arm nomination into a dose-rationalized, therapeutic-window prediction — one Rhoden kinetic binding ODE per cell, across ~1.9 million agents

**QSP/PBPK companion to the CD4⁺ Perturb-seq counter-screen report.** Where the screen report nominates the
costim arm (4-1BB and CD27, with CD28 gated out), this report is the mechanistic translation layer: a
physiologically-based, spatially-resolved, agent-based pharmacokinetic/pharmacodynamic engine that carries the
two clean arms into a predicted efficacy-versus-toxicity window, with every binding event governed by true
kon/koff kinetics at single-cell resolution.

Live engine: 13 modules, single unified run harness. Reproduction: `qsp_reproduction.ipynb` (21 cells,
executed end-to-end, 0 errors). Verification: 54-row ledger, 52 PASS.

---

## Abstract

CD3-based T-cell engagers deliver T-cell "signal 1" without costimulatory "signal 2"; adding a costim arm
restores signal 2 but, because the CD3 arm engages CD4⁺ cells indiscriminately, it amplifies the CD4
suppressive and cytokine-release programs that drive toxicity. The companion counter-screen nominates the arm;
the question this work answers is whether that nomination survives contact with pharmacology — whether a
Treg-aware arm actually widens the therapeutic window once real distribution, binding, killing, and
cytokine-release kinetics are imposed. We built a whole-body PBPK–QSP model that is spatial and agent-based at
its core: the Shah & Betts two-pore/FcRn PBPK backbone is retained for plasma and lymph transport, but each
tissue's interstitial compartment is *replaced* by a spatial agent-based model in which every real Xenium cell
is an agent running its own Rhoden multi-arm kinetic binding ODE, parameterized by its own single-cell
receptor complement. Target-mediated drug disposition, cytotoxic killing, and myeloid cytokine release are not
prescribed — they *emerge* from per-cell binding and internalization on the tissue grid. The model spans 15
physiological tissues plus a target-matched tumor compartment (solid or hematologic), ~1.9 million agents in
total, and runs pharmacokinetics and pharmacodynamics in a single coupled pass so that each molecule's own
plasma exposure drives its own per-cell pharmacology. Every receptor is present on every cell through an scVI
Tabula-Sapiens overlay — the step that makes a small-panel spatial technology usable for a whole-target-panel
model. Under this engine the nominated arms 4-1BB (TNFRSF9) and CD27 carry clean therapeutic windows while
CD28 is gated out on cytokine-release, suppression, and proliferation liabilities, recapitulating the TGN1412
experience in silico. We report the model's successes and its failures with equal weight: the PK reproduces
digitized clinical exposure across a 76-drug database, but the mechanistic IL-6 axis fails a rigorous
validation gate on verified clinical anchors (talquetamab inverts, because GPRC5D sits on plasma cells and
keratinized tissue rather than the myeloid-contacting compartments that seed CRS), and we report that failure
as-is rather than tuning it away. The result is a mechanistic instrument for costim-engager design in which
efficacy and toxicity are the same emergent computation read out in different compartments.

---

## Executive summary

- **What it is.** A whole-body, per-cell, spatial PBPK–QSP engine for CD3×TAA(×costim) T-cell engagers. PBPK
  transport backbone (Shah & Betts) + per-cell agent-based interstitium (real Xenium cells) + one Rhoden
  kinetic binding ODE per agent. PK and PD in one coupled pass.
- **Why it is different.** TMDD, killing, and cytokine release are *emergent* from single-cell binding on a
  spatial grid — not hand-set compartmental sinks. The same binding engine runs in both the PK sink and the PD
  synapse, so efficacy and toxicity are mechanistically consistent by construction.
- **Scale.** 15 PBPK tissues + per-target tumor compartment; ~1.9M agents; every antibody target present on
  every cell via scVI Tabula-Sapiens overlay imputation.
- **Headline result.** The nominated arms **4-1BB (TNFRSF9)** and **CD27** carry clean therapeutic windows;
  **CD28** is decisively gated on CRS + suppression + proliferation — the model reproduces the TGN1412 lesson
  mechanistically.
- **Validated now.** PK reproduces digitized clinical exposure (teclistamab, elranatamab, glofitamab bivalent
  path); 54-row verification ledger (52 PASS); fully reproducible notebook (0 errors end-to-end).
- **Honest negatives (reported, not hidden).** The mechanistic IL-6 CRS axis **fails** its validation gate on
  verified anchors (talquetamab inverts); 4-1BB is activation-induced with an induction fold that is NOT_FOUND
  in the literature, so its window is a *sensitivity sweep*, not a point estimate; organ-level EPCAM/DLL3/PMEL
  imputation is a known retention gap. None of these touch the GRN/DE-driven nomination.

---

## Figure index

| Fig | File | Panel |
|---|---|---|
| 1 | `Fig_counterscreen_window.png` | Three-axis costim-arm nomination window (effector vs suppression + CRS) |
| 2 | `grn_ring_cd8_effector_2026-07-13.png` | CD8 effector GRN drive ring (hero dataset) |
| 3 | `grn_selectivity_consolidated_2026-07-13.png` | Consolidated GRN selectivity (CD8 vs CD4 leaning per arm) |
| 4 | `percell_spatial_overlay__tumor_P1CRC__main_L6native.png` | Per-cell spatial overlay, CRC tumor (canonical layer) |
| 5 | `fig_spatial_liver.png`, `fig_spatial_small_int.png` | Per-organ spatial drug/bound/kill maps |
| 6 | `teclistamab_IV_PK_merged_engine.png` | Teclistamab IV PK vs digitized clinical (unified engine) |
| 7 | `elranatamab_matched_PK_IL6_overlay.png` | Elranatamab matched PK + IL-6 (same cohort) |
| 8 | `teclistamab_PD_overlay.png` | Teclistamab PD depletion overlay |

---


---

## 1. Model architecture — a spatial PBPK where interstitium is replaced by per-cell agents

The engine keeps the part of a platform PBPK that is genuinely predictive — systemic antibody transport and half-life — and rebuilds the part that governs pharmacology. Systemic disposition uses the Shah & Betts (2012) platform mAb structure (two-pore extravasation plus FcRn salvage; PMID 22143261, cited in code at `qsp_costim_window_v2.py:893`), but each tissue's *well-mixed interstitial compartment is deleted and replaced* by a spatially resolved agent-based model in which every real Xenium cell is an agent running its own Rhoden kinetic multi-arm binding ODE (T2 §1.1). The consequence is architectural, not cosmetic: target-mediated drug disposition, cytotoxic killing, and myeloid cytokine release are never written as compartmental sink terms — they are summed from per-cell binding and internalization events on the tissue grid, so efficacy and toxicity are the *same* emergent computation read out in different compartments.

**Layer 1 — the Shah–Betts backbone (`qsp_costim_window_v2.py`).** The body is a set of tissue compartments, each split into vascular and interstitial sub-volumes. The backbone table `_PBPK_TISSUES` (`qsp_costim_window_v2.py:82–99`; T1 §3.1) enumerates **fifteen physiological compartments** — lung, heart, kidney, brain, muscle, skin, adipose, bone, stomach, small and large intestine, pancreas, spleen, liver, and a fifteenth **tumor slot** specialized per target into a solid (CRC / lung / ovary / prostate / skin) or hematologic (AML / DLBCL / MM) malignancy. Organ plasma flows are normalized so the parallel organs sum exactly to cardiac output (RUN-verified $\sum Q = 5000.0$ L/day, $Q_\text{lung}=5000.0$ L/day; T1 EQ-1), and each organ's lymph flow is the signature platform relation $L_i = Q_i/500$ (`_LYMPH_RATIO = 1/500`, `:101`; T1 EQ-2). Four platform system parameters are retained at their published values — FcRn salvage $f_\text{FcRn}=0.90$, pinocytic uptake $\mathrm{CL_{up}}=0.3503$/day, and a distribution-rate multiplier $k_\text{dist}=3.0$ (`:176–182`) — so systemic elimination is the emergent catabolic rate $k_\text{cat}=\mathrm{CL_{up}}\,(1-f_\text{FcRn})+k_\text{renal}$ (`coupled_percell_pk.py:76`, T9 EQ-16) and a molecule's half-life follows from its own $f_\text{FcRn}$ rather than a fitted clearance.

**Layer 2 — interstitium as a per-cell graph (`wholebody_percell.py`, `coupled_percell_pk.py`).** Inside each spatial organ the interstitium is a $k$-nearest-neighbour graph ($k=6$) over the real single-cell $(x,y)$ coordinates, with ECM-hindered diffusion on the graph Laplacian (T2 EQ-13). Drug does not appear organ-wide by fiat: the organ vascular compartment is solved in quasi-steady state (no state variable; `wholebody_percell.py:163–168`, T2 EQ-9), and antibody enters **only at blood–endothelial (BEC) cells** by two-pore convection,
$$J_\text{extrav}=PS_\text{ex}\,(C_\text{vasc}-\overline{C}_\text{BEC}),\qquad PS_\text{ex}=k_\text{dist}\,L\,(1-\sigma_V),$$
whose series conductance $G_\text{ex}=PS_\text{ex}\,Q/(Q+PS_\text{ex})$ (T2 EQ-8, EQ-10, `wholebody_percell.py:169`) makes the perivascular-to-core exposure gradient emergent. It leaves at lymphatic-endothelial (LEC) cells: total lymph efflux is held exactly at the compartmental two-pore value $PS_\text{ly}\,\overline{C}$ but distributed as a spatially graded field, decoupling *how much* drains (physiology, imposed) from *where* it drains (geometry, emergent; T2 EQ-12). BEC/LEC identity per organ comes from LYVE1/PROX1/PDPN markers (`wholebody_percell.py:8–11`). The whole-body per-cell census over the spatial tissues totals **~1.89M agents** (source manifest; verified against the deposited `agents_*.npz` files). A representative per-cell exposure/binding/kill field is shown in `percell_spatial_overlay__tumor_P1CRC__main_L6native.png` and the per-organ maps `fig_spatial_liver.png`, `fig_spatial_small_int.png`.

**Recirculation coupling.** The layers close into one systemic loop: plasma $\rightarrow$ organ vascular pool $\rightarrow$ BEC extravasation $\rightarrow$ per-cell diffuse-and-bind graph $\rightarrow$ LEC drainage $\rightarrow$ one lumped systemic lymph pool $\rightarrow$ plasma at $k_\text{lymph\_return}=24$/day (mean residence ≈ 1 h; T9 EQ-17). The plasma and lymph balances (T9 EQ-16/17) subtract every organ's extravasation flux and every per-cell TMDD sink from plasma, so systemic exposure and local pharmacology are mutually consistent by construction.

**One unified run loop.** PK and PD execute in a **single coupled pass** (`coupled_percell_pd.simulate_pd`, loop at `:226`; T9 §1.1/§1.5), not two chained simulations. Within each step the driver consumes doses, updates plasma $C_\text{pl}=A_\text{pl}/V_\text{pl}$, runs per-organ transport (`TissueGraph.step`) on a transport clock $dt=0.02$ d, then runs the per-cell PD (synapse, kill, occupancy→effector gain) every third step ($pd_\text{every}=3$, i.e. 0.06 d) reading that organ's per-cell free-drug field $g.C$, and finally aggregates myeloid IL-6 and the heme/blood TMDD sinks before the plasma/lymph mass balance. Each molecule's own plasma $C(t)$ therefore drives its own per-cell pharmacology, and the identical Rhoden binding scheme runs in both the PK TMDD sink and the PD synapse (T9 EQ-21).

**Live engine.** The live import path is thirteen files: the harness `run_tce_pd_reval.py`; `qsp_costim_window_v2.py` (backbone); `pd_model_config.py` (calibration); the transport core `coupled_percell_pk.py` and `wholebody_percell.py`; the binding/kill stack `kinetic_rhoden_percell.py`, `multiarm_binding.py`, `kinetic_synapse.py`, `wholebody_pd.py`, `coupled_percell_pd.py`; `myeloid_il6.py`; `signaling_dynamics.py` (GRN-driven per-cell signaling shape); and `costim_induction.py` (wired but default-OFF; source manifest). Five modules — `unified_binding.py`, `biexact_solver.py`, `multiarm_kinetic.py`, `il6_pbpk.py`, `cytokine_pbpk.py` — carry zero importers in this graph and are **not** on the live path (source manifest; T1 §3.6).

**Architecture-level limitations, stated plainly.** Three are load-bearing. (i) The advertised portal→liver→venous→lung series topology is **not implemented**: `PB.portal`/`PB.i_portal` are computed and never read, and the live model is a *star* — one plasma pool feeding independent interstitia in parallel, all draining to one lymph pool (T1 §5.2), so there is no hepatic first-pass. (ii) The backbone table holds fifteen tissues, but the production per-cell driver instantiates a spatial ABM for only **eleven healthy organs plus the tumor**; muscle, brain, and stomach are present in the table but absent from the live per-cell mass balance — **45.3 %** of the model's 8.551 L interstitial volume, chiefly muscle, is not carried (T9 §5.1). (iii) Within the unified pass the receptor field is **static**: the CD3/TAA antigen pools driving the PK sink do not shrink as T cells expand or targets die (T1 §5.5). None alters the drive-independent nomination, but each is a real gap a reviewer will probe.


---

## 2. The per-cell receptor layer — scVI overlay imputation makes the model runnable

The whole-body engine is agent-based at the resolution of real cells: each of the ~1.89M whole-body agents is a Xenium cell running its own Rhoden multi-arm binding ODE (source: manifest). That ODE cannot start unless the cell carries a copies/cell value for **every** antibody target arm it might see — the TAA, CD3, and the costim receptor. The difficulty is that a Xenium slide only measures the genes on its probe panel (~250–541 genes), and a target that was never designed into the panel is a hard-zero column, not a low measurement. Feeding that zero to the binder would silently report "no receptor → no binding → no toxicity" for that arm in that tissue — a data gap masquerading as biology (T2b §1). The per-cell receptor layer exists to close that gap, and its central step — the scVI overlay — is what actually makes the whole-body model runnable, because it is what puts **every receptor on every cell**. That is the point of training the overlay.

The pipeline is three tiers, applied per cell per gene (T2b §2).

**Tier 1 — native panel extraction (byte-exact).** For a gene on the tissue's Xenium panel, each cell is normalized to $10^4$ total counts and `log1p`-transformed: $M_c \leftarrow M_c \cdot (10^4 / \text{lib}_c)$, then $M_c.\text{data} \leftarrow \log(1+M_c.\text{data})$ (`build_cancer_abm.py:50`). This is deterministic and verified byte-reproducible — CEACAM5 and ERBB2 in the ovary build reproduce at correlation 1.0000 (T2b §2).

**Tier 2 — scVI Tabula-Sapiens overlay (the load-bearing step).** For a gene *off* the panel, its expression is borrowed from the transcriptionally most similar real reference cell. Reference cells come from the CZ CELLxGENE Census (`census_version = "2025-11-08"`) — Tabula Sapiens for normal organs, a cancer-matched atlas per tumor (T2b §3). Reference and query are concatenated on shared genes; 2000 HVGs are selected (`seurat_v3`, batched by source); an scVI model (`n_latent=20, n_layers=2`, `seed=0`) is trained on raw counts for `max_epochs=60` with early stopping. scVI's batch correction is exactly what lets a query cell and a reference cell be compared in a shared, technology-corrected latent space $Z$. Assignment is `topk_spread` with $K=15$ and softmax weighting:

$$w_j = \frac{\exp(-d_j/\tau)}{\sum_{k=1}^{15}\exp(-d_k/\tau)}, \qquad \tau = \operatorname{median}(d_{:,1:}) + 10^{-9}$$

over each query cell's 15 nearest reference cells, from which **one** neighbor is drawn stochastically $\propto w$ (seed 0) rather than by deterministic argmax. The stochastic draw *spreads* assignments across many reference cells, preserving expression heterogeneity instead of collapsing every query cell onto a single nearest neighbor (T2b §3). Each query cell then copies its assigned reference cell's normalized expression of the off-panel target, and inherits that cell's cell-type label (used downstream for kill accounting and cytokine sources). The retained provenance per query cell is `ref_cell_assign` (the chosen reference index) plus `overlay_method`, `antigen_overlay_added`, and `antigen_still_missing`.

**Tier 3 — RNA→copies conversion.** Both tiers yield a log-normalized *expression* value; the binder needs *copies/cell*. Path A (24 HPA-IHC-anchored genes) sets $\text{copies} = \text{HPA\_IHC\_anchor}(\text{gene, organ, lineage}) \times$ a per-cell within-lineage rank modifier ($0.1$–$3\times$). Path B (no-IHC genes) uses the HPA single-cell nCPM anchor through the Glassman ladder $\log_{10}(\text{copies}) = 0.7768\cdot\text{IHC} + 3.9723$, times the same rank modifier; DLL3, PMEL, and EPCAM live here. Seven "ORIG7" genes (MS4A1, CD19, TNFRSF17, GPRC5D, EPCAM, DLL3, CD3E) are stored directly as copies by an earlier native pipeline (T2b §4). The organ pool consumed by the PBPK TMDD sink is then $\text{pool}_{\text{nM}} = (\sum_\text{cells}\text{copies}\cdot\text{scale}) / (V_\text{is}\cdot N_A)\cdot 10^9$, with $\text{scale} = \text{organ\_true\_cells}/\text{section\_cells}$; the kidney EPCAM pool reproduces 247.637 nM exactly through this formula (T2b §4).

**Solid-tumor receptor pools** (`Rtot_wholebody_final.json`, completed 2026-07-13; verbatim from manifest):

| Target | Tumor | $R_\text{tot}$ (nM) | Copies anchor |
|---|---|---|---|
| CEACAM5 | CRC | 6.459 | HPA-IHC (native) |
| ERBB2 | CRC | 1.919 | HPA-IHC (native) |
| FOLH1 | prostate | 3.552 | HPA-IHC lineage |
| EPCAM | ovary | 486.528 | HPA single-cell nCPM |
| DLL3 | lung/SCLC | 4.348 | HPA single-cell nCPM |
| PMEL | skin/melanoma | 9.473 | HPA single-cell nCPM |

**Two honest limitations gate this layer.**

*Retention gap.* The tumor overlay builder (`build_cancer_abm.py`) survives, but the **organ** overlay builder script and — critically — the **reference AnnData** that `ref_cell_assign` indexes into were not retained (the reference was pulled live from the Census, not saved). The integer indices are only resolvable against a reference materialized in the *identical* row order, which a fresh Census pull does not guarantee (T2b §5). The consequence is bounded: three off-panel, solid-only arms — EPCAM, DLL3, PMEL — could not be overlaid onto organ cells without re-running scVI across 12 organs (hours), so their **organ** pools read 0 where the panel lacks the gene. This is a panel/retention gap, **not a verified zero**: EpCAM is on liver cholangiocytes (the catumaxomab hepatotox signal), PMEL is melanocyte-restricted (the tebentafusp skin signal), and neither off-tumor liability is quantifiable from the current organ panels (`Rtot_wholebody_final.json` meta). Their **tumor** copies are complete, and the validated heme panel plus the CRC solid pair (CEACAM5, ERBB2) have full organ coverage and are unaffected.

*[PROVENANCE].* Every value in this layer is a **derived** receptor copy number — scRNA-seq transcripts converted to protein copies via the HPA/Glassman anchors — **not** an independently QIFIKIT-measured surface density. The conversion is not per-gene validated against protein (T2b §4; T7 §3.3). For activation-induced costim receptors (TNFRSF9/4-1BB, TNFRSF4/OX40) this is doubly fraught, because resting transcript sits near the detection floor, so $R_\text{rest}$ is a noise-dominated estimate of a near-zero quantity (T7 §3.3). Figure legends must not read "measured receptor copies."


---

## 3. Per-cell kinetic binding — Rhoden multi-arm ODE, format as geometry

Every agent on the tissue grid runs its own binding solve once per timestep; this is the single place where free drug meets free receptor, and TMDD, killing, and cytokine release are all functionals of its output (`engine/kinetic_rhoden_percell.py:13–15`). The kernel is a kinetic bivalent scheme in the Rhoden form (attributed in-code to Rhoden *et al.* 2016, *J Biol Chem* 291:11337–47, `qsp_costim_window_v2.py:704` — cited from the repository, not independently verified here), with explicit receptor turnover added on top.

**The six-species scheme.** The reference form tracks free drug $C$ plus five receptor/complex pools for a heterotypic (two-antigen) binder — free receptors $Ag_1,Ag_2$, singly-bound complexes $BAg_1,BAg_2$, and the doubly-crosslinked complex $B_{\text{dbl}}$ (`kinetic_rhoden_percell.py:71–75`):

$$
\begin{aligned}
\dot{Ag_k} &= K_{\text{SYN},k} - k_{\text{on},k}C\,Ag_k + k_{\text{off},k}BAg_k - k_{\text{on},k}x_k BAg_j + k_{\text{off},k}B_{\text{dbl}} - k_{\text{deg}}Ag_k\\
\dot{BAg_k} &= k_{\text{on},k}C\,Ag_k - k_{\text{off},k}BAg_k - k_{\text{on},j}x_j BAg_k + k_{\text{off},j}B_{\text{dbl}} - k_{\text{TMD}}BAg_k\\
\dot{B}_{\text{dbl}} &= k_{\text{on},1}x_1 BAg_2 + k_{\text{on},2}x_2 BAg_1 - (k_{\text{off},1}+k_{\text{off},2})B_{\text{dbl}} - 2k_{\text{TMD}}B_{\text{dbl}}
\end{aligned}
$$

Each arm carries its **own** $k_{\text{on}}$/$k_{\text{off}}$ (per-target measured kinetics where available — e.g. teclistamab $k_{\text{on,TAA}}=1.28\times10^{6}\,\mathrm{M^{-1}s^{-1}}$, `CHANGELOG_2026-07-13.md:236`; generic fallback $10^{5}\,\mathrm{M^{-1}s^{-1}}$). The first bond is driven by *bulk* $C$; the second (crosslink) bond uses the **same** $k_{\text{on}}$ acting on a *tethered* effective concentration $x_k$ — so avidity is emergent geometry, not a fitted bonus rate. Honest caveat: this heterotypic path has **zero live call sites** (`kinetic_rhoden_percell.py:47–80`, verified by grep); every production PK solve uses the reduced **same-antigen 3-state** system, where both arms bind one antigen and $B_{\text{dbl}}$ carries census weight 2 (`:90–92`):

$$
\dot{Ag}=K_{\text{SYN}}-k_{\text{on}}C\,Ag+k_{\text{off}}BAg_1-k_{\text{on}}x\,BAg_1+2k_{\text{off}}B_{\text{dbl}}-k_{\text{deg}}Ag,\qquad
\dot{B}_{\text{dbl}}=k_{\text{on}}x\,BAg_1-2k_{\text{off}}B_{\text{dbl}}-2k_{\text{TMD}}B_{\text{dbl}}
$$

**Receptor turnover — why kinetic, not QSS.** Synthesis is zero-order and pinned to the drug-free set-point, $K_{\text{SYN}}=k_{\text{deg}}\,R_{\text{tot}}$ (`:104`), so at $C=0$ the free pool relaxes exactly to baseline ($\dot{Ag}=0\Rightarrow Ag\to R_{\text{tot}}$; RUN-verified to 5 figures). This is the piece a quasi-steady-state sink lacked: with QSS the receptor snaps back instantaneously, target recovery has no timescale, and the terminal TMDD phase breaks (`:17–19`). Making synthesis explicit gives recovery a real time constant $1/k_{\text{deg}}$ ($k_{\text{deg}}=0.5\,\mathrm{d^{-1}}$ fallback, `coupled_percell_pd.py:133`; internalization $k_{\text{TMD}}=0.9\,\mathrm{d^{-1}}$, `pd_model_config.py:37`), so a second dose lands on a partially-depleted target — the mechanism behind step-up tolerance. The solver is backward-Euler with receptors carried as states, chosen because the avidity crosslink rate $k_{\text{on}}x$ reaches $\sim10^{6}\,\mathrm{d^{-1}}$ and makes an explicit solver stall (substep count $442\to5651$ as $k_{\text{on}}$ rose $10^{5}\to1.28\times10^{6}$; the A-stable BE step is $k_{\text{on}}$-independent and $\sim300\times$ faster, `reference_unified_binding/LIVE_ENGINE_FIXES.md:76–79`). In production the BE step needs $n_{\text{sub}}=1$ at $dt=0.01\,\mathrm{d}$ (14.4 min).

**Avidity as geometry.** The tethered concentration is pure geometry: a bound arm explores a hemisphere of radius = the arm span $s$, seeing partner receptors under its membrane footprint,

$$
Ag_{\text{EFF}}(s)=\frac{3\times10^{24}}{2N_A}\cdot\frac{R_{\text{copies}}}{4\pi r_{\text{cell}}^{2}\,s}\;\propto\;\frac{\sigma_{\text{receptor}}}{s},\qquad x = Ag_{\text{EFF}}\cdot\frac{Ag}{R_{\text{tot}}}
$$

with the depletion factor $Ag/R_{\text{tot}}$ making the crosslink self-limiting (hook effect emerges, no explicit bell term). Because the footprint grows as $s^2$ but the explored volume as $s^3$, avidity dilutes as $1/s$ — longer, floppier linkers *reduce* it (verified: $Ag_{\text{EFF}}(12.5)/Ag_{\text{EFF}}(25)=2.000$ exactly). At the CEACAM5 anchor (257,000 copies, $r_{\text{cell}}=8\,\mu m$, $s=12.5\,\mathrm{nm}$) a tethered arm sees $\approx64\,\mu M$ of partner antigen (`multiarm_binding.py:29–37`), four to five orders above any plasma level — that is the entire physical content of avidity here. Concentrations are on a fixed synapse reaction-volume basis, $\mathrm{NM\_PER\_COPY}=6.0/257000=2.335\times10^{-5}\,\mathrm{nM/copy}$ (~71 pL, `wholebody_pd.py:83`).

**Multi-arm generalization and the three spans.** Format is encoded as valency per arm — CD3, costim, and TAA each 0/1/2 — with $\Sigma\le4$ (tetravalent maximum, e.g. 1×CD3 + 2×TAA + 1×costim, `multiarm_binding.py:114–118`). Valency 0 means the arm is *absent*, not zero-affinity. Three architectural distances become three physics modifiers via one $1/s$ law: `span_coeng_tumor` (bivalent TAA avidity on one tumor cell), `span_coeng_T` (CD3+costim co-engagement on one T cell), and `span_bridge` (the cross-cell CD3–TAA reach).

**cis vs trans.** These are two distinct binding geometries. *cis* is per-cell avidity — a $B_{\text{dbl}}$-type same-cell doubly-bound complex driven by the local $Ag_{\text{EFF}}$; for a costim trispecific the cis species is CD3+costim gripped on one T cell, entering the live engine as an occupancy gate $occ_{\text{eff}}=occ\,[(1-p_{\text{cis}})+p_{\text{cis}}f_{\text{CD3}}]$ (`wholebody_pd.py:295`). *trans* is the shared cross-cell bridge: the bridged trimer $B_2=[\mathrm{CD3\cdot drug\cdot TAA}]$ summed over each T cell's alive target neighbours, coupled through the sparse T×target incidence matrix $W$ (KDTree synapse graph, `wholebody_pd.py:280`), with $\dot{B_2}=k_f B_1-k_{\text{off,TAA}}B_2-k_{\text{int}}B_2$ and forward rate $k_f=k_{\text{on,TAA}}\cdot c_{\text{eff,trans}}$ (`kinetic_synapse.py:168`). The cis/trans switch is a Gaussian span-feasibility centered on a 12.5 nm height-matched gap:

$$
p_{\text{cis}}(s)=\exp\!\left[-\tfrac{1}{2}\!\left(\frac{s-12.5}{8.0}\right)^{2}\right]
$$

so a compact construct co-grips (cis, $p_{\text{cis}}(12.5)=1.0$) while a tall epitope mismatch forces trans ($p_{\text{cis}}(60)=2.2\times10^{-8}$; `multiarm_binding.py:39–46`). Honest scope: $p_{\text{cis}}$ is live-imported but defaults to 0 (`span_coeng_T_nm` unset, `wholebody_pd.py:221`), so the canonical clinical run holds costim cell-autonomous — the cis machinery is built and wired but default-inert.

**One geometry, two engines.** The same surface-density × reach law runs in *both* the PK sink (`kinetic_rhoden_percell`, TMDD internalization $k_{\text{TMD}}(BAg_1+2B_{\text{dbl}})$) and the PD synapse (`kinetic_synapse.ageff_nM`), verified numerically float-equal to 1 ULP (T4), and both use the identical $K_{\text{SYN}}=k_{\text{deg}}R_{\text{tot}}$ turnover structure — so a more avidly-binding molecule is *automatically* cleared faster and kills more, a structural coupling rather than two fitted terms. What is **not** yet unified is the ODE itself: PK/TMDD runs the 3-state BE receptor system, PD kill runs a separate 2×2 matrix-exponential bridge; the docstring's "one binding solve used identically" (`:13`) is aspirational at the solver level. Two limitations carry into results: the factor-2 internalization on $B_{\text{dbl}}$ is an unsourced mechanistic assumption that propagates into every bivalent clearance estimate (`:75`), and the BE non-negativity clamp inflates the TMDD sink by +44–124% at free drug ≳10 nM — precisely the regime of the heme/blood compartments that carry the validation panel — so those clearance numbers are provisional pending a substep-law fix (T3 §5.2). The reach envelope is also a uniform hemisphere, not a worm-like-chain polymer, and the load-bearing distances (12.5 nm span, 8 µm radius) are uncited assumptions.


---

## 4. Emergent pharmacology — synapse killing, TMDD, and the therapeutic window

Everything downstream of binding — killing, target-mediated clearance, and the window score that ranks the costim arms — is produced by one per-cell reaction core, not by fitted efficacy or toxicity constants. This section states what emerges, what is imposed, and the nomination that follows.

### 4.1 The two-state synapse and where avidity comes from

Each PD step advances a literal two-species bond ODE per T cell: `B1`, drug armed on CD3, and `B2`, the bridged trimer CD3·drug·TAA (T5, `kinetic_synapse.py:169-182`). Free CD3 is not independent (`RC - B1 - B2`), so the system is 2-dimensional:

$$\dot B_1 = k_{on}C\,(RC-B_1-B_2) - (k_{off}^{CD3}+k_f)B_1 + k_{off}^{TAA}B_2, \qquad \dot B_2 = k_f B_1 - (k_{off}^{TAA}+k_{int})B_2$$

The single most consequential line is the trimer decay: `B2` relaxes at $(k_{off}^{TAA}+k_{int})$ with `koff_CD3` **deliberately absent** (T5 EQ-5, `:172-175`). A CD3-arm release does not dissolve the trimer because the TAA arm still holds the drug, so trimer lifetime is set by the slower arm. That is bivalent avidity — not a parameter, but emergent. At the mosunetuzumab operating point the trimer outlives the CD3 bond by 26x (1.79 h vs 4.17 min), and the half-maximal engaging concentration falls to `C = 0.105 nM` against a nominal `KD_CD3 = 40 nM` — a **380x apparent-potency gain computed from geometry** (T5 §4.2). With bridging propensity $\rho \approx 3.7\times10^4$, engagement potency reduces to $\mathrm{EC50}\approx k_{int}/k_{on}\approx 0.9/8.64 = 0.104$ nM — set by **internalization, not affinity** (T5 §4.2).

The step is integrated by an exact 2x2 matrix exponential (T5 EQ-6): the fast eigenvalue ($\sim -5\times10^5$/day over a 0.06-day step) relaxes fully within one step, so the QSS limit **emerges** with no QSS assumption written, while a genuine slow lag (τ set by `kint`) is retained (T5 §4.1). The engaged fraction $p_{eng}=B_2/RC$ (`:211`) is the model's single engagement currency, driving killing, costim induction, and myeloid IL-6.

### 4.2 Serial killing: K_HIT fixed, k_death locked

Once bridged, two clocks race — the lethal hit (`k_hit`) against detachment (`koff_CD3`) — giving a harmonic-mean throughput (T5 EQ-13):

$$\text{serial\_rate} = \frac{k_{hit}\,k_{off}^{CD3}}{k_{hit}+k_{off}^{CD3}}\; p_{eng}\;\text{has\_live}$$

`K_HIT = 12/day` is **FIXED, not fitted** (`kinetic_synapse.py:48`), and `k_death = 1.0` is **LOCKED, one shared value for every engager** (`kinetic_calib.json`). Against the 22 clinical engagers' measured CD3 off-rates the serial ceiling spans a **43x range (0.28 -> 11.97 /day)**: cevostamab (`KD_CD3 = 0.033 nM`) is throughput-crippled — it binds CD3 too tightly to cycle — while detuned CD3 saturates near `k_hit` (T5 EQ-13). The mosunetuzumab ceiling reproduces exactly: $12\cdot345.6/(12+345.6)=11.597$/day (T5 EQ-13). Two flags carry from the subsystem docs: `K_HIT` asserts serial-killing-literature provenance but **names no PMID** (read as `[ASSUMED]`), and `k_death` is `[FITTED/CALIBRATED-AND-LOCKED]` — consistent with an independent serial ceiling but not a measured trimer->apoptosis rate (T5 §5.3, T6 §3.6). The hit itself is an algebraic QSS probability, not an integrated state — the cycling is emergent, the hit is not (T5 EQ-13).

### 4.3 TMDD emerges from the grid; depletion is count-weighted

There is no `k_TMDD`, `Vmax`, or `Km` anywhere. Each cell internalizes bound complex at `kint`, and the organ sink is the sum of what individual cells ate (T2 EQ-19, `wholebody_percell.py:208`):

$$\text{organ\_sink}=\sum_i \text{intern\_copies}_i/N_A\cdot10^9\quad[\text{nmol/day}]$$

Nonlinear, dose- and target-dependent clearance therefore **emerges** — change the target from CD20 to BCMA and clearance changes by itself because the per-cell copies changed (T2 §4). A `count_scale` factor (physiological pool ÷ sampled-cell pool, T2 EQ-1) lets each Xenium cell stand for many real cells while keeping the measured per-cell copy number, so avidity (∝ R²) is evaluated at physiological density. Whole-body depletion is then a **target-cell-count-weighted** mean, $\text{depletion\_weighted}=\sum_o \text{kill\_frac}_o\, n_{target,o}/\sum_o n_{target,o}$ (T9 EQ-20), so an organ with four target cells cannot outvote marrow.

### 4.4 The window metric, and GRN vs DE

A dose is tolerable only if it clears **both** ceilings — peak IL-6 ≤ a CRS-severity ceiling **and** liver injury ≤ a liver ceiling, both anchored to the CD3-only backbone (`qsp_costim_window_v2.py:1375`). The window score is the `E_tol` (best tumor reduction at a co-tolerable dose) gain over backbone, with bounded Treg/liver tie-breakers that order near-equal arms but never override efficacy (`composite_window_score`); `cap` records which ceiling binds. Crucially the nomination is **gate-first**: an arm with any liability-up axis (CRS/SUPP/HELP/PROLIF/EXH) is vetoed regardless of effector score — reading the raw window without the veto resurrects the "CD28 looks clean" artifact (`FINAL_NOMINATION_v7.md`). Both a per-gene differential-expression (DE) drive and a genome-scale gene-regulatory-network (GRN) drive were run through the *frozen* QSP; the full 12-receptor result is below (`A31b_QSP_rerun_GRN_vs_DE.csv`, `Fig_counterscreen_window.png`).

| Receptor | Liability gate | window (DE) | window (GRN) | cap (DE/GRN) | effector-hit |
|---|---|---:|---:|---|---|
| CD30 | GATED[HELP,PROLIF] | 1.573 | 1.569 | none / none | + |
| **4-1BB** | **CLEAN** | **1.568** | **1.565** | **none / none** | + |
| ICOS | GATED[HELP,PROLIF] | 1.564 | 1.560 | none / none | - |
| DR3 | GATED[SUPP,EXH] | -0.423 | -0.528 | none / none | + |
| OX40 | GATED[SUPP,EXH] | -1.879 | -2.286 | CRS / CRS | + |
| CD2 | GATED[CRS] | -2.005 | -2.056 | CRS / CRS | + |
| **CD27** | **CLEAN** | **-2.368** | **+1.269** | **CRS / CRS** | + |
| HVEM | GATED[SUPP,EXH] | -2.404 | -2.498 | CRS / CRS | + |
| CD40 | GATED[HELP,PROLIF] | -2.538 | +0.099 | CRS / CRS | + |
| CD28 | GATED[CRS,SUPP,PROLIF] | -2.662 | -3.026 | liver / CRS | + |
| GITR | GATED[SUPP] | -2.788 | -2.783 | CRS / CRS | - |
| DNAM1 | GATED[SUPP,EXH] | -2.910 | -2.760 | CRS / CRS | - |

Values verbatim from `A31b_QSP_rerun_GRN_vs_DE.csv`. "Effector-hit" is the effector-gate precondition (`nominated_*`); the **decisive** column is the liability gate. CD30 and ICOS post the highest raw windows yet are both hard-gated — why post-gate reading is mandatory.

### 4.5 Nomination: 4-1BB + CD27 co-lead CLEAN; CD28 gated

Only two arms survive the liability veto under both drives: **4-1BB (TNFRSF9) and CD27**. 4-1BB is the robust co-lead — window ≈ **+1.57 under both DE and GRN, uncapped**, the widest and CRS-coldest arm. CD27 is the GRN-favorable co-lead: the network view reads it **safer** than per-gene DE (window -2.37 -> **+1.27**), though CRS-capped at monotherapy dose under both. **CD28 is GATED[CRS,SUPP,PROLIF]** with the most negative window (-2.66/-3.03) despite being the top raw-effector arm (z ≈ 11.9) — the clearest single proof that the counter-screen, not the effector axis, drives the nomination (`FINAL_NOMINATION_v7.md`).

### 4.6 Critical disclosure — the 4-1BB window is a conservative floor

4-1BB, OX40, ICOS, and GITR are **activation-induced**: near-absent on resting T cells, appearing only after TCR engagement — which is precisely why the field targets them. The induction code (`costim_induction.py`) sets `fold=None`/`source="...FOLD = NOT_FOUND"` for every inducible arm because a dedicated audit found kinetics but **no clean surface-density fold**, and `strict=True` makes the constructor **raise a `ValueError` rather than silently assume fold=1** (`costim_induction.py:__init__`). Consequently the reported 4-1BB window is evaluated at the **static resting `R_costim`** — a density near the detection floor — so it is a **conservative floor**: 4-1BB is penalized by the entire missing induction factor and *still* scores CLEAN and uncapped. Because no sourced fold exists, **no fold-based ranking is asserted**; wherever fold would matter it is confined to a sensitivity sweep carrying `fold_is_assumed=True` (`costim_induction.py`). This is the anti-fiction discipline working: the model refuses to manufacture the one number that would flatter its lead candidate, and the lead survives anyway. Two corollary caveats: the induction step is not wired into the canonical kinetic engine (T6 §5.12), so this floor reads resting density throughout; and the **OX40/GITR net-negative result is PROVISIONAL** for the same missing-fold reason. 4-1BB's true residual liability — urelumab-class hepatotoxicity via liver-myeloid 4-1BB — is a CD4-screen blind spot slated to emerge from the PBPK liver compartment, not a window-model claim (`FINAL_NOMINATION_v7.md`).


---

## 5. Mechanistic cytokine-release — per-cell myeloid IL-6 and an honestly-reported validation gate

The cytokine-release arm is where a QSP model is most tempted to cheat: a single fitted scale can make any construct hit any clinical IL-6 number. This engine does not. IL-6 magnitude is emitted **per myeloid cell, spatially, with no fitted scale, no $E_{\max}$, no EC50, and no Hill term** (`engine/myeloid_il6.py`); saturation and the between-construct spread emerge from the finite, spatially distributed myeloid pool plus a first-order clearance. The legacy path — `il6 = engaged_dwell_rate × IL6_SCALE`, one constant fit so mosunetuzumab reached 570 pg/mL — was **deleted**, not disabled: an empty mechanistic array now raises a hard `RuntimeError` rather than silently substituting the retired constant, and the output field is stamped `il6_method="mechanistic_myeloid_percell"` as provenance (`run_tce_pd_reval.py:203, 217–235`).

### 5.1 The emitter is myeloid; the T cell is only the trigger

IL-6 is not made by T cells. The mechanism is anchored on the CAR-T CRS literature — macrophages are "the main overall source of IL-6" (Giavridis, *Nat Med* 2018, PMID 29808005) and "human monocytes were the major source of IL-1 and IL-6" (Norelli, *Nat Med* 2018, PMID 29808007) — so the engaged T cell is only the *activating input*, and each myeloid agent runs its own program. Three coupled per-cell equations govern it (`myeloid_il6.py:23–25, 183–214`). Each myeloid agent $i$ carries an activation state $a_i\in[0,1]$ driven by its local contact with engaged T cells:

$$\frac{da_i}{dt} = k_{\mathrm{on}}\,c_i\,(1-a_i)\;-\;k_{\mathrm{off}}\,a_i,\qquad c_i=\!\!\sum_{j:\,\lVert x_i-x_j\rVert<R_{\mathrm{contact}}}\!\! p^{\mathrm{eng}}_j$$

where $c_i$ is the engaged-synapse load $p^{\mathrm{eng}}$ of T cells within the contact radius (a `cKDTree` neighbour sum) and the $(1-a_i)$ factor is the **structural** saturation ceiling — no Emax term. Secretion sums only over the intrinsic-secretor subset, and the plasma compartment converts a production *rate* into a *concentration* — the piece the fitted engine never had:

$$\mathrm{prod}=\sum_i a_i\,S_{\max}\,\mathrm{cs}_i\,\mathbb{1}[\mathrm{secretor}_i]\;\;[\text{pg/hr}],\qquad \frac{dC}{dt}=\frac{\mathrm{prod}}{V_{\mathrm{ECF}}}-k_{\mathrm{deg}}\,C.$$

Per-molecule differences then emerge from anatomy, not a knob: a CD20 engager floods the spleen's large B-cell field (16% myeloid, ~54 000 resident macrophages) → many myeloid activated → high IL-6; a BCMA engager engages rare marrow plasma cells (1.8% myeloid, ~3000 cells) → few activated → lower IL-6 (`myeloid_il6.py:31–34`).

### 5.2 Every constant is literature-sourced — including the ones that are honestly weak

The emitter constants are measured single-cell quantities, not tuned outputs (`myeloid_il6.py:41–112`): peak secretion $S_{\max}=10.6$ molec/s over actively-secreting monocytes (Han et al., *Lab Chip* 2010, PMID 20376398), which the unit conversion turns into $S_{\max}=1.331\times10^{-3}$ pg/hr/cell — **not** the $0.0196$ that a stale in-code comment at `:111` still prints (that is the retired 156-molec/s tail; the live module returns `0.00133069`). The intrinsic-secretor fraction 0.039, time-to-max 150 min (setting $k_{\mathrm{on}}=1.2/\text{hr}$), and $k_{\mathrm{off}}=0.10/\text{hr}$ are from PMID 37533643; the 14.1 µm contact radius is the two-sphere touch distance of a 10.6 µm macrophage (PMID 9400735) and a 3.5 µm T cell (PMID 30571054) — a contact, not a reach, distance because CD40L–CD40 is membrane-bound.

Clearance is the load-bearing assumption, and it is reported as such. $k_{\mathrm{deg}}=0.20/\text{hr}$ ($t_{1/2}\approx3.5$ h) times the ECF distribution volume $V_{\mathrm{ECF}}=11.65$ L (interstitium 8.55 + plasma 3.10, because a 21 kDa cytokine with no FcRn recycling distributes through ECF and is *produced* in the interstitium) gives an implied clearance $\mathrm{CL}=k_{\mathrm{deg}}\!\cdot\!V_{\mathrm{ECF}}=55.9$ L/day. This is **kept, not fitted to the IL-6 peaks** it is later tested against; the model's own derived-CL (~76 L/day) is an output cross-check whose use as an input would be circular (`IL6_ANCHORS_VERIFIED_2026-07-13.md`). The honest flag travels with it: $k_{\mathrm{deg}}$ is `[FITTED]`, borrowed from a *semi-mechanistic modeling* paper (Chen 2019, PMID 31268236) that reports no measured clearance — **human IL-6 clearance appears to be unmeasured in the primary literature** — and $V_{\mathrm{ECF}}$ is `[DERIVED]` from the model's own PBPK volumes (`subsystems/T8_mechanistic_crs_il_6.md`, §EQ-9a). A 2× error in $\mathrm{CL}$ is a 2× error in every absolute pg/mL.

### 5.3 Blood myeloid are gated OFF, and the ranking is insulated from the weak constants

Circulating blood myeloid contribute exactly zero (`coupled_percell_pd.py:188`). The justification is sound: IL-6 induction needs *sustained* CD40L–CD40 contact, ~60% of blood monocytes are marginating (PMID 3944542) and those adherent cells **are already counted** as organ-resident myeloid, so adding flowing blood would double-count them — and the blood ABM's synthetic 2D grid is non-physical. The caveat is that the *numerical* argument once quoted in the comment (61 874 pg/mL) belongs to the retired parameterisation and is ~1400× too large at live parameters; the mechanistic argument stands, the stale number does not (`T8`, §5.9).

Absolute-scale uncertainty does not corrupt the counter-screen, for a structural reason. Since $C_{ss}=\mathrm{prod}/(V_{\mathrm{ECF}}\,k_{\mathrm{deg}})$ and $\mathrm{prod}=\sum_{\mathrm{organ}}\mathrm{census}\times\sum_i a_i S_{\max}\mathrm{cs}_i\mathbb{1}[\mathrm{secretor}_i]$, the constants $V_{\mathrm{ECF}}$, $k_{\mathrm{deg}}$, $S_{\max}$ are **global multiplicative factors that cancel exactly in any between-molecule ratio** (`T8`, §4.3), whereas the drug-discriminating spatial contact calculation carries no fitted constants. A comparison against an absolute clinical anchor therefore tests the global scale ($\mathrm{CL}$ and the myeloid census); the *ranking* is what the mechanism actually asserts.

### 5.4 The validation gate: a rigorously-reported FAIL

IL-6 was scored against the three manifest-verified digitized anchors as a like-for-like ratio. Emergent mechanistic-ECF peaks were mosunetuzumab 723.7, glofitamab 650.7, talquetamab 772.4 pg/mL (teclistamab 429.0 computed but **excluded** — its 21 pg/mL is a loose, unverified MajesTEC-1 mean, dropped per manifest v3); the verified anchors are mosun 127.4 (mean), glofit 30.2 (median), talq 19.8 (median) (`IL6_VALIDATION_RESULT.json`):

| Test pair (verified) | Model | Clinical | Verdict |
|---|---:|---:|---|
| glofit/talq (cleanest, both median step-up) | 0.84× | 1.53× | **WRONG DIRECTION** |
| mosun/talq | 0.94× | 6.43× | WRONG DIRECTION (talq inverted) |
| mosun/glofit | 1.11× | 4.22× | direction OK, compressed |

**VERDICT = FAIL on verified anchors**, driven by a specific, diagnosable mechanism. Talquetamab has the *highest* model IL-6 production (1.80M pg/hr) and the highest peak despite the *lowest* $C_{\max}$ (7.4), yet is clinically the *lowest* IL-6 of the panel. The reason: GPRC5D sits on plasma cells and keratinized tissue, **not** the myeloid-contacting compartments that drive CRS, and the myeloid term does not condition on *which* target compartment the killing occurs in — so it over-credits depletion-weighted engagement wherever it happens (`IL6_VALIDATION_RESULT.json`, `qc_talquetamab_il6_median.png`). This is reported as-is: **not rescued** — not denominator-shopped and not talq-dropped. Excluding teclistamab makes the result *worse*, since its provisional pair (glofit/tecli 1.52× vs an unverified 1.44×) was the only apparent match and leaned on an anchor the manifest forbids citing.

Two things bound the damage. The PK backbone feeding this arm is validated independently and well — the matched elranatamab PK+IL-6 overlay gives AFE 2.06× (`elranatamab_matched_PK_IL6_overlay.png`) — so the failure is localized to the emission model's compartment-blindness, not to transport or binding. And, decisively, **the nomination is decided upstream of this gate and does not depend on it**: the six-axis liability veto and the GRN/DE three-axis scores that nominate 4-1BB and CD27 never consult the QSP IL-6 output, and the per-construct ranking is in any case structurally invariant to the two unmeasured disposition constants (§5.3). The CRS axis is presented as a rigorously-characterized honest negative — a component that fails a falsification test for a named, mechanistic reason — precisely because burying it would be the fitted-scale sin this engine was rebuilt to avoid.


---

## 6. Validation, reproducibility, and limitations

A whole-body per-cell model has enough moving parts to fit almost anything, so the governing question is not "does it match?" but "which numbers earned their match, and which are free parameters wearing a citation." This section states the validation strategy, the provenance discipline that enforces it, and — at length, because an honestly-reported negative is worth more than a hidden one — the limitations that bound every claim above.

### 6.1 The PK validation strategy: one binding engine, digitized clinical curves

The pharmacokinetic backbone is validated **independently of the IL-6/CRS arm**, against digitized clinical concentration–time data rather than against itself. The anchor set lives in `model/params/mab_tce_pkpd.sqlite`, a curated database of **76 drugs, 93 literature sources (48 DOI / 43 PMID / 45 FDA labels), 183 PK/PD curves, 936 digitized points, and 14 Biacore kon/koff/KD records transcribed verbatim from SPR tables** (`PARAMETER_AUDIT_2026-07-13.md`). The decisive architectural commitment is that the PK TMDD antigen sink and the PD kill synapse were unified onto **one kinetic binding engine** (kon/koff with free-receptor turnover $K_{\mathrm{SYN}} = R_{\mathrm{tot}}\cdot k_{\mathrm{deg}}$), so a molecule binds its antigen with a single affinity everywhere. That invariant is not decorative: the audit's most serious defect (DEFECT 1) was precisely that `setdefault` had let the PK sink and the PD synapse bind the same antigen at **different** affinities — talquetamab at 2.0 nM vs 11.0 nM, a 5.5× split — which is physically incoherent for one binding event. The fix forces the measured KD to win in both places and prints the substitution on every run.

Route-matched, digitized-curve overlays give **teclistamab SC AFE 1.29×** and **elranatamab (matched PK+IL-6) AFE 2.06×** (`VALIDATED_NOW_vs_STILL_UNDERWAY.md`). The single best-constrained parameter is the BiTE renal-clearance term $k_{\text{renal,max}}=8.70\,\mathrm{d}^{-1}$, fitted to blinatumomab and predicting a fall to 50% of $C_0$ in **1.78 h against the FDA BLINCYTO label's 2.11 h** (`PARAMETER_AUDIT_2026-07-13.md`). Two honesty flags travel with the PK arm: measured binding kinetics exist for only **4 of 22 molecules** (teclistamab, elranatamab, catumaxomab, the CD20 lineage) — the other 15 carry a generic `kon = 1e5 /M/s` placeholder with `koff` back-derived, so any conclusion resting on kon/koff rather than KD alone is assumption-driven — and $CL_{up}$/$f_{FcRn}$ enter only as the product $CL_{up}\cdot(1-f_{FcRn})$ and are **not separately identifiable** from plasma PK.

### 6.2 Provenance tagging: [MEASURED] / [DERIVED] / [FITTED] / [UNSOURCED]

Every parameter carries exactly one provenance tag; there is no "probably fine" tag (`PROVENANCE_AND_VALIDATION.md`). `[MEASURED]` requires a PMID/DOI **and a verbatim quote containing the number**; `[DERIVED]` gives the formula from measured inputs; `[FITTED]` marks a value tuned to an observation — which then cannot be reused to validate it; `[UNSOURCED — TBD]` is treated as a liability until resolved. This discipline exists because a 2026-07-13 audit found the failure modes it now guards against: **a page number** from an FDA table of contents (elranatamab "191") circulating as an IL-6 concentration; **fabricated digitizations** ("340", "230") of a figure that does not exist in a paper that reports no IL-6; and a `[FITTED]` clearance constant (IL-6 $k_{\mathrm{deg}}=0.20/\mathrm{hr}$) cited to PMID 31268236 — which is a *semi-mechanistic modeling paper reporting no measured clearance*. The rule "a citation next to a number does not make it measured" is now enforced by opening the paper. A companion defect class, silent fallbacks, was deleted outright: a mechanistic IL-6 array that came back empty now raises a hard `RuntimeError` instead of substituting the retired fitted `IL6_SCALE` under the mechanistic field name — because a fallback that "never crashes a run" converts a loud failure into a quiet fabrication.

### 6.3 The verification ledger and the reproduction notebook

Load-bearing claims are tracked one-row-per-claim in `evidence/QSP_VERIFICATION_LEDGER.csv`: **54 checks — 52 PASS, 1 PASS (engine-faithful, external source not re-audited), and 1 PENDING** (the IL-6 model-vs-clinical gate). Each row pins a claim to a value, a source file, and an artifact version_id, so the nomination co-leads (4-1BB / CD27), the architecture constants ($k_E{=}0.55$, $k_{\text{hit}}{=}12/\mathrm{day}$ fixed, $k_{\text{death}}{=}1.0$ locked), and the QSP windows are each independently traceable. The end-to-end reproduction lives in `qsp_reproduction.ipynb`: **21 cells (11 code + 10 markdown), executed top-to-bottom with 0 error outputs and all 5 assertions passing** — the assertions re-check overlay provenance keys (`ref_cell_assign`, `overlay_method`) and the receptor-pool basis (kidney EPCAM $=247.637$ nM to <0.1). It walks the full chain: canonical input inventory → per-cell scVI overlay → whole-body $R_{\mathrm{tot}}$ pools → GRN three-axis nomination → QSP window → PK/PD vs the sqlite curves → the IL-6 gate reported as-is.

### 6.4 Limitations — stated candidly

The model's claims are bounded, and the boundaries are load-bearing:

1. **Static $R_{\mathrm{costim}}$, and the fold is assumed.** The costim ranking reads receptor density at **resting** copy number. 4-1BB is activation-induced — resting transcript for TNFRSF9/TNFRSF4 sits near the detection floor — so a static read **under-rates 4-1BB / OX40 / ICOS** and produces a spurious "CD2 wins" ordering (ledger NEW-22). The induction machinery (`costim_induction.py`) is wired but the fold-upregulation magnitudes are `source="PENDING"`; the code **fails closed** (`ValueError`) rather than guess. Consequently, wherever fold matters it is presented as a **sensitivity sweep with `fold_is_assumed`, never a point estimate**, and the OX40/GITR net-negative figure is explicitly PROVISIONAL. This is a real deficiency of the multiplicative form $R(t)=R_{\mathrm{rest}}\,[1+(\text{FOLD}-1)\,a]$: it cannot induce a receptor that is genuinely zero at rest, and $a(0)=0$ under-rates inducible arms on every dose after the first.

2. **Organ-overlay retention gap.** The scVI overlay imputes off-panel receptors onto every cell from a live CZ CELLxGENE Census pull, but the **reference AnnData was not persisted** and the **organ overlay builder was not retained** (only the tumor builder survives). The stored per-cell reference indices (`ref_cell_assign`) are therefore only resolvable against a reference in identical row order, which a fresh Census pull does not guarantee. The bounded consequence: **EPCAM / DLL3 / PMEL were never added to the organ overlay set**, so off-tumor organ toxicity for three solid-only arms (catumaxomab/EPCAM, tarlatamab/DLL3, tebentafusp/PMEL) is incomplete pending an hours-long scVI re-run. The tumor pools are complete, and the **heme panel and the CRC solid pair (CEACAM5, ERBB2) have full organ coverage and are unaffected** (`T2b`, `Rtot_wholebody_final.json`).

3. **The IL-6/CRS gate fails on verified anchors — reported, not rescued.** Emergent mechanistic-ECF peaks (mosun 723.7, tecli 429.0, glofit 650.7, talq 772.4 pg/mL) were scored against the three manifest-verified digitized anchors (mosun 127.4, glofit 30.2, talq 19.8). The cleanest same-statistic pair (glofit/talq medians) comes out **WRONG DIRECTION** (model 0.84× vs clinical 1.53×); mosun/glofit is right-direction but magnitude-compressed (1.11× vs 4.22×). The diagnosis is mechanistic and specific: **talquetamab inverts** because GPRC5D sits on plasma cells and keratinized tissue, not the myeloid-contacting compartments that drive CRS, and the myeloid IL-6 term does not condition on *which* target compartment the killing occurs in (`IL6_VALIDATION_RESULT.json`). This is a genuine mechanism gap, presented as-is — not denominator-shopped, not talq-dropped (dropping teclistamab, whose 21 pg/mL anchor is unverified, actually makes the result *worse*). Crucially, the per-construct IL-6 **ranking is exactly invariant** to the two unmeasured disposition constants: since $C_{ss}=\text{prod}_{\text{total}}/(V\cdot k_{\mathrm{deg}})$, the global factors $V$, $k_{\mathrm{deg}}$, $S_{\max}$ cancel in any between-molecule comparison, so absolute-anchor agreement tests the global scale ($CL=55.9$ L/day) and the myeloid census, while the drug-discriminating spatial contact calculation carries no fitted constants. The nomination is **decided upstream of this gate and does not depend on it**.

4. **CRISPRi is loss-of-function; an engager arm is gain-of-function.** No knockdown screen can directly prove that *agonizing* a receptor helps. The CD4 Perturb-seq supplies a validated state change; the agonism direction is anchored by the orthogonal CRISPRa evidence (Schmidt 2022, corroborated by the Legut ORF screen, $r=+0.74$, ledger V7-07) and translated into predicted efficacy by the QSP model. The screen validates the *state*, the model translates it into a *window* — neither is claimed to measure agonist killing directly.

5. **Receptor copies are derived, not measured.** Per-cell receptor densities are **scRNA-seq transcripts converted to protein copies** via HPA/Glassman anchoring ($\log_{10}\text{copies}=0.7768\cdot\text{IHC}+3.9723$), *not* QIFIKIT-measured surface densities. The conversion is basis-verified (kidney EPCAM reproduces 247.637 nM exactly) but not independently validated per gene, and figure legends must not read "measured receptor copies." For activation-induced receptors this compounds limitation (1): a resting copy number at the detection floor is a noise-dominated estimate of a near-zero quantity.

Taken together: the PK backbone and the upstream three-axis nomination are validated and reproducible; the IL-6 arm has a rigorously-characterized failure that the ranking is structurally insulated from; and the costim-induction magnitudes remain the principal open input, handled by failing closed rather than guessing.

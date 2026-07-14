# Claude Design Brief — Costim Counter-Screen Deck
*Single source of truth for building the presentation. Every number here is verified against the current v7 artifacts (6-track audit, 0 refuted findings). Where a number is provisional or a common mistake exists, it is flagged. Build from THIS file, not from DECK_NARRATIVE.md (that one is v6-era / stale).*

> **MAJOR UPDATE (latest cycle): the de novo AF3 binder campaign is COMPLETE.** What §5–§8 previously flagged as "deferred to cloud / future work" has now been executed: a **3,964-design AlphaFold3 screen → 152 refold → 39 ten-seed finalists**, across all four arms (4-1BB, CD27, CD3 redirector, CEA5 TAA). Both costim co-leads now have confirmed de novo binders with sub-Ångström fold self-consistency. This upgrades stage 7 (protein design) from "aimed" to "designed + interface-confirmed in silico." New numbers in §3.5, new assets in §4.5, full detail in §5. The older "not yet run" caveats are corrected in-place below.

---

## 0. THE FOUR INTAKE ANSWERS (what Claude Design asks first)

**1. What we built — the biological question + finding.**
Question: CD3 T-cell engagers activate CD4⁺ cells too, and adding a costimulatory arm (CD28, 4-1BB…) amplifies whatever those CD4 cells do — driving cytokine-release syndrome and Treg expansion. *Which costim receptor arm boosts CD8 effector killing WITHOUT feeding the CD4 toxicity programs?*
Finding: A genome-scale CD4⁺ Perturb-seq **counter-screen** scored every costim arm on one benefit axis (CD8 effector) and five CD4 liability axes (CRS, Treg-suppression, help-erosion, exhaustion, proliferation). A liability veto collapses 11 arms to exactly two clean co-leads: **4-1BB (TNFRSF9) and CD27**. CD28 — the strongest effector arm — is gated out on CRS + suppression + proliferation. The nomination feeds a whole-body PBPK-QSP model that translates the three axes into a predicted therapeutic window, and de novo antibody binders are designed against the winners.

**2. Track:** Researcher. Judges have domain expertise. Rubric: Impact 25 / Claude-Use 25 / Depth 20 / Demo 30.

**3. Form:** Animated, storyline-driven slide deck (screen-recorded), **headline-first with a progressive reveal** — open with a big attention-grabbing headline (the "complete discovery pipeline in one week" hook, §1.5), then progressively reveal what was done via a 7-stage pipeline bar that lights up stage-by-stage as the science unfolds. Embedded RFdiffusion binder-design **timelapse animation** as the visual centerpiece. Design-consistent high-quality figures throughout.

**4. Assets available (all as artifacts — IDs in §4/§5):** ~15 publication figures already rendered; the audit digest (verified numbers + sources); the binder diffusion animation package (hero stills + movie); **the completed AF3 binder campaign — 39 ten-seed finalists across 4 targets, graded testable panels, funnel + CDR-property analysis figures**; the QSP model architecture map; the full verified data tables. The binder campaign that was "in progress" is now finished end-to-end.

---

## 1. THE ONE-SENTENCE STORY
**"The wrong-subset dataset is the right instrument: a CD4 screen can't measure CD8 killing, but it's the one place the toxic CD4 sub-program can be separated from effector help — so we use it as a counter-screen to nominate a costim arm (4-1BB / CD27) that boosts the killers while starving the storm."**

## 1.5 THE HOOK & THE THROUGHLINE (the Claude-as-multiplier story — this frames the whole deck)

**The headline (open with this — it's the attention-grab before any biology):**
This is a **complete antibody-therapeutic discovery campaign** — the same end-to-end process a drug-discovery lab runs to nominate and de-risk a candidate — **compressed into one week, run by one scientist directing Claude across every discipline it normally takes a multidisciplinary team to cover.**

The pipeline that was executed, end to end:
**Literature review & problem framing → genome-scale target screening → data digitization & integration → genome-scale GRN inference → mechanistic whole-body PBPK-QSP modeling → 3-layer validation → de novo protein/binder design.**

Each of those is normally a separate specialist (computational biology, systems pharmacology, structural/protein engineering, clinical/regulatory lit). Here one expert spanned all of them at team-scale throughput because Claude carried the execution — writing and running the analysis code, pulling and aligning the deposited matrices, building the 1.65M-edge network and the QSP model, running the RFdiffusion campaigns, and continuously auditing its own numbers.

**Headline options (pick one, Claude Design):**
- *"A drug-discovery pipeline that normally takes a team a quarter — target to designed binder — run in one week."*
- *"One scientist. One week. Target nomination → mechanistic QSP → de novo binder. This is what an expert + Claude does to the timeline."*
- *"The whole antibody-development stack, compressed."*

**The progressive-reveal device (the visual throughline — use it to carry the whole talk):**
Put a **7-stage pipeline bar** on screen from slide 1. Start it dark/greyed. As each section of the talk is reached, that stage **lights up** (use the accent gold `#e8a33d`). By the closing slide the entire pipeline glows — and *that* is the payoff image: "every one of these was done, in one week, expert-directed + Claude-executed." The science slides slot *underneath* their pipeline stage, so the audience always sees where they are in the campaign.

Stage → talk-section mapping for the reveal:
1. **Literature review & framing** → the problem slide (costim toxicity, name the sub-program)
2. **Target screening** → the counter-screen, 3 axes, the veto (11→2)
3. **Data digitization & integration** → the "unbiased confirmation" / genome-wide scan (pulling + aligning the full deposited matrices, Schmidt CD8, TAA cohorts, CITE-seq)
4. **GRN inference** → the mechanism/topology rings (1.65M-edge signed network)
5. **Mechanistic QSP** → the therapeutic-window translation
6. **Validation** → the rigor slide (3 layers, 0 refuted, self-retraction)
7. **Protein design** → the molecule / RFdiffusion animation

**Honest framing (keep it credible to expert judges — see guardrail #9):** the claim is **expert-directed augmentation**, not autonomous AI. A domain expert made every scientific decision, corrected course repeatedly, and supplied the hypothesis and house methods; Claude was the force multiplier on *execution breadth and speed*. That framing is both true and far more persuasive to a skeptical expert panel than "AI did drug discovery."

---

## 2. NARRATIVE ARC (headline-first, progressive-reveal)
0. **THE HOOK** — the headline (§1.5) + the 7-stage pipeline bar previewed (all dark). One line: *"Here's a complete discovery campaign — watch it get built."* Then the one-sentence science story (§1). This is the attention draw before any detail.
1. **The problem** *(lights stage 1)* — CD3 engagers hit CD4 too; costim amplifies CRS + Treg. Name the enemy: not the CD4 lineage, but a CD4 *sub-program*.
2. **The instrument** *(stage 2)* — the Marson/Pritchard genome-scale CD4 Perturb-seq. Why a "wrong subset" screen is the right tool (it uniquely resolves the suppressive + cytokine programs).
3. **The 3 axes** — benefit (CD8 effector, Schmidt CRISPRa) vs liability (CRS, Treg-suppression) → the efficacy/toxicity landscape.
4. **The veto** — 6-axis gate collapses 11 → 2. CD28 falls despite top effector. *Co-leads: 4-1BB + CD27.*
5. **Unbiased confirmation** *(stage 3 — data digitization/integration)* — genome-wide effector-first scan re-derives the same two; novel hits deselect on mechanism.
6. **Mechanism/topology** *(stage 4 — GRN)* — GRN rings: 4-1BB is a CD8-selective effector hub (~19× CD8 vs CD4).
7. **QSP translation** *(stage 5)* — three axes → one predicted therapeutic window; Treg-aware arm beats pan-costim.
8. **The molecule** *(stage 7 — protein design)* — TAA slate (CEACAM5/6 lead) + the **completed** de novo binder campaign: 3,964 AF3-screened designs → 39 confirmed finalists, co-lead binders for both 4-1BB and CD27 (§3.5). The RFdiffusion agonism-epitope animation lands here; the funnel/panel figures (§4.5) show the campaign was real and unbiased.
9. **Rigor** *(stage 6 — validation)* — 3 verification layers, 0 refuted, and we retracted our own best number (see §6).
10. **THE PAYOFF** — the full pipeline bar now fully lit. Restate the multiplier: *"target → designed binder, one week, one expert + Claude."* This is the Claude-Use (25%) close.

*(Note: stage 6/validation is threaded throughout — the continuous auditor caught and forced fixes on live numbers during the build — so you can also light stage 6 incrementally rather than only at slide 9. That "it audited itself as it went" detail is a strong Claude-Use point.)*

---

## 3. VERIFIED HEADLINE NUMBERS (quotable — cite ONE source per slide)

**Nomination (source: COSTIM_FINAL_3AXIS_SCORE_v7.csv):**
- CLEAN set = **{4-1BB, CD27} only**; 9 of 11 arms GATED.
- Gate tags: CD28 GATED[CRS,SUPP,PROLIF]; ICOS/CD30/CD40 GATED[HELP,PROLIF]; DNAM1/OX40/HVEM/DR3 GATED[SUPP,EXH]; GITR GATED[SUPP].
- Effector benefit (E_schmidt_z): **CD28 = 12.11 (panel max, but gated)**; CD27 = 4.28; 4-1BB = 3.74.
- Monte-Carlo nomination: joint P(both co-leads top-2) = **0.84** (expression-weighted); CD27 alone P(top-1) = **0.74**.
- Expression selectivity (CD8−CD4conv): **CD27 = +0.46 (CD8-selective)**; 4-1BB = +0.04 (neutral, rests on network argument); CD28 = −0.09 (CD4-biased).

**Discovery / GRN:**
- Genome-scale signed GRN backbone = **1,653,594 edges** (coef −22.9…+44.7; 41% negative).
- 4-1BB network effector-hub selectivity = **+54.4** (CD8 57.4 vs CD4 3.0, ≈19×).

**QSP (source: A31b_QSP_rerun_GRN_vs_DE.csv):**
- 4-1BB window = **+1.57** under BOTH drive models, uncapped (robust co-lead); TI up to 62 (DE) / 744 (GRN).
- CD27 window is **drive-dependent**: −2.37 (DE) → +1.27 (GRN), CRS-capped both. ← footnote this if presenting as equal co-lead.
- CD28 window = **−2.66 (DE) / −3.03 (GRN)** — best effector, erased by liver/CRS cap.
- Deck worked-example scale (qsp_worked_example_3way.png): CD28→Treg-aware swap buys **+5 window units**.

**TAA (source: TAA_finalists_6.csv):**
- 6 surface-valid CRC finalists: CEACAM6, CEACAM5, ITGB4, TSPAN6, LY6E, DPEP1 (say "6 selected finalists," NOT "top 6").
- Leads: **CEACAM5 (z=2.39) & CEACAM6 (z=2.37)**, both clinically precedented. Binder target = CEACAM5 single **B3 Ig-domain (res 593–675)**.

### 3.5 BINDER CAMPAIGN — de novo AF3
_(source: T1_funnel_attrition.csv, panel39_final.csv, finalist_combined_verdict.csv v2)_
- Funnel: **3,964 designs** AF3-screened (1-seed) → **152** cleared the refold gate → **39** finalists (10-seed), across 4 targets. (Of the 152, 150 returned 3-seed scores; 2 4-1BB folds dropped on a compute node — the same straggler batch that left the screen at 99.75%. Use **152** as the funnel number; it matches T1 and the funnel figure.)
- Scoring = 3 orthogonal gates: **ipSAE** (interface, MIN-of-directions), **fold-scRMSD** (self-consistency, binder-on-binder Kabsch, ChimeraX-validated), **10-seed consistency**.
- Testable panels (LEAD+CONFIRM, fold<2Å): **4-1BB 2 · CD27 5 · CD3 8 · CEA5 12 = 27 designs**.
- **4-1BB co-lead binder `41bb_r0867`:** ipSAE **0.49**, fold-scRMSD **0.86 Å**, pDockQ 0.43. VHH.
- **CD27 co-lead binder `cd27_r0011`:** ipSAE **0.55**, fold-scRMSD **0.60 Å**. VHH. (backup `cd27_r0846` ipSAE 0.45, fold 0.44 Å, seed-stable — safest single pick).
- CD3 redirector lead `cd3_r0478` ipSAE 0.49; CEA5 TAA lead `cea5_r0075` ipSAE 0.59 (12/12 LEAD-class — deepest panel).
- **Selection is unbiased:** retained-vs-filtered designs are separated by AF3 interface metrics by ~90 orders of magnitude more than by sequence (bt_pae_min p=1.6e⁻⁹⁴, interface-iPTM p=9.8e⁻⁹⁷). CDR sequence properties are essentially non-discriminating — GRAVY/aromatic/length/liabilities all p>0.05; the one exception, cdr_frac_charged, reaches nominal significance (p=0.013) but with a trivial effect size (Δmedian −0.014, 0.136 vs 0.150) — not a compositional bias. MPNN score itself does not predict AF3 success (p=0.78).
- **Developability clean:** 0 unpaired cysteines across all 3,964 designs; **all 35 fold-passing finalists have net-negative paratopes** (−0.8 to −4.0, complementary to acidic TNFR/CEACAM surfaces) — the lone net-positive design (41bb_r0022, +1.0) is a dropped fold-failure, not a testable binder.
- Screen completeness: 3,964 / 3,974 = **99.75%** (10 stragglers dropped on one node, never scored, none are CD27/CEA5).

**Verification:**
- 30/30 engineering ledger PASS · 25-hypothesis register **0 refuted** · 10/10 primary-source re-derivation PASS.
- Binder-metric integrity: fold-scRMSD formula bug (double-sqrt) caught by continuous audit and corrected; verdict unchanged (co-leads pass either way) — another self-audit credibility point.

---

## 4. FIGURE / ASSET MANIFEST (artifact IDs → what they show → slide)
*All are current latest_version_ids. Load with the {{artifact:ID}} embed.*

| Figure | Artifact / version | Shows | Slide |
|---|---|---|---|
| Efficacy-toxicity landscape | `dbff8d33-6cf9-41a0-90a3-89fb9e9ed129` | CD8 benefit vs worst CD4 liability; co-leads below veto, CD28 gated far-right | 3–5 |
| A17 genome-wide discovery | `788514eb-17c8-47f6-8cac-175a6b94e740` | Effector-first scan; known leads survive, novel hits deselect | 5 |
| GRN ring CD8 effector | `78cbc55c-f767-4c99-8811-e2d6de092e4d` | 4-1BB dark-red hub in CD8 | 6 (pair) |
| GRN ring CD4 effector | `b126804c-077d-4a00-acbd-62e2b2e4095f` | 4-1BB collapses in CD4 (colorbar ~7× smaller) — selectivity | 6 (pair) |
| QSP worked example (3-way) | `9e819406-84b4-4d25-8c98-c91fd088a615` | Window bars + CD28→Treg-aware swap = +5 units | 7/8 (in_deck) |
| TAA funnel | `04e5a5b8-23dc-49df-8940-084c45df0d0f` (v `e37773af…`) | 33,694 genes → 15-antigen slate; surface-QC caught ER decoys | 8 |
| **Binder hero still (41BB_23)** | `404998c8-afd0-4476-9eed-d670c486933d` | de novo VHH caps 4-1BB CRD1 agonist site, over the gold hotspots | 8/9 |
| **Epitope comparison** | `0f7b086b-5dc1-4eac-90c5-45ddad9303e5` | our binder = urelumab agonist site, NOT ligand/utomilumab | 8/9 |

**Master data (for any custom figure Claude Design wants to build):**
- audit_digest.json `aa453899-61fb-4e59-9b10-7597240ff5ca` — all numbers + sources + flags
- COSTIM_FINAL_3AXIS_SCORE_v7.csv `2490744f-4d2a-4af5-bfbd-3141f3d7bdba`
- A31b_QSP_rerun_GRN_vs_DE.csv `88c264d7-1ec3-431c-ae1a-710a7827411f`
- TAA_finalists_6.csv `470be6af-9d51-4be4-8be0-6a84cb923362`
- MODEL_ARCHITECTURE_MAP.md `18b7aa7a-462f-46f6-b10c-b0f5fc59b4bf` (current row; ignore older sibling)

### 4.5 BINDER-CAMPAIGN ASSETS (new — the completed AF3 screen)
| Figure / table | Artifact / version | Shows | Slide |
|---|---|---|---|
| **Funnel + interface separation** | `b0472532-f11c-4725-a646-b975051adf6d` | 3,964→39 attrition per target; retained-vs-filtered bimodal split at the refold gate (AF3 selects on interface, not sequence) | 8 |
| **CDR properties (4-panel)** | `38627667-8533-40b5-a44a-f2ba45b8e71f` | fold-passing paratopes net-negative (dropped fold-fails shown as open markers); finalists inside design envelope; CDR-H3 length; developability liability burden | 8/9 |
| Testable panels doc | `e5ecdad6-4984-49c9-8b30-4206e5b248c7` | graded per-target panels, selection logic, key reads | ref |
| State-recovery index | `ce0dcbc1-17fa-40a1-a2be-9f464ef9cb45` | master artifact index for the whole AF3 campaign | ref |

**Binder data tables (for custom figures):**
- panel39_final.csv `326c8d3f-37a1-486a-810c-b38b97c1a818` — 39 finalists, all metrics, panel_class
- af3_analysis_universe.parquet `a70c2829-1de8-493a-8007-a3c0b3679827` — 3,964 designs × 48 cols (every metric + CDR seq + properties + liabilities + funnel stage); slice any further comparison off this
- finalist_combined_verdict.csv (v2, corrected) `03829d25-502f-4a4e-8297-0b7d71c112c4`
- T1–T4 analysis tables: `7148f1dc…` (funnel) · `7bb770f4…` (retained-vs-filtered stats) · `a78b5548…` (finalist head-to-head) · `c435d071…` (per-target CDR profile)
- Co-lead structures (render in Mol*): 41bb_r0867 `27817fb8-b678-4811-b197-2332f798948e` · cd27_r0011 `51f0769e-4d52-4f23-ad63-5ae57889ec8d` · cd3_r0478 `78c46746-7621-4add-aab3-c2fb74cfb2d4` · cea5_r0075 `ad4c7f5d-9bc8-45c6-b68e-63d532da7898`

---

## 5. THE ANIMATION PACKAGE (visual centerpiece)
De novo RFdiffusion binder-design timelapse for the nominated 4-1BB arm.
- **What it shows:** a VHH binder diffuses from noise into a folded domain that docks the 4-1BB CRD1 agonist epitope, with the receptor anchored in a lipid bilayer — the actual therapeutic geometry.
- **The scientific payload (verified):** the de novo binder independently recapitulates the **urelumab-class agonist epitope** (ligand-independent, clustering-compatible: 8.9 Å from the 4-1BBL trimer, 0 clashes) and avoids the utomilumab blocking site and the ligand groove entirely. That's the agonism rationale, structurally.
- **Animation hero design = 41BB_23** (chosen from 100 RFdiffusion backbones by direct-overlay scoring; body sits over the hotspots, 3.1 Å lateral). This is the design shown recapitulating the urelumab epitope in the animation/stills.
- **Status:** hero stills + epitope figure DONE (§4); membrane-anchored movie + CD3-side-by-side DONE.
- **IMPORTANT — two different 4-1BB designs, don't conflate them (see guardrail #7):**
  - *Animation hero `41BB_23`* = the epitope-overlay showcase (RFdiffusion backbone, urelumab-site recapitulation). Chosen for the *agonism-geometry story*, not from the AF3 screen.
  - *AF3 screen co-lead `41bb_r0867` (backbone 62)* = the top interface-confirmed binder from the completed 3,964-design AF3 campaign (ipSAE 0.49, fold 0.86 Å). Chosen for *binding confidence*.
  - These are separate designs. The epitope-recapitulation claim is verified for `41BB_23`; the AF3 interface-confidence claim is verified for `r0867`. Keep the two stories on their own designs.
- **AF3 confidence is now DONE, not deferred:** the "binding-confidence (AF3 iPTM) is the deferred cloud step" caveat from earlier drafts is **superseded** — the full AF3 screen ran (3 cloud A100 pods) and produced the 39-finalist panel with ipSAE + fold self-consistency (§3.5). Present binders as *designed AND in-silico interface-confirmed*. The one remaining honest limit: no wet-lab binding yet (SPR/BLI KD is the bench next step).

---

## 6. TRUTHFULNESS GUARDRAILS (do-NOT-say list — from the audit)
These are the traps. Every one is a real inconsistency the audit caught.
1. **OX40/GITR "net-negative kill in Treg-rich tumor"** — NOT artifact-backed (in-chat self-report only; the saved model can't produce negative kill). **Do not present as an established result.** Provisional at most.
2. **Two QSP window scales exist** (A31b vs deck worked-example). Same ordinal story, different absolute numbers. **Quote ONE source per slide and label it.**
3. **CD27 is drive-dependent** (fails window on DE, passes on GRN). If shown as equal co-lead to 4-1BB, footnote it. 4-1BB is the robust one.
4. **H09 patient-stratification:** cite **CD137 ADT z = 11.2**, NEVER z = 37.28 (the latter was inflated by circular stratification and self-retracted — the retraction is itself a credibility asset, tell that story).
5. **CD30 is NOT a lead.** The GRN nomination panel-b ranks CD30 #1 — that's a PRE-VETO exploratory rank; CD30 is gated (HELP+PROLIF). Caption any CD30-topped figure as exploratory.
6. **TAA:** say "6 selected finalists" (not "top 6" — they're ranks 1,2,3,8,11,12). TSPAN6 is user-retained (single-cohort). LY6E is retained for clinical precedent, NOT tumor-restricted (negative z). CEACAM5 target = single-domain **B3**, not A3B3.
7. **Binders (UPDATED — campaign now complete):** ProteinMPNN + the full AF3 screen **ARE now done** (3,964 → 39 finalists, §3.5). Present binders as *designed and in-silico interface-confirmed*. Two remaining honest limits: (a) **no wet-lab binding data yet** — SPR/BLI KD is the explicit bench next step, don't imply measured affinity; (b) **the animation-hero design `41BB_23` and the AF3 co-lead `41bb_r0867` are different designs** — the urelumab-epitope recapitulation is verified for `41BB_23`, the AF3 interface confidence for `r0867`; don't merge the two claims onto one design. Also: 4-1BB's panel is genuinely thin (2 designs pass both gates — hard VHH-vs-cysteine-rich-TNFR target; only 23 of its 981 screened designs cleared the refold gate = 2.3%, the lowest per-target pass rate vs CD27 5.2% / CD3 3.6% / CEA5 4.2%), so present 4-1BB as "hard target, 2 confirmed binders," not a deep panel; CEA5 is the deep one (12/12).
8. **DECK_NARRATIVE.md is stale (v6).** Do not lift numbers from it.
9. **The Claude-multiplier claim = expert-directed augmentation, NOT autonomous AI.** Say "one expert directed Claude across every discipline" / "team-scale throughput," never "Claude did drug discovery by itself." A domain expert framed the hypothesis, supplied house methods, and made every scientific call. Overclaiming autonomy will read as naïve to expert judges and undercut the (true, stronger) efficiency-multiplier point. Also don't put a hard number on the speedup (e.g. "10×", "a quarter's work in a week") as a measured fact — it's an illustrative comparison; phrase as "normally a multidisciplinary team over months," not a benchmarked figure.

---

## 7. DESIGN SYSTEM
**Palette (established across existing figures — keep consistent):**
- CLEAN / nominated arm: `#c0392b` (deep red)
- GATED arm: `#95a5a6` (muted grey)
- CD3 / redirector: `#2c3e50` (dark slate)
- TAA / delivery: `#2e86c1` (blue)
- Accent / hotspot / benefit highlight: `#e8a33d` (gold)

**Structural-figure style (from the user's ChimeraX house style — match for any new 3D render):**
- Membrane: CG beads, **PO4 headgroups orange, NC3 goldenrod, glycerol grey, tails white**.
- Lighting: soft + shadows (intensity 0.7), `material dull`, silhouettes width 1.5, white background.
- Binder during diffusion: **ball form (atoms), never tube/ribbon**. Avoid the red/white/blue rainbow — use a cleaner single-family gradient.
- Camera: membrane horizontal, ectodomain up, binder epitope facing camera.

**Typography/plots:** publication-grade; data fills ≥40% of axes; no overlapping labels; one idea per panel.

---

## 8. ROADMAP (what's finished vs in progress)
**DONE & verified:** 3-axis nomination (v7) · genome-wide effector-first confirmation · GRN backbone + rings · QSP window translation · TAA slate · verification (3 layers, 0 refuted) · ~15 figures · binder hero stills + epitope figure · membrane-anchored binder animation (4-1BB + CD3) · **the full de novo binder campaign: RFdiffusion backbones + ProteinMPNN + AF3 3,964-design screen → 39 ten-seed finalists across all 4 targets, with graded testable panels + funnel/CDR-property analysis (§3.5, §4.5).**
**BENCH NEXT STEP (explicit future work):** wet-lab expression + binding validation (SPR/BLI KD, agonism cell assay) of the 27-design testable panel; the in-silico campaign hands off a ranked, developability-clean shortlist.
**OPTIONAL / DEFERRED:** structure-derived interface metrics (SASA, contact-residue epitope mapping) from the finalist .cif files; additional Xenium spatial tissues for the healthy-tissue safety map; 4-1BB CRD1 agonist-epitope compatibility check on the AF3 co-lead r0867 (only checked for animation-hero 41BB_23 so far).

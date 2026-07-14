## 5. Discovery modules: redirector co-expression & kinetic onset

The three-axis nomination (effector benefit vs. suppression and CRS liability)
is a static, magnitude-based scoring of each costim receptor. Two discovery
modules add orthogonal, single-cell and temporal context that the static score
cannot see: **A22** asks *where*, at single-cell resolution, a CD3-redirector
arm and each candidate costim receptor physically co-occur — and on which
compartment (CD8 killer vs. CD4-conv/Treg); **A24** asks *when* each receptor's
CRS and suppression liabilities appear — early (within the first dose interval)
vs. late (accruing over a multi-dose course). Both are hypothesis-generating
layers that sit beside, and support rather than establish, the nomination. The
CLEAN co-lead call — 4-1BB (TNFRSF9) and CD27 — is unchanged; these modules
refine *how each co-lead earns its place* and translate the liability profile
into dosing-relevant terms.

### 5.1 Redirector × costim co-expression (A22)

**Design.** A CD3 engager co-engages whatever T cell its anti-CD3 arm binds;
adding a costim arm then amplifies that cell. The delivery question is therefore
single-cell: on the same cells where the redirector (CD3) sits, is the costim
receptor also present, and in which compartment? A22 answers this on the RTCC
CITE-seq dataset (GEO GSE292621), a surface-proteome map of engager-activated
primary T cells — a different, orthogonal dataset from the hero CD4 CRISPRi
screen. Cells were gated on CLR-normalized surface markers into CD8 (n = 5,412),
CD4-conv (n = 3,156), Treg (n = 2,262) and other (n = 14,518), with a CD3⁺
threshold of CLR = 0.9214. Two readouts are reported per receptor: (i) a
**co-localization bias**, `copos_CD8_minus_tox` = (% of CD3⁺ CD8 cells that are
also receptor⁺) − mean(same fraction in CD4-conv, Treg), where positive means
the redirector and receptor preferentially co-occur on killers; and (ii) a
**graded** within-compartment Spearman correlation between CD3 and receptor
surface signal.

**Result.** CD27 is the panel standout on co-localization bias: **+17.6
percentage points** — 84.4% of CD3⁺ CD8 cells are CD27⁺, versus 60.5%
(CD4-conv) and 73.1% (Treg). It is the most CD8-biased receptor in the panel,
and the measurement is protein–protein (surface CD3 × surface CD27), the
higher-confidence modality. 4-1BB, by contrast, is essentially **compartment-flat
(−1.4 pp**; 57.5% CD8 vs. 57.1% CD4-conv / 60.8% Treg) — co-expressed with the
redirector across all three compartments, neither CD8-preferential nor
toxicity-biased. At the opposite pole sit the receptors whose co-localization is
skewed *toward* the toxicity compartments: DNAM1/CD226 is the most CD4/Treg-biased
at **−18.1 pp**, followed by TNFRSF25/DR3 (−9.3), OX40/TNFRSF4 (−8.4) and CD28
(−7.6). RNA-only arms (CD30, HVEM, GITR, DR3, DNAM1-panel RNA measurements) use a
noisier cross-modality (CD3-protein × receptor-RNA) estimate and are flagged
accordingly.

![A22 — Redirector (CD3) × costim receptor co-expression across RTCC-engaged CD8 / CD4-conv / Treg compartments. (a) Co-localization bias (CD8 minus mean CD4-conv/Treg, percentage points): CD27 +17.6 (most CD8-biased), DNAM1 −18.1 (most CD4/Treg-biased); 4-1BB near-neutral (−1.4). (b) Graded within-compartment Spearman co-expression with the redirector arm.]({{artifact:art_baa10775-cf7a-459a-b63b-43d57fa3f0c0}})

**Interpretation, with the honest nuance.** The co-localization axis is a
delivery argument, not a liability gate. For **CD27**, the redirector and the
costim target co-occur preferentially on CD8 killers, so a CD27 costim arm would
be delivered *cis* to the redirected effector — a single-cell rationale layered
on top of CD27's liability-cold profile. For **4-1BB**, co-expression is
compartment-flat, so co-localization neither helps nor hurts; 4-1BB's case rests
entirely on its liability-cold wiring (axes 2–3) and CD8 effector benefit
(Schmidt CRISPRa z = 3.74), not on preferential killer co-localization. The two
co-leads therefore clear the nomination by *different* logic, and the report does
not claim otherwise. One further caveat must be read directly off panel (b): the
**graded** CD3 × CD27 correlation is highest in Treg (Spearman 0.42), not CD8
(0.22). This is not a contradiction of the co-positivity result but a different
quantity — CD27 is so nearly ubiquitous on CD8 (84%) that little binary variance
remains there, whereas in Treg its *surface level* covaries more tightly with
CD3. The binary co-positivity (prevalence) favors CD8; the graded intensity
covariation favors Treg. CD27's co-localization advantage is a prevalence
statement, and should be cited as such.

### 5.2 Genome-wide kinetic onset of CRS and suppression (A24)

**Design.** A24 adds the time axis to the two liability programs. Using the hero
CD4 CRISPRi Perturb-seq log-fold-change layer (Zhu et al. 2025), an agonism proxy
= −1 × mean(log_fc over the program gene set) was computed at two stimulation
timepoints — Stim8hr (n = 11,415 regulators) and Stim48hr (n = 11,281) — and
z-scored across all regulators *within each timepoint separately*. The CRS
program is the storm-cytokine set TNF · IL2 · IFNG; the suppression program is
SUPP_full (15 genes). All **11,210 regulators present at both timepoints** are
scored; `onset_velocity` = z(48hr) − z(8hr), where positive means the program
builds late. Because the agonism proxy is derived by sign-flipping a knockdown
(loss-of-function) effect, it is a directional surrogate for agonism, consistent
with the screen's stated GOF/LOF limitation; cross-timepoint reading uses the
velocity term rather than raw z.

**Result — co-leads stay cold across time.** Both co-leads remain below the
fast-onset CRS gate (8hr z > 1.5) at **both** timepoints. 4-1BB is frankly
CRS-cold — CRS z = −0.76 (8hr) → −1.62 (48hr), below the genome mean at both times
(≈19th → 5th percentile) and declining further late. CD27 sits modestly above the
genome mean at 8 hr (z = +0.67, ≈77th percentile) but declines by 48 hr (z = +0.24)
and never approaches the fast-onset line. On suppression, both co-leads are below
the genome mean at both timepoints (4-1BB −1.10 → −0.80; CD27 −0.54 → −0.58) with
near-flat velocity — i.e., **no late-building suppression**.

**Result — contrast arms separate by onset shape.** CD2 (CRS z = +4.43 at 8 hr,
≈99.9th percentile) and CD28 (z = +2.32, ≈98.8th percentile) are the **only two
finalists above the fast-onset CRS gate**, and both spike early: CD2's CRS partly
decays by 48 hr (velocity −2.57) but stays high (+1.85); CD28 stays high (+1.49).
GITR is the sustained suppression hazard — the highest SUPP z at both timepoints
(+2.42 → +1.96, ≈98th/97th percentile). Two arms show a **late-onset** signature
an early PD sample would miss: OX40 builds on *both* axes (CRS −1.29 → +0.13,
velocity +1.41; SUPP +0.82 → +1.52, velocity +0.70), and CD40's CRS builds late
(−0.68 → +0.46, velocity +1.14).

![A24 — Genome-wide CRS/suppression onset kinetics across 11,210 regulators scored at Stim8hr and Stim48hr. Co-leads 4-1BB and CD27 stay below the fast-onset CRS gate (8hr z > 1.5) at both timepoints and below the genome mean on suppression; CD28 and CD2 are the only finalists that spike early (8 hr) on CRS; OX40 and CD40 build late.]({{artifact:art_e64c8a17-922a-4995-ac28-fb736bafad76}})

**Interpretation — onset maps to dose scheduling.** Onset timing is directly a
dosing variable. A liability already present at 8 h — within the first dose
interval — behaves like the TGN1412 precedent, in which an anti-CD28 superagonist produced
a systemic inflammatory response with rapid proinflammatory cytokine induction
within 90 minutes of a single intravenous dose in all six healthy volunteers
(Suntharalingam et al. 2006). This is the liability class that
constrains the priming dose and mandates conservative step-up / step-fractionated
dosing — the field-standard CRS-mitigation for CD3 bispecifics, where initial
step-fractionated dosing limits systemic T-cell activation and cytokine release
without compromising tumor response (Hosseini et al. 2020). By that
logic, CD28 and CD2 carry a first-dose penalty. A liability that instead builds
late (OX40 and CD40 on CRS; OX40, GITR, HVEM on suppression) would be under-read
by an early cytokine sample and accrue over a multi-dose course, eroding the
therapeutic window with cumulative exposure rather than at first dose. The two
co-leads have **neither** profile — no early CRS spike (no priming-dose penalty)
and no late-building suppression (window stable across repeat dosing). A24 is thus
the kinetic complement to the static liability-VETO gate: the co-leads are
liability-cold not only in magnitude but across the dosing-relevant time axis.

### 5.3 Scope and caveats

These two modules are exploratory layers. A22 is cross-sectional co-expression in
an engager-activated context on a separate CITE-seq dataset (GSE292621);
co-expression is association, not causation, and RNA-only arms are the noisier
cross-modality estimate. A24's agonism proxy is a sign-flipped CRISPRi
(loss-of-function) surrogate, z-standardized within each timepoint, so onset is
read through the velocity term. Neither module alters the nomination gate, which
remains the static three-axis liability-VETO score; they add a single-cell
delivery rationale (favoring CD27) and a dosing-kinetic rationale (favoring both
co-leads, penalizing the early spikers CD28/CD2 and the late-builders OX40/CD40).

---
*Data sources: A22_redirector_coexpression.csv/.png (RTCC CITE-seq GSE292621,
25,348 gated cells); A24_kinetic_onset_genomewide.csv/.png (hero CD4 CRISPRi
Perturb-seq, Zhu et al. 2025, 11,210 regulators at Stim8hr/Stim48hr). Key
references: Zhu et al. 2025 CD4⁺ Perturb-seq (bioRxiv, doi:10.64898/2025.12.23.696273);
Schmidt et al. 2022 Science (doi:10.1126/science.abj4008); Suntharalingam et al.
2006 N Engl J Med 355:1018–28 (doi:10.1056/NEJMoa063842, PMID 16908486);
Hosseini et al. 2020 npj Syst Biol Appl 6:28
(doi:10.1038/s41540-020-00145-7).*

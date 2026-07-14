# AF3 Screening Protocol for Costim-Agonist VHH Campaign
## (basis: Baker-lab de novo antibody lineage + our clustering-agonist filter)

## Provenance of the filter criteria — why AF3-iPTM, why these numbers

1. **Bennett et al. preprint (bioRxiv 585103v2, Mar 2024)** — the RFantibody paper.
   AF3 was NOT available at project outset; designs were filtered on a fine-tuned
   **RoseTTAFold2** (RF2) self-consistency check (design vs re-predicted complex).
   AF3 appears only retrospectively.

2. **Bennett et al., Nature (published Nov 5 2025)** — adds the AF3 result:
   AF3 **retrospectively discriminates VHH binders from non-binders, AUC = 0.86**,
   using **interface iPTM** as the metric; authors conclude AF3 outperforms
   RFdiffusion as a predictor of binding and should be an additional filter.
   => For antibodies/VHHs specifically, the validated AF3 discriminator is **iPTM**.

3. **Poddiakov et al., Sci Rep 2025 (10.1038/s41598-025-10241-5)** — the on-target
   4-1BB agonist VHH paper. Step-1 selection = **PAE_interaction < 10 AND Rosetta
   interaction score < -20 REU**; reference 7D4B VHH = -46.888 REU / 1313 Anstrom^2 interface.
   Confirms hotspot-conditioned VHH-on-4-1BB is a real, published route.

## The screen (cloud AF3 = folding; local = RFdiffusion+ProteinMPNN)

STAGE A  (LOCAL, A2000)  RFdiffusion_Ab: 100 backbones/target, VHH framework, CRD1 (4-1BB) / CD70-site (CD27) hotspots
STAGE B  (LOCAL, A2000)  ProteinMPNN: 8 seqs/backbone, temperature 1e-6, omit_AAs=CX, loops H1,H2,H3  -> ~800 seqs/target
STAGE C  (CLOUD, AF3)    Re-fold each binder+receptor complex. GATE:
                          - iPTM (interface)         primary rank + gate  [Bennett/Nature, AUC 0.86; >0.85 stringent / >0.6 permissive]
                          - RMSD(AF3 vs design) low   self-consistency      [Bennett/Nature: AF3 flags mis-designed VHH]
                            (threshold set from the design-set RMSD distribution, not an external cutoff)
STAGE D  (LOCAL)         CLUSTERING-COMPATIBILITY (agonist filter, ours, validated):
                          superimpose each AF3-passing binder onto all 3 receptors of the
                          6MGP 4-1BBL:4-1BB 3:3 assembly; require ligand-clash==0 AND inter-binder-clash==0.
                          VALIDATED on controls: urelumab (potent) PASS 0 clashes; utomilumab (weak) FAIL 13191.
                          => AF3 says "binds, confident, designed pose"; 6MGP says "will agonize".
STAGE E  (LOCAL->QSP)    finalists ranked by iPTM x clustering-margin -> feed KD_costim / n_costim knobs.

## CD27 vs 4-1BB asymmetry (mechanism-level)
- 4-1BB: target CRD1 (urelumab-class, non-ligand-blocking) -> clustering preserved + augmented.
- CD27 : target CD70 site (varlilumab-class, ligand-MIMETIC IS the agonist mode).
  => Stage D clustering filter is 4-1BB-specific; CD27 uses ligand-mimicry + multivalency.


## EXACT thresholds — Bennett et al., Nature 649:183 (1 Jan 2026), verbatim

- AF3 **iPTM** (interface predicted TM) "is predictive of binding success (area
  under the curve = 0.86)" for designed VHHs. iPTM is the PRIMARY rank metric.
- VHH enrichment anchor: "only 9% of our ordered VHH designs have an ipTM > 0.6,
  suggesting that success rate will be improved by incorporation of an ipTM filter."
  => **iPTM > 0.6 = permissive VHH gate** (top ~9%).
- scFv enrichment anchor: "Only 4% of the initial design library has ipTM > 0.85,
  whereas 5 out of 6 experimentally confirmed designs pass this threshold."
  => **iPTM > 0.85 = stringent gate** (finalist tier).
- AF3 setup for the fold: "MSA and templates for the target and only a template
  for the VHH (as CDRs are de novo, MSA would be of limited utility)."
- The RF2 fine-tuned filter they used at design time had "limited filtering power" —
  AF3 iPTM is the upgrade; use AF3, not RF2, for the confidence screen.
- Self-consistency rescue case: AF3 correctly predicted the ONE inaccurately-designed
  SARS-CoV-2 VHH's true (non-designed) structure; "had AlphaFold3 been used as an
  initial filter, this design would have been rejected" -> RMSD(AF3 vs design) filter.

## FINAL tiered gate for our campaign (every threshold traces to a read source)
  TIER-1 (order-worthy):   AF3 iPTM > 0.85           [Bennett/Nature: 5/6 confirmed designs pass; top ~4%]
                           AND clustering-compatible (4-1BB)   [ours, urelumab-validated]
  TIER-2 (expand pool):    AF3 iPTM > 0.6            [Bennett/Nature: top ~9% of ordered VHHs]
                           AND design self-consistency (low AF3-vs-design RMSD)  [Bennett/Nature]
  ROSETTA cross-check:     PAE_interaction < 10 AND Rosetta interaction score < -20 REU
                           [Poddiakov/Sci Rep 2025, on-target 4-1BB VHH Step-1 selection]
  AGONIST layer (4-1BB):   6MGP clustering-compat (ligand-clash==0, inter-binder-clash==0)  [ours, urelumab-validated]
  AF3 folding = CLOUD (or Boltz-2/Chai-1 iPTM analog, LOCAL); RFdiffusion+ProteinMPNN + clustering filter = LOCAL A2000.

## NOTE on removed content
  A prior draft cited a paper 'RFdiffusion3 (bioRxiv Jan 2026)' for numeric thresholds
  (PAE<1.5, pTM>0.8, RMSD<2.5A, 4-8 seqs/backbone). That reference was NOT retrieved or
  verified in-session and has been REMOVED. The 4-8 ProteinMPNN-seqs/backbone choice is a
  standard RFantibody-pipeline setting (we use 8, from the RFantibody nanobody example),
  not a claim from that paper. All remaining thresholds trace to Bennett/Nature or Poddiakov/Sci Rep.

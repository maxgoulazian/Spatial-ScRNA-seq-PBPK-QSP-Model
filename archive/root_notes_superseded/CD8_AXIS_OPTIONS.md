# CD8+ Effector Axis — Matched-Dataset Options
(counter-screen project: costim T-cell-engager arm nomination)

## The gap we are filling
The hero dataset (Zhu et al. 2025, CZI billion-cell-project, DOI 10.64898/2025.12.23.696273)
is a genome-scale Perturb-seq in primary human **CD4+** T cells — 12 files, 4 donors x
{Rest, Stim8Hr, Stim48Hr}. It supplies the CD4 SUPPRESSION and CRS liability axes.
It does NOT contain CD8 cells, so it cannot supply the CD8 EFFECTOR-BENEFIT axis on its own.

The CZI Virtual Cells Platform has only two namespaces: `cellxgene` (55,696 OBSERVATIONAL
atlases, incl. many "CD8 T" health-atlas sets — none perturbational) and
`billion-cell-project` (exactly our 12 CD4 files). => There is NO CD8 Perturb-seq companion
on the same platform. A matched CD8 dataset must come from GEO/literature.

## Why a CRISPRi/CRISPRa knockdown screen can't give agonism direction
An engager costim ARM is gain-of-function (agonism). CD4 CRISPRi is loss-of-function.
So the CD8 anchor should ideally be GAIN-of-function (CRISPRa / ORF overexpression) to
supply agonism-direction evidence. Two same-modality options below are gain-of-function.

================================================================================
## RANKED CD8 EFFECTOR-AXIS OPTIONS
================================================================================

### OPTION A (TOP PICK) — Legut et al. 2022, Nature "A genome-scale screen for synthetic
    drivers of T-cell proliferation"  (Sanjana lab, NYU)
  DOI: 10.1038/s41586-022-04494-7 | PMID 35296855 | GEO: GSE193736 (public, 2022-01-21)
  Modality: **GAIN-OF-FUNCTION** — ~12,000 barcoded human ORFs overexpressed.
  Readout: proliferation (CFSE-low) sorted **separately in CD8 AND CD4** primary human T cells,
    plus activation & cytokine secretion (IL-2, IFN-gamma).
  Single-cell arm: **OverCITE-seq** (5 modalities: sc_GEX + sc_ADT surface protein +
    sc_ORF perturbation call + sc_TCR + sc_HTO hashtag).
  Supplementary (directly usable, small):
    - GSE193736_ORF_screen_counts.csv.gz  (2.1 MB; cols: CD4_CFSE_low, CD4_presort,
      CD4_prestim, CD8_CFSE_low, CD8_presort, CD8_prestim, plasmid)  -> 11,610 distinct ORFs
    - GSE193736_bulkRNAseq_counts.csv.gz  (CD8 & CD4 bulk RNA-seq, LTBR vs tNGFR, rest/stim)
    - GSE193736_ATAC_counts.csv.gz ; GSE193736_RAW.tar (single-cell matrices)
  COSTIM PANEL COVERAGE in ORF library (16/22 present, multi-barcode):
    PRESENT: TNFRSF9(4-1BB), CD27, CD28, ICOS, TNFRSF18(GITR), CD40LG, TNFRSF14(HVEM),
             CD2, TNFRSF25(DR3), CD40, LTBR, PDCD1, HAVCR2, LAG3, TIGIT, TNFRSF1A, TNFRSF1B
    ABSENT : TNFRSF4(OX40), TNFRSF8(CD30), CD226(DNAM1), CTLA4, BTLA
  WHY IT FITS: directly measures whether AGONIZING (overexpressing) each costim receptor
    increases CD8 effector output — the exact direction an engager arm acts, and the exact
    evidence a CD4 CRISPRi screen structurally cannot provide. Same functional-genomics
    framework, CD8-specific sort. NOT the same lab as the hero (Sanjana, not Marson), but
    the strongest scientific match on modality + panel coverage.

### OPTION B (SAME GROUP) — Schmidt et al. 2022, Science "CRISPRa and CRISPRi screens decode
    stimulation responses in primary human T cells"  (Marson lab)  *** already downloaded ***
  DOI: 10.1126/science.abj4008 | GEO SuperSeries GSE174292
  CD8-specific readout = genome-wide pooled FACS screen **GSE174255**: IFN-gamma-sorted arm
    run in "Primary CD4- T cells" (= CD8), CRISPRa (Calabrese) + CRISPRi (Dolcetto) libraries.
    (matches proposal's "screening CD4+ on IL-2 and CD8+ on IFN-gamma").
  Single-cell CRISPRa Perturb-seq (GSE190604) is **bulk T cells, not CD8-sorted** — 74 targets,
    only 4/19 costim panel (TNFRSF9, CD27, CD28, CD2). So Schmidt's CD8 signal is POOLED-screen
    (guide enrichment), NOT single-cell. Good as a same-lab cross-check; weaker per-receptor
    resolution than Legut for the sc panel.

### OPTION C (ORTHOGONAL, CD8-DEDICATED) — McCutcheon et al. 2023, Nature Genetics
    "Transcriptional and epigenetic regulators of human CD8+ T cell function identified through
    orthogonal CRISPR screens"  (Gersbach lab, Duke)
  DOI: 10.1038/s41588-023-01554-0 | PMID 37945901
  Modality: pooled CRISPRa + CRISPRi (gain + loss) over 120 transcriptional/epigenetic
    regulators, in primary human **CD8+** T cells; phenotype = memory/effector state (BATF3 hit).
  Scope: focused (120 regulators, TFs/epigenetic) rather than genome-scale or costim-receptor-
    centric — complementary, not a costim-panel screen. Different lab.

### OPTION D (in-silico fallback) — subset CD8 in-silico from Schmidt sc arm (GSE190604)
  Features include CD8A/CD8B/CD4/IFNG/GZMB/PRF1 => CD8 subsetting IS technically possible, but
  the arm is bulk (mixed CD4/CD8) and only 4/19 costim panel present. Lowest-value path.

================================================================================
## RECOMMENDATION
================================================================================
Primary CD8 effector anchor  = OPTION A (Legut GSE193736, gain-of-function ORF, CD8 sort,
    16/22 costim panel, sc OverCITE-seq surface-protein layer bonus).
Same-lab cross-check          = OPTION B (Schmidt GSE174255 IFN-gamma / CD4- pooled screen)  [in hand].
Orthogonal regulator check    = OPTION C (McCutcheon, CD8 TF/epigenetic).
=> Legut best satisfies "perturb-seq or equally useful"; Schmidt best satisfies "same group".
   Using BOTH gives same-lab continuity AND gain-of-function per-receptor CD8 resolution.

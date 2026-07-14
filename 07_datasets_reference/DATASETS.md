# Costim T-cell-engager toxicity counter-screen — datasets

Project dir: /media/balthasar-lab/RAID4/costim_engager_counterscreen
Populated: 2026-07-07

## Three-axis nomination — which dataset feeds which axis

| Axis | Dataset | Accession / source | Status |
|------|---------|--------------------|--------|
| HERO — CD4 suppression + CRS counter-screen (IL-10/Treg + storm cytokines) | Marson/Pritchard CD4+ T-cell CRISPRi Perturb-seq (Zhu et al. 2025) | VCP "Primary Human CD4+ T Cell Perturb-seq"; bioRxiv 10.64898/2025.12.23.696273 | PENDING vcp login |
| CD8 EFFECTOR — gain-of-function IFN-g anchor | Schmidt et al. 2022 (Science 10.1126/science.abj4008) | GEO GSE174292 (SuperSeries) | DOWNLOADED |
| KILLING readout (supporting) | Shifrut et al. 2018 (Cell 10.1016/j.cell.2018.10.024), SLICE/CROP-seq | GEO GSE119450 | DOWNLOADED |

## Schmidt 2022 deposit (GSE174292 SuperSeries, BioProject PRJNA729110)

data/schmidt2022_tcell_perturbseq/
  GSE174255_genomewide_screen_readcounts/   <- CD8 EFFECTOR ANCHOR (the axis-3 file)
      GSE174255_sgRNA-Read-Counts.xlsx        genome-wide pooled CRISPRa (Calabrese) + CRISPRi (Dolcetto)
                                              IFN-g sorted in Primary CD4- (=CD8) T cells; IL-2 sorted in CD4+ T cells
                                              4 sheets: Calabrese SetA/B, Dolcetto SetA/B
  GSE190604_CRISPRa_Perturbseq/             <- CRISPRa Perturb-seq (scRNA); "Primary bulk T cells", NOT CD8-sorted
      GSE190604_matrix.mtx.gz (950M), barcodes.tsv.gz, features.tsv.gz,
      GSE190604_cellranger-guidecalls-aggregated-unfiltered.txt.gz
  GSE190846_supplemental_CRISPR_screens/
      GSE190846_supp_CD4_CRISPR_screens_read_counts.tsv.gz
  GSE174284_bulk_RNAseq/
      GSE174284_gene_counts_raw.txt.gz

NOTE ON "CD8 Perturb-seq": Schmidt's single-cell (Perturb-seq) arm was run on BULK T cells,
not CD8-sorted. The CD8-specific gain-of-function signal (the IFN-g CRISPRa hits: 4-1BB/CD27/
CD40/OX40) comes from the genome-wide POOLED FACS screen (GSE174255), where IFN-g was sorted
in CD4-depleted (CD8) T cells. That pooled CRISPRa screen is the correct axis-3 effector anchor.

## Shifrut 2018 deposit
data/shifrut2018_slice_killing/GSE119450_CROPseq/
    GSE119450_RAW.tar (328M)  -- CROP-seq in primary human T cells

## Provenance
All GEO files pulled from ftp.ncbi.nlm.nih.gov/geo/series/... on 2026-07-07.
MD5 checksums: manifests/geo_md5.txt


================================================================================
## CD8 EFFECTOR AXIS — RESOLVED (session 2)
================================================================================
Q: does the CD8 gain-of-function axis have gaps for the costim panel?
A: NO panel-membership gaps. 14 costim-arm candidates, all covered:
   - Legut GSE193736 ORF (single-cell OverCITE-seq + CD8 proliferation sort): 11/14
       gaps = CD226(DNAM1), TNFRSF4(OX40), TNFRSF8(CD30)
   - Schmidt GSE174255 genome-wide CRISPRa (Calabrese), IFN-gamma / CD4- (=CD8) sort: 14/14
       => covers the 3 Legut misses (OX40, CD30, DNAM1) at pooled-screen (guide-enrichment) level
   - COMBINED CD8 gain-of-function gaps: NONE.
   Checkpoints missing from Legut (CTLA4, BTLA) are irrelevant — inhibitory, you block not agonize.
Residual (nice-to-have, not blocker): OX40/CD30/DNAM1 have POOLED-screen effector evidence
   but not single-cell per-receptor transcriptomic resolution (absent from Legut sc arm).
Decision: CD8 effector axis is dataset-complete for nomination. Defer chasing single-cell
   resolution on those 3 unless one becomes a top nominee.
Coverage matrix: manifests/costim_cd8_axis_coverage.csv


================================================================================
## CD8 EFFECTOR-AXIS DATASETS — DOWNLOADED (session 2, "pull it all")
================================================================================
### Legut 2022 (Sanjana) GSE193736  -> data/legut2022_orf_overcite_gse193736/
  Genome-scale ORF overexpression (GAIN-OF-FUNCTION) screen, CD8+CD4 proliferation sort,
  + single-cell OverCITE-seq (5 modalities). Files (~18.6 MB total):
    GSE193736_ORF_screen_counts.csv.gz  (11,610 ORFs; cols CD8_CFSE_low/CD4_CFSE_low/presort/prestim/plasmid)
    GSE193736_bulkRNAseq_counts.csv.gz  (CD8 & CD4 bulk, LTBR vs tNGFR, rest/stim)
    GSE193736_ATAC_counts.csv.gz
    GSE193736_RAW.tar -> RAW/ : GEX(12M), ADT surface-protein(88K), ORF(40K), TCR(1.4M), HTO(55K)
  Costim panel in ORF library: 11/14 (misses OX40, CD30, DNAM1).

### McCutcheon 2023 (Gersbach) GSE218985-988 + GSE241933  -> data/mccutcheon2023_cd8_regulators_gse218991/
  Orthogonal CRISPRa+CRISPRi over 120 TF/epigenetic regulators in primary human CD8+ T cells.
  Functional screen tables (GSE241933, small): TFome-KO with/without BATF3 (6,999 gRNAs, Low/High x 2 donors),
    CRISPRa 1x/2xVP64 IL2RA, CRISPRa/CRISPRi TF CCR7, CRISPRi B2M, CRISPRi CD2.
  RNA-seq TPMs (GSE218986). RAW tars: GSE218985(2.1G RNA-seq), GSE218987(1.5G CUT&RUN), GSE218988(3.6G).

### Schmidt 2022 (Marson, SAME GROUP) GSE174255 [already in hand]
  Genome-wide CRISPRa(Calabrese)+CRISPRi(Dolcetto), IFN-gamma / CD4-(=CD8) sort. Costim panel 14/14.

=> CD8 gain-of-function coverage of the 14-member costim panel: COMPLETE (no gaps).
   Legut = per-receptor single-cell resolution (11/14); Schmidt CRISPRa = full panel at pooled level.
   OX40/CD30/DNAM1: pooled effector evidence only (no single-cell) — acceptable for QSP scalar score.

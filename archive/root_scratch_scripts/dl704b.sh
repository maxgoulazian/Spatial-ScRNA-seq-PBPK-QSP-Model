#!/bin/bash
cd /media/balthasar-lab/RAID4/costim_engager_counterscreen/data/exhausted_cd8_gse295704
for F in GSE295704_IPIHNSC126_T1_scrna_CD45_enriched_filtered_feature_bc_matrix.h5 GSE295704_IPIPOOL004_P1_scrna_CD45_enriched_filtered_feature_bc_matrix.h5 GSE295704_feature_reference.csv.gz; do
  aria2c -x8 -s8 -c --file-allocation=none "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295704/suppl/$F"
done
echo DL704_DONE

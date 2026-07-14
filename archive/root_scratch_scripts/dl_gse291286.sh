#!/bin/bash
D=/media/balthasar-lab/RAID4/costim_engager_counterscreen/data/activated_cd4_citeseq_gse291286
mkdir -p $D; cd $D
echo "START $(date -u +%H:%M:%S)"
curl -s --max-time 900 -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291286/suppl/GSE291286_RAW.tar"
curl -s --max-time 60 -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291286/suppl/GSE291286_feature_reference.csv.gz"
echo "TAR $(stat -c%s GSE291286_RAW.tar 2>/dev/null) bytes"
tar xf GSE291286_RAW.tar && echo "EXTRACTED $(ls | wc -l) files"
ls -la | head -20
echo "DL_DONE $(date -u +%H:%M:%S)"

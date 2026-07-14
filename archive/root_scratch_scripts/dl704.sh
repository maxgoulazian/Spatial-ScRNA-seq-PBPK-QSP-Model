#!/bin/bash
mkdir -p /media/balthasar-lab/RAID4/costim_engager_counterscreen/data/exhausted_cd8_gse295704 && cd /media/balthasar-lab/RAID4/costim_engager_counterscreen/data/exhausted_cd8_gse295704
aria2c -x8 -s8 -c --file-allocation=none -o RAW.tar 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295704/suppl/GSE295704_RAW.tar'
tar -xf RAW.tar
echo EXTRACTED_DONE

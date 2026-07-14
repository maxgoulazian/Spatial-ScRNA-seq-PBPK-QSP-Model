#!/bin/bash
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT=/media/balthasar-lab/RAID4/costim_engager_counterscreen/analysis/citeseq_rna2protein/cd30_sweep
mkdir -p $OUT; cd $OUT
srch () { curl -s --max-time 30 -G "$BASE/esearch.fcgi" --data-urlencode "term=$1" -d "db=gds&retmax=$2&retmode=json" | grep -oE '"20[0-9]{7}"' | tr -d '"'; }
> uids.txt
srch "CITE-seq" 500 >> uids.txt
srch "TotalSeq-C" 300 >> uids.txt
srch "antibody capture AND T cell" 300 >> uids.txt
srch "CITE-seq AND (lymphoma OR tumor OR activated OR stimulated)" 300 >> uids.txt
srch "AbSeq single cell" 150 >> uids.txt
srch "cellular indexing transcriptomes epitopes" 200 >> uids.txt
sort -u uids.txt | grep -E '^20[0-9]{7}$' > uids_u.txt
echo "UNIQUE_SERIES=$(wc -l < uids_u.txt)"
> gse.txt; while read u; do echo "GSE${u:3}"; done < uids_u.txt > gse.txt
echo -e "GSE\tpanel_file\thas_CD30\thas_CD137\thas_both" > sweep_results.tsv
N=0
while read GSE; do
  N=$((N+1))
  STUB=$(echo $GSE | sed -E 's/[0-9]{3}$/nnn/')
  LIST=$(curl -s --max-time 15 "https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/" 2>/dev/null | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"//')
  for AB in $(echo "$LIST" | grep -iE 'feature|antibody|adt|totalseq|abseq|panel|tag_?ref|hashtag' | grep -ivE 'matrix.mtx|barcodes.tsv.gz$' | head -3); do
    URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/$AB"
    CONTENT=$(curl -s --max-time 20 "$URL" 2>/dev/null | zcat 2>/dev/null | head -400)
    C30=$(echo "$CONTENT" | grep -icE 'CD30[^0-9A-Za-z]|TNFRSF8|Ber-?ACT35|BerH2')
    C137=$(echo "$CONTENT" | grep -icE 'CD137|TNFRSF9|4-1BB|4.1BB')
    if [ "$C30" -gt 0 ] || [ "$C137" -gt 0 ]; then
      BOTH="no"; [ "$C30" -gt 0 ] && [ "$C137" -gt 0 ] && BOTH="YES"
      echo -e "${GSE}\t${AB}\t${C30}\t${C137}\t${BOTH}" >> sweep_results.tsv
      break
    fi
  done
  [ $((N % 50)) -eq 0 ] && echo "...checked $N series"
done < gse.txt
echo "SWEEP_DONE checked=$N"
echo "=== HITS with CD30 ==="; awk -F'\t' '$3>0' sweep_results.tsv
echo "=== HITS with BOTH ==="; awk -F'\t' '$5=="YES"' sweep_results.tsv

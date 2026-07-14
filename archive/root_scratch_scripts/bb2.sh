#!/bin/bash
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
srch () { curl -s --max-time 30 -G "$BASE/esearch.fcgi" --data-urlencode "term=$1" -d "db=gds&retmax=40&retmode=json" | grep -oE '"20[0-9]{7}"' | tr -d '"'; }
> /tmp/bb2.txt
srch "T cell stimulation time course CITE-seq" >> /tmp/bb2.txt
srch "in vitro T cell activation multimodal single cell" >> /tmp/bb2.txt
srch "4-1BB CD137 T cell CITE-seq" >> /tmp/bb2.txt
srch "TotalSeq activation induced markers T cell" >> /tmp/bb2.txt
sort -u /tmp/bb2.txt | grep -E '^20[0-9]{7}$' > /tmp/bb2u.txt
echo "unique UIDs: $(wc -l < /tmp/bb2u.txt)"
> /tmp/bb2gse.txt; while read u; do echo "GSE${u:3}"; done < /tmp/bb2u.txt > /tmp/bb2gse.txt
echo -e "GSE\tCD137recAb\tpanel_file"
while read GSE; do
  STUB=$(echo $GSE | sed -E 's/[0-9]{3}$/nnn/')
  LIST=$(curl -s --max-time 15 "https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/" 2>/dev/null | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"//')
  for AB in $(echo "$LIST" | grep -iE 'feature|antibody|adt|totalseq|panel|tag_?ref' | grep -ivE 'matrix.mtx|barcodes.tsv.gz$' | head -2); do
    C=$(curl -s --max-time 18 "https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/$AB" 2>/dev/null | zcat 2>/dev/null | grep -ivE 'Gene Expression')
    REC=$(echo "$C" | grep -icE 'CD137|4-1BB|TNFRSF9')
    RECr=$(echo "$C" | grep -iE 'CD137|4-1BB|TNFRSF9' | grep -ivE 'CD137L|Ligand|TNFSF9|Gene' | head -1)
    if [ -n "$RECr" ]; then echo -e "${GSE}\tYES: ${RECr:0:40}\t${AB}"; break; fi
  done
done < /tmp/bb2gse.txt
echo DONE

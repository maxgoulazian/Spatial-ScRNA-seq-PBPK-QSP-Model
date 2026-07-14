#!/bin/bash
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
srch () { curl -s --max-time 30 -G "$BASE/esearch.fcgi" --data-urlencode "term=$1" -d "db=gds&retmax=40&retmode=json" | grep -oE '"20[0-9]{7}"' | tr -d '"'; }
> /tmp/bb_uids.txt
srch "CITE-seq AND (stimulated OR activated) AND T cell AND (CD137 OR 4-1BB)" >> /tmp/bb_uids.txt
srch "T cell activation CITE-seq TotalSeq-C CD137" >> /tmp/bb_uids.txt
srch "anti-CD3 CD28 stimulation CITE-seq single cell surface" >> /tmp/bb_uids.txt
srch "CAR T CITE-seq 4-1BB antibody" >> /tmp/bb_uids.txt
srch "tumor infiltrating T cell CITE-seq activation 4-1BB" >> /tmp/bb_uids.txt
sort -u /tmp/bb_uids.txt | grep -E '^20[0-9]{7}$' > /tmp/bb_u.txt
echo "unique: $(wc -l < /tmp/bb_u.txt)"
> /tmp/bb_gse.txt; while read u; do echo "GSE${u:3}"; done < /tmp/bb_u.txt > /tmp/bb_gse.txt
echo -e "GSE\tpanel\tCD137_recept_Ab\tCD30_Ab\tstim_context"
while read GSE; do
  STUB=$(echo $GSE | sed -E 's/[0-9]{3}$/nnn/')
  LIST=$(curl -s --max-time 15 "https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/" 2>/dev/null | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"//')
  for AB in $(echo "$LIST" | grep -iE 'feature|antibody|adt|totalseq|abseq|panel|tag_?ref' | grep -ivE 'matrix.mtx|barcodes.tsv.gz$' | head -2); do
    C=$(curl -s --max-time 18 "https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/$AB" 2>/dev/null | zcat 2>/dev/null | grep -ivE 'Gene Expression')
    REC=$(echo "$C" | grep -iE 'CD137|4-1BB|TNFRSF9' | grep -ivE 'CD137L|4-1BBL|Ligand|TNFSF9' | head -1 | wc -l)
    C30=$(echo "$C" | grep -iE '(^|[,\t| ])CD30([,\t| ]|$)|anti-human_CD30[,_ ]' | grep -ivE 'CD30L|CD300|CD303|CD153|CD304|CD305|CD301' | head -1 | wc -l)
    if [ "$REC" -gt 0 ]; then echo -e "${GSE}\t${AB}\t${REC}\t${C30}\t-"; break; fi
  done
done < /tmp/bb_gse.txt
echo DONE

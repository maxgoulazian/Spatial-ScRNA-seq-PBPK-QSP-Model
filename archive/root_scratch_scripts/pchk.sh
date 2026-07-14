#!/bin/bash
for GSE in GSE295601 GSE295704 GSE314596; do
  STUB=$(echo $GSE | sed -E 's/[0-9]{3}$/nnn/')
  LIST=$(curl -s --max-time 20 "https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/" 2>/dev/null | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"//')
  PANEL=$(echo "$LIST" | grep -iE 'feature_ref|antibody|totalseq|abseq|panel|feature.*csv' | grep -ivE 'matrix|barcodes.tsv' | head -1)
  echo "=== $GSE panel=$PANEL ==="
  if [ -n "$PANEL" ]; then
    C=$(curl -s --max-time 25 "https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/$PANEL" 2>/dev/null | zcat 2>/dev/null)
    echo "  n_antibodies: $(echo "$C" | grep -icE 'totalseq|abseq|anti|_ADT|pAbO|CD[0-9]')"
    echo "  CD137-receptor: $(echo "$C" | grep -iE 'CD137|4-1BB|TNFRSF9' | grep -ivE 'ligand|CD137L|TNFSF9|4-1BBL' | head -2 | tr '\n' '|')"
    echo "  CD30: $(echo "$C" | grep -iE '(^|[,\t ])CD30([,\t ]|$)|_CD30_|anti.*CD30' | grep -ivE 'CD300|CD303|CD304|CD305|CD301|CD30L|CD153' | head -1)"
    echo "  isotype-ctrl: $(echo "$C" | grep -iE 'isotype|IgG1|IgG2a|IgG2b|IgG-|mIgG|rIgG|control.*Ig' | head -3 | tr '\n' '|')"
  fi
done
echo DONE

#!/bin/bash
TSV=/media/balthasar-lab/RAID4/costim_engager_counterscreen/analysis/citeseq_rna2protein/cd30_sweep/sweep_results.tsv
echo -e "GSE\tCD30_is_antibody\tpanel_kind"
awk -F'\t' '$5=="YES"{print $1"\t"$2}' $TSV | while IFS=$'\t' read GSE AB; do
  STUB=$(echo $GSE | sed -E 's/[0-9]{3}$/nnn/')
  URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${STUB}/${GSE}/suppl/$AB"
  C=$(curl -s --max-time 20 "$URL" 2>/dev/null | zcat 2>/dev/null)
  # does an Antibody-Capture row name CD30? (TotalSeq 3-col tsv)
  ISAB_TSV=$(echo "$C" | awk -F'\t' 'toupper($3) ~ /ANTIBODY/ && ($2=="CD30"||$2 ~ /^CD30[:_ ]/){f=1} END{print f+0}')
  # BD AbSeq / feature_reference: CD30 named + NO "Gene Expression" rows at all -> pure antibody file
  HAS_GEX=$(echo "$C" | grep -icE 'Gene Expression|Gene_Expression')
  HAS_CD30NAME=$(echo "$C" | grep -icE 'CD30[:|,_ ]|CD30\||^CD30,|,CD30,|Ber-?ACT35|BerH2')
  if [ "$ISAB_TSV" = "1" ]; then echo -e "$GSE\tYES-totalseq\t$AB";
  elif [ "$HAS_GEX" = "0" ] && [ "$HAS_CD30NAME" -gt 0 ]; then echo -e "$GSE\tYES-abseq\t$AB";
  else echo -e "$GSE\tNO-gene-artifact\t$AB"; fi
done
echo "VERIFY_DONE"

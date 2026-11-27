#!/usr/bin/env bash
# Downloads the per-gene XML entries of the Human Protein Atlas v25 for the genes
# of 9_7_0_HPA_genes.tsv (symbol, Ensembl gene ID), pausing between requests.
# Files already present are skipped.
# Output: $OVARY_SIGNATURES/9_Revision_analyses/hpa/xml/<ENSG>.xml, download.log
set -u
BASE="${OVARY_SIGNATURES:-$HOME/Ovary_signatures}"   # project data root
OUT="$BASE/9_Revision_analyses/hpa/xml"
GENES="$(dirname "$0")/9_7_0_HPA_genes.tsv"
UA="Mozilla/5.0 (X11; Linux x86_64) research-script/1.0"
mkdir -p "$OUT"
while read -r gene ensg; do
  f="$OUT/${ensg}.xml"
  [ -s "$f" ] && continue
  code=$(curl -sS -L -A "$UA" -o "$f" -w '%{http_code} %{url_effective}' "https://v25.proteinatlas.org/${ensg}.xml")
  echo -e "$gene\t$ensg\t$code\t$(date -u +%FT%TZ)" | tee -a "$OUT/download.log"
  sleep 3
done < "$GENES"

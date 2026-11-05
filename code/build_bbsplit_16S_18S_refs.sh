#!/bin/bash
# Make small 16S and 18S reference sets for splitting reads of primer pairs that
# amplify both (e.g. 515F-Y/926R) with bbsplit.
#
# 16S and 18S are very different, so for binning reads one sequence per genus is
# plenty and keeps the bbsplit index small:
#   16S: SILVA DADA2 training set, Bacteria and Archaea, plus PR2's organelle
#        16S (plastid, mitochondrion, apicoplast, chromatophore; host labeled, so
#        one per host genus covers far more chloroplast diversity than SILVA's
#        "Chloroplast" lineage).  Organelle reads should stay with the 16S.
#   18S: PR2 DADA2 training set, nuclear (and nucleomorph) 18S only.  PR2 marks
#        organelle sequences by their domain, e.g. Eukaryota:plas
#
# Usage: build_bbsplit_16S_18S_refs.sh <silva_trainset.fa.gz> <pr2_SSU_dada2.fasta.gz> <out_16S.fasta.gz> <out_18S.fasta.gz>
set -euo pipefail

silva=$1
pr2=$2
out_16S=$3
out_18S=$4

# keep the first sequence for each lineage truncated to the given number of ranks
one_per_lineage() {
    local ranks=$1 domains=$2
    seqkit fx2tab \
        | awk -F'\t' -v ranks="$ranks" -v domains="$domains" '
            BEGIN { split(domains, d, ","); for (i in d) keep_domain[d[i]] = 1 }
            {
                n = split($1, tax, ";")
                if (!(tax[1] in keep_domain)) next
                key = tax[1]
                for (i = 2; i <= ranks && i <= n; i++) if (tax[i] != "") key = key ";" tax[i]
                if (key in seen) next
                seen[key] = 1
                printf ">%s_%d %s\n%s\n", substr(tax[1], 1, 3) substr(tax[1], 11, 4), ++count, key, $2
            }' \
        | gzip
}

mkdir -p "$(dirname "$out_16S")" "$(dirname "$out_18S")"
{
    one_per_lineage 6 Bacteria,Archaea < <(zcat "$silva")
    one_per_lineage 8 Eukaryota:plas,Eukaryota:mito,Eukaryota:apic,Eukaryota:chro < <(zcat "$pr2")
} > "$out_16S"
one_per_lineage 8 Eukaryota,Eukaryota:nucl < <(zcat "$pr2") > "$out_18S"

echo "16S references: $(zcat "$out_16S" | grep -c '>')"
echo "18S references: $(zcat "$out_18S" | grep -c '>')"

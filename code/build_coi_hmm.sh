#!/bin/bash
# Build a nucleotide profile HMM covering the full-length COI (cox1) gene.
#
# Sequences come from a MIDORI2 DADA2-format CO1 file (headers are the
# taxonomy: phylum;class;order;family;genus;species;).  Only full-length
# sequences are used so that the model spans the whole gene, and with that the
# binding sites of as many COI primers as possible (Folmer, Leray, Zeale, ...).
# The sampling is balanced across taxa: at most PER_ORDER sequences per order
# and PER_PHYLUM per phylum, so that well sequenced groups (insects,
# vertebrates) don't dominate the model.
#
# Usage: build_coi_hmm.sh <midori2_CO1_DADA2.fasta.gz> <out.hmm> [threads] [per_order] [per_phylum]
#
# Besides <out.hmm> this writes <out>.seqs.fasta (sampled sequences),
# <out>.aln.fasta (MAFFT alignment) and <out>.taxa.tsv (sequence ID to taxonomy).
set -euo pipefail

in_fasta=$1
out_hmm=$2
threads=${3:-8}
per_order=${4:-5}
per_phylum=${5:-150}

MIN_LEN=1400
MAX_LEN=1699
SEED=1

out_base=${out_hmm%.hmm}
mkdir -p "$(dirname "$out_hmm")"

# full-length, (nearly) unambiguous sequences, in random but reproducible order,
# then keep up to per_order per order and per_phylum per phylum
seqkit seq --min-len $MIN_LEN --max-len $MAX_LEN "$in_fasta" \
    | seqkit shuffle --rand-seed $SEED \
    | seqkit fx2tab \
    | awk -F'\t' -v per_order="$per_order" -v per_phylum="$per_phylum" '
        BEGIN { OFS = "\t" }
        {
            seq = toupper($2)
            n_acgt = gsub(/[ACGT]/, "&", seq)
            if (n_acgt < 0.99 * length(seq)) next
            split($1, tax, ";")
            phylum = tax[1]; order = tax[1] ";" tax[2] ";" tax[3]
            if (n_phylum[phylum] >= per_phylum || n_order[order] >= per_order) next
            n_phylum[phylum]++; n_order[order]++
            id = sprintf("COI_%05d", ++n)
            print id, $1 > "/dev/stderr"
            print ">" id "\n" seq
        }' \
    2> "${out_base}.taxa.tsv" \
    > "${out_base}.seqs.fasta"

echo "Sampled $(grep -c '>' "${out_base}.seqs.fasta") sequences from $(cut -f2 "${out_base}.taxa.tsv" | cut -d';' -f1 | sort -u | wc -l) phyla"

mafft --auto --thread "$threads" "${out_base}.seqs.fasta" > "${out_base}.aln.fasta"

# RNA alphabet to match the rRNA models in the combined amplicon HMM database
# (hmmpress refuses mixed alphabets, for nhmmscan searches of reads it makes no difference)
hmmbuild --rna --cpu "$threads" -n COI_mito "$out_hmm" "${out_base}.aln.fasta"

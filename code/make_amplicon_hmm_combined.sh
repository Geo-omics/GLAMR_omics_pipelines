#! /bin/bash
# Build the combined amplicon HMM database used for amplicon target detection.
# Models are listed explicitly (not *.hmm) so that other HMMs in this directory,
# like lgt_516.hmm, are not pulled in.  COI_mito.hmm is made by rule build_coi_hmm.

source ~/miniconda3/bin/activate hmmer

cd ~/GLAMR

hmm_dir=data/reference/hmm_amplicons

rm -f $hmm_dir/combined.hmm*

cat $hmm_dir/arc.hmm \
    $hmm_dir/bac.hmm \
    $hmm_dir/euk.hmm \
    $hmm_dir/mito.hmm \
    $hmm_dir/coi/COI_mito.hmm \
    > $hmm_dir/combined.hmm
hmmpress $hmm_dir/combined.hmm

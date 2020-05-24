#!/bin/bash
#$ -q mic,bio
#$ -N perl
#$ -o sum_cov.out
#$ -e sum_cov.err
#$ -m beas
##$ -ckpt restart

cd $HOME

perl ~/diel_RNAseq/join_count_ann_tables.pl $BIODIR/martiny_diel_seqs/ read_counts.GH.phylodist


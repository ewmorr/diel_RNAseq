#!/bin/bash
#$ -q mic,bio
#$ -N perl
#$ -o sum_cov.out
#$ -e sum_cov.err
#$ -m beas
##$ -ckpt restart

cd $HOME

perl ~/bin/join_count_ann_tables.pl /dfs3/bio/morrise1/martiny_diel_seqs/ read_counts.KO.phylodist
perl ~/bin/join_count_ann_tables.pl /dfs3/bio/morrise1/martiny_diel_seqs/ read_counts.KO
perl ~/bin/join_count_ann_tables.pl /dfs3/bio/morrise1/martiny_diel_seqs/ read_counts.phylodist

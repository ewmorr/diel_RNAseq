#!/bin/bash
#$ -q mic,bio
#$ -N perl
#$ -o sum_cov.out
#$ -e sum_cov.err
#$ -m beas
##$ -ckpt restart

READSDIR=$BIODIR/martiny_diel_seqs
for i in $READSDIR/*metaT
do
	echo $i	
	perl $HOME/diel_RNAseq/sum_read_counts_by_GH_annotation.pl $i/reads_per_contig.txt $BIODIR/P1P2P3_GH_mapping.txt > $i/read_counts.GH
done

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
	perl $HOME/bin/sum_read_counts_by_single_annotation.pl $i/reads_per_contig.txt $BIODIR/martiny_diel_metaG/P1P2P3.KO 2 > $i/read_counts.KO
    perl $HOME/bin/sum_read_counts_by_single_annotation.pl $i/reads_per_contig.txt $BIODIR/martiny_diel_metaG/P1P2P3.phylodist 4 > $i/read_counts.phylodist
done


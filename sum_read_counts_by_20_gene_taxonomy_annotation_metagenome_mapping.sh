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
	perl $HOME/diel_RNAseq/sum_read_counts_by_20_gene_taxonomy_annotation.pl $i/reads_per_contig.txt $BIODIR/metagenome_markers2taxonomy_03192020.txt > $i/read_counts.20_gene_taxonomy.metaG
done

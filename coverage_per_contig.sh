#!/bin/bash
#$ -q mic
#$ -N perl
#$ -o cov_by_seq.out
#$ -e cov_by_seq.err
#$ -m beas

READSDIR=$BIODIR/martiny_diel_seqs
for i in $READSDIR/*metaT
do
	echo $i	
	perl $HOME/assembly_coverage_taxonomic_id/coverage_per_contig_samtools.pl $i/coverage.txt > $i/coverage_by_sequence.txt
done


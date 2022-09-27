#!/bin/bash
#$ -q mic
#$ -N perl
#$ -o sum_cov.out
#$ -e sum_cov.err
#$ -m beas

READSDIR=$BIODIR/martiny_diel_seqs
for i in $READSDIR/*metaT
do
	echo $i	
	echo "phylodist"
	perl $HOME/bin/sum_coverage_by_annotation.pl $i/coverage_by_sequence.txt $BIODIR/martiny_diel_metaG/P1P2P3.phylodist 4 > $i/coverage.phylodist
	echo "KO"
	perl $HOME/bin/sum_coverage_by_annotation.pl $i/coverage_by_sequence.txt $BIODIR/martiny_diel_metaG/P1P2P3.KO 2 > $i/coverage.KO
	echo "COG"
	perl $HOME/bin/sum_coverage_by_annotation.pl $i/coverage_by_sequence.txt $BIODIR/martiny_diel_metaG/P1P2P3.COG 1 > $i/coverage.COG
done


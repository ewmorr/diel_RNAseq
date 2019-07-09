#!/usr/bin/bash
#$ -N usearch
#$ -q mic
#$ -o derep.out
#$ -e derep.err
#$ -m beas
#$ -pe openmp 8
module purge
module load usearch/7.0.1090

cd $BIODIR/martiny_diel_metaG

usearch-64b -derep_fulllength genes_P1P2P3.fna -minseqlength 0 -output genes_P1P2P3.derep.fna

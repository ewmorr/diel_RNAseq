#!/bin/bash
#$ -q mic,bio
#$ -N mapping
#$ -o mapping.out
#$ -e mapping.err
#$ -pe openmp 24
#$ -m beas
module purge
module load bwa/0.7.8
module load samtools/1.8-11

READSDIR=$BIODIR/martiny_diel_seqs

for i in $READSDIR/*_metaT/*.fastq.gz
do
	FILE=${i##*/}
	
	SUBDIR=${i%/*}
	SUBDIR=${SUBDIR##*/}
	echo "starting $SUBDIR"

    SUBJECT=$BIODIR/martiny_diel_metaT_assemblies/$SUBDIR/genes.fna

    bwa index $SUBJECT

	bwa mem -M -t 24 $SUBJECT $i > $READSDIR/$SUBDIR/raw_mapped.sam

	samtools view -@ 24 -Sb -F 4 -o $READSDIR/$SUBDIR/mapped.bam $READSDIR/$SUBDIR/raw_mapped.sam
	samtools sort -@ 24 -T $READSDIR/$SUBDIR/sorted_mapped -o $READSDIR/$SUBDIR/sorted_mapped.bam $READSDIR/$SUBDIR/mapped.bam
	samtools index $READSDIR/$SUBDIR/sorted_mapped.bam
	samtools idxstats $READSDIR/$SUBDIR/sorted_mapped.bam > $READSDIR/$SUBDIR/reads_per_contig_SELF.txt
	rm $READSDIR/$SUBDIR/raw_mapped.sam
	rm $READSDIR/$SUBDIR/mapped.bam
	rm $READSDIR/$SUBDIR/sorted_mapped.bam
	rm $READSDIR/$SUBDIR/sorted_mapped.bam.bai

done


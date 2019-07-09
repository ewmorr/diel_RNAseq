#!/bin/bash
#$ -q mic
#$ -N mapping
#$ -o mapping.out
#$ -e mapping.err
#$ -pe openmp 24
#$ -m beas
module purge
module load bwa/0.7.8
module load samtools/1.8-11

SUBJECT=$BIODIR/martiny_diel_metaG/genes_P1P2P3.derep.fna
READSDIR=$BIODIR/martiny_diel_seqs

bwa index $SUBJECT
#SUBJECTDIR=${SUBJECT%/*}
#perl $HOME/assembly_coverage_taxonomic_id/get_seq_lens.pl $SUBJECT > $SUBJECTDIR/genome.txt

for i in $READSDIR/*_metaT/*.fastq.gz
do
	FILE=${i##*/}
	
	SUBDIR=${i%/*}
	SUBDIR=${SUBDIR##*/}
	echo "starting $SUBDIR"

	bwa mem -M -t 24 $SUBJECT $i > $READSDIR/$SUBDIR/raw_mapped.sam

	samtools view -@ 24 -Sb -F 4 -o $READSDIR/$SUBDIR/mapped.bam $READSDIR/$SUBDIR/raw_mapped.sam
	samtools sort -@ 24 -T $READSDIR/$SUBDIR/sorted_mapped -o $READSDIR/$SUBDIR/sorted_mapped.bam $READSDIR/$SUBDIR/mapped.bam
	samtools index $READSDIR/$SUBDIR/sorted_mapped.bam
	samtools depth -a $READSDIR/$SUBDIR/sorted_mapped.bam > $READSDIR/$SUBDIR/coverage.txt
	samtools view -@ 24 -c -F 4 $READSDIR/$SUBDIR/sorted_mapped.bam > $READSDIR/$SUBDIR/num_mapped_reads.txt
	rm $READSDIR/$SUBDIR/raw_mapped.sam
	rm $READSDIR/$SUBDIR/mapped.bam
	rm $READSDIR/$SUBDIR/sorted_mapped.bam
	rm $READSDIR/$SUBDIR/sorted_mapped.bam.bai
	#perl $HOME/assembly_coverage_taxonomic_id/coverage_per_contig_samtools.pl $READSDIR/$SUBDIR/coverage.txt >> $READSDIR/$SUBDIR/coverage_by_sequence.txt

done


#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="bwa"
#SBATCH --output=bwa_diel.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

SUBJECT=$HOME/martiny_diel_metagenomes/genes_P1P2P3.derep.fna
READSDIR=$HOME/martiny_diel_metagenomes

bwa index $SUBJECT

for i in $READSDIR/*.fastq.gz
do
	FILE=${i##*/}
	SUBDIR=${file%%.filter-METAGENOME.fastq.gz}
    echo "starting $plot"

	srun bwa mem -M -t 24 $SUBJECT $i > $READSDIR/$SUBDIR/raw_mapped.sam

    srun samtools view -@ 24 -Sb -F 4 -o $READSDIR/$SUBDIR/mapped.bam $READSDIR/$SUBDIR/raw_mapped.sam
    srun samtools sort -@ 24 -T $READSDIR/$SUBDIR/sorted_mapped -o $READSDIR/$SUBDIR/sorted_mapped.bam $READSDIR/$SUBDIR/mapped.bam
    srun samtools index $READSDIR/$SUBDIR/sorted_mapped.bam
    srun samtools idxstats $READSDIR/$SUBDIR/sorted_mapped.bam > $READSDIR/$SUBDIR/reads_per_contig.txt

	rm $READSDIR/$SUBDIR/raw_mapped.sam
	rm $READSDIR/$SUBDIR/mapped.bam
	rm $READSDIR/$SUBDIR/sorted_mapped.bam
	rm $READSDIR/$SUBDIR/sorted_mapped.bam.bai

done


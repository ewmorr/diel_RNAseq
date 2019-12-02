#!/usr/bin/bash
#$ -N perl
#$ -q mic
#$ -o extract_genes.out
#$ -e extract_genes.err
#$ -m beas

for i in $BIODIR/martiny_diel_metaT_assemblies/*metax_FD
do
	cd $i
	perl $HOME/bin/extract_gene_seqs_w_gff.pl $i/gene_calls.CDS.gff $i/contigs.fna $i/genes.fna
done	

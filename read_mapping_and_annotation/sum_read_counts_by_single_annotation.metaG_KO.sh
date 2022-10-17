#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="diel_anns"
#SBATCH --output=diel_anns.out
#SBATCH --partition=shared

reads_dir=$HOME/martiny_diel_metagenomes/read_counts
ann_dir=$HOME/martiny_diel_metagenomes
repo_dir=$HOME/repo/diel_RNAseq/read_mapping_and_annotation

#summarize counts over annotations
for i in $reads_dir/*
do

	echo $i
	perl $repo_dir/sum_read_counts_by_single_annotation.pl $i/reads_per_contig.txt \
    $ann_dir/P1P2P3.KO 2 > $i/read_counts.KO.metaG
    
    perl $repo_dir/sum_read_counts_by_GH_annotation.pl $i/reads_per_contig.txt \
    $ann_dir/P1P2P3_GH_mapping.txt > $i/read_counts.GH.metaG
    
    perl $repo_dir/diel_RNAseq/sum_read_counts_by_20_gene_taxonomy_annotation.pl $i/reads_per_contig.txt \
    $ann_dir/metagenome_markers2taxonomy_03192020.txt > $i/read_counts.20_gene_taxonomy.metaG
    
done

#join the tables from the different samples

for i in $reads_dir/*
do

    perl $repo_dir/join_count_ann_tables.pl /mnt/home/garnas/ericm/martiny_diel_metagenomes/read_counts read_counts.KO.metaG
    perl $repo_dir/join_count_ann_tables.pl /mnt/home/garnas/ericm/martiny_diel_metagenomes/read_counts read_counts.GH.metaG
    perl $repo_dir/join_count_ann_tables.pl /mnt/home/garnas/ericm/martiny_diel_metagenomes/read_counts read_counts.20_gene_taxonomy.metaG

done

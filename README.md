### This document describes workflow for processing RNAseq data for assessment of diel patterns in gene expression within microbial communities decomposing grass litter in southern California Mediterranean ecosystem.

#### Illumina shotgun sequencing of 108 eRNA samples (metatranscriptomes) and 3 eDNA samples (metagenomes) was performed by DOE-JGI. Assembly, gene prediction, and taxonomic/functional annotation (phylodist, KEGG, etc) were performed by JGI.

##### Workflow on UCI HPC
This repo is cloned to HPC
##### processing three assembled metagenomes from JGI as target for read mapping
#Rename IMG files; original .gff files are stored as gene_calls.gff
```
cp P1T30_metag/3300028588/3300028588.a.fna P1T30_metag/contigs.fna
cp P1T30_metag/3300028588/3300028588.a.gff P1T30_metag/gene_calls.gff
cp P3T30_metag/3300028586/3300028586.a.fna P3T30_metag/contigs.fna
cp P3T30_metag/3300028586/3300028586.a.gff P3T30_metag/gene_calls.gff
cp P2T30_metag/3300028585/3300028585.a.fna P2T30_metag/contigs.fna
cp P2T30_metag/3300028585/3300028585.a.gff P2T30_metag/gene_calls.gff
```

#Extract only calls for CDS to remove ribosomal sequenecs etc

```
grep "CDS" P1_metaG/gene_calls.gff > P1_metaG/gene_calls.CDS.gff
grep "CDS" P2_metaG/gene_calls.gff > P2_metaG/gene_calls.CDS.gff
grep "CDS" P3_metaG/gene_calls.gff > P3_metaG/gene_calls.CDS.gff
```

#Extract nucleotides sequences of CDS from contigs

```
qsub ~/diel_RNAseq/extract_gff_seqs.sh
```
###### output is written to genes.fna

```cat P1_metaG/genes.fna P2_metaG/genes.fna P3_metaG/genes.fna > genes_P1P2P3.fna```

#Also cat annotations
```
cat $BIODIR/martiny_diel_metaG/P1_metaG/P1T30_metag/*KO $BIODIR/martiny_diel_metaG/P2_metaG/P2T30_metag/*KO $BIODIR/martiny_diel_metaG/P3_metaG/P3T30_metag/*KO > $BIODIR/martiny_diel_metaG/P1P2P3.KO

cat $BIODIR/martiny_diel_metaG/P1_metaG/P1T30_metag/*phylodist $BIODIR/martiny_diel_metaG/P2_metaG/P2T30_metag/*phylodist $BIODIR/martiny_diel_metaG/P3_metaG/P3T30_metag/*phylodist > $BIODIR/martiny_diel_metaG/P1P2P3.phylodist

cat $BIODIR/martiny_diel_metaG/P1_metaG/P1T30_metag/*EC $BIODIR/martiny_diel_metaG/P2_metaG/P2T30_metag/*EC $BIODIR/martiny_diel_metaG/P3_metaG/P3T30_metag/*EC > $BIODIR/martiny_diel_metaG/P1P2P3.EC

```

#Dereplicate gene calls

```
qsub ~/diel_RNAseq/derep_usearch.sh
````

### mapping for read counts. Reads mapped per target length is input for downstream analyses

#output mapped RNA reads per contig with samtools idxstats

```
qsub ~/diel_RNAseq/map_diel_RNA_to_genes_read_counts.sh
```

#Calculate mapped reads by annotations. As above except with read counts per contig not by per base coverage. This will be input to statistical analyses

```
qsub ~/diel_RNAseq/sum_read_counts_by_annotations.sh
qsub ~/diel_RNAseq/sum_read_counts_by_single_annotation.sh
```

#total read count mapped, total length of target reference contigs that were mapped to, and number of genes within an annotation are joined across samples and then reported within three different files (these are the files on Eric's laptop)

```
qsub ~/diel_RNAseq/join_ann_tables.sh
```

#### This same workflow or similar applies for both coverage based and read count based mapping

#Get read mapping totals and total reads

#number of input reads

```
for i in *metaT;do grep "Output" $i/*report.txt | cut -f 2 -d "|" | sed 's/ //g' | sed 's/,//g' > $i/num_input_reads.txt; done
```
#number of mapped reads can be summed from the files `read_per_contig.txt`
```
perl ~/bin/sum_mapped_read_by_sample.pl /dfs3/bio/morrise1/martiny_diel_seqs/ reads_per_contig.txt > mapped_reads_read_counts_112119.txt
```
#Number of input bases

```
for i in *metaT;do grep "Output" $i/*report.txt | cut -f 4 -d "|" | sed 's/ //g' | sed 's/,//g' > $i/num_bases_mapped.txt; done
```

#get metagenome wide coverage totals (length, total bases mapped, avg. cov.). This comes from coverage based mapping and the workflow is no longer included here

```
for i in *metaT;do grep "genome" $i/coverage_by_sequence.txt > $i/total_cov.txt; done
```
#combine to a single file

```
for i in ./*;do echo $i >> totals_mapped_input.txt; cat $i/num_input_reads.txt >> totals_mapped_input.txt; cat $i/num_mapped_reads.txt >> totals_mapped_input.txt; cat $i/num_bases_mapped.txt >> totals_mapped_input.txt; cat $i/total_cov.txt >> totals_mapped_input.txt; done
```

#total counts file can be reformatted with ```reformat_total_counts.pl```

##### Running locally

#Reformat taxonomic strings. The strain field does not play nice in R, and there are special characters in some species fields
```
perl ~/repo/diel_RNAseq/split_and_clean_taxonomic_strings_phylodist.pl read_counts.phylodist.readCount.join > read_counts.phylodist.readCount.join.cleanStrings
perl ~/repo/diel_RNAseq/split_and_clean_taxonomic_strings_phylodist.pl read_counts.phylodist.readCountNumGenes.join > read_counts.phylodist.readCountNumGenes.join.cleanStrings
perl ~/repo/diel_RNAseq/split_and_clean_taxonomic_strings_phylodist.pl read_counts.phylodist.readCountRefLen.join > read_counts.phylodist.readCountRefLen.join.cleanStrings

perl ~/repo/diel_RNAseq/split_and_clean_taxonomic_strings_KO_phylodist.pl read_counts.KO.phylodist.readCount.join > read_counts.KO.phylodist.readCount.join.cleanStrings
perl ~/repo/diel_RNAseq/split_and_clean_taxonomic_strings_KO_phylodist.pl read_counts.KO.phylodist.readCountNumGenes.join > read_counts.KO.phylodist.readCountNumGenes.join.cleanStrings
perl ~/repo/diel_RNAseq/split_and_clean_taxonomic_strings_KO_phylodist.pl read_counts.KO.phylodist.readCountRefLen.join > read_counts.KO.phylodist.readCountRefLen.join.cleanStrings

```

#### read mapping summaries
#Running this interactively
```
read_mapping_stats.r
```

#### delta-delta plots
#First arg is cateogry to summarize counts by, second arg is output file for graphs
#input files are hardcoded in script and currently loaded from working dir

```
Rscript ~/repo/diel_RNAseq/delta_delta_plots_summarize_by_category.r KO delta-delta_plots.KO.pdf
Rscript ~/repo/diel_RNAseq/delta_delta_plots_summarize_by_category.r Genus delta-delta_plots.Genus.pdf
```

#Filter counts table by genus then summarize by KO and run delta-delta for top ten genera
```
topTenGen=(Curtobacterium Massilia Pseudomonas Alternaria Rhizobium Parastagonospora Frigoribacterium Pyrenophora Bipolaris Sphingomonas)
for i in ${topTenGen[@]}
do(
    echo $i
    Rscript ~/repo/diel_RNAseq/delta_delta_plots_filter_then_summarize_by_category.r Genus $i KO delta-delta_plots.${i}-KO.pdf
    
)
done
```

#Run Cyanobacteria
```
Rscript ~/repo/diel_RNAseq/delta_delta_plots_filter_then_summarize_by_category.r Phylum Cyanobacteria KO delta-delta_plots.Cyanobacteria-KO.pdf
```
#Minimum sequence depth is 17...
#modify script to take user defined minDepth
#last arg is desired rarefaction depth. If min depth is below that prints

```
Rscript ~/repo/diel_RNAseq/delta_delta_plots_filter_then_summarize_by_category_minDepthLimit.r Phylum Cyanobacteria KO delta-delta_plots.Cyanobacteria-KO.pdf 1000

    "minimum sequence depth: "                       "17"
    [1] "minimum sequence depth is lower than minimum allowed."
    [1] "Please enter a new rarefaction depth based on the depth per sample below."
```
#also a list of seq depth by sample..
```
100
```


### uploading and processing IMG metatranscriptome assembly data for analysis of taxonomic content
#### the original downloaded data is in `rna_assemblies_IMG_data.tar.gz`

```
tar -xzvf rna_assemblies_IMG_data.tar.gz

#P2T1 has two files per filetype. Using second file (redo of sequencing effort), to avoid probalems with overlapping contig IDs

#rename misnamed sample
for i in rna_assemblies_IMG_data/*/P4T35metax_FD
do(
    path=${i%/*}
    mv $i ./$path/P1T0metax_FD
)
done


#organizing assembled contigs, gene calls, annotations,
for i in $BIODIR/rna_assemblies_IMG_data/rna_assemblies/*/IMG_Data/*.assembled.fna
do(
    file=${i##*/}
    dir=${i%/*}
    cd $dir
    mv $file ../contigs.fna
    cd $BIODIR
    rm $i
    rmdir $dir
)
done

#gff
for i in $BIODIR/rna_assemblies_IMG_data/rna_assemblies_gene_calls/*/IMG_Data/*.assembled.gff
do(
    dir=${i%/*}
    sampleDir=${dir%/IMG_Data}
    sampleDir=${sampleDir##*/}
    mv $i $BIODIR/rna_assemblies_IMG_data/rna_assemblies/$sampleDir/gene_calls.gff
    rmdir $dir
)
done

#names_map
for i in $BIODIR/rna_assemblies_IMG_data/rna_assemblies_names_map/*/IMG_Data/*.assembled.names_map
do(
dir=${i%/*}
sampleDir=${dir%/IMG_Data}
sampleDir=${sampleDir##*/}
mv $i $BIODIR/rna_assemblies_IMG_data/rna_assemblies/$sampleDir/contigs-scf.names_map
rmdir $dir
)
done

#KO
for i in $BIODIR/rna_assemblies_IMG_data/rna_assemblies_KO/*/IMG_Data/*.assembled.KO
do(
dir=${i%/*}
sampleDir=${dir%/IMG_Data}
sampleDir=${sampleDir##*/}
mv $i $BIODIR/rna_assemblies_IMG_data/rna_assemblies/$sampleDir/annotations.KO
rmdir $dir
)
done

#EC
for i in $BIODIR/rna_assemblies_IMG_data/rna_assemblies_EC/*/IMG_Data/*.assembled.EC
do(
dir=${i%/*}
sampleDir=${dir%/IMG_Data}
sampleDir=${sampleDir##*/}
mv $i $BIODIR/rna_assemblies_IMG_data/rna_assemblies/$sampleDir/annotations.EC
rmdir $dir
)
done

mv rna_assemblies_IMG_data/rna_assemblies ./martiny_diel_metaT_assemblies
```

#extract CDS from gff files
```
for i in $BIODIR/martiny_diel_metaT_assemblies/*
do(
    cd $i
    grep "CDS" gene_calls.gff > gene_calls.CDS.gff
)
done
```

#extract nuc acid seqs
```
qsub extract_gff_seqs_metatranscriptome_assemblies.sh
```

#### Mapping reads to self assembly. This is for Alex's tax ID mapping. Will use annotation of genomes against tax DB with original mapping
```
qsub ~/diel_RNAseq/map_diel_RNA_to_SELF_asmb_genes_read_counts.sh



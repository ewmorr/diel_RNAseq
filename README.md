### This document derscribes workflow for processing RNAseq data for assessment of diel patterns in gene expression within microbial communities decomposing grass litter in southern California Mediterranean ecosystem.

#### Illumina shotgun sequencing of 108 eRNA samples (metatranscriptomes) and 3 eDNA samples (metagenomes) was performed by DOE-JGI. Assembly, gene prediction, and taxonomic/functional annotation (phylodist, KEGG, etc) were performed by JGI.

##### Workflow on UCI HPC

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
qsub ~/batch_scripts/extract_gff_seqs.sh
```
###### output is written to $DIR/genes.fna

```cat P1_metaG/genes.fna P2_metaG/genes.fna P3_metaG/genes.fna > genes_P1P2P3.fna```

#Dereplicate gene calls

```
qsub ~/batch_scripts/derep_usearch.sh
````

#perform mapping of RNA reads. This script outputs mapped RNA reads per contig with samtools idxstats

```
qsub ~/batch_scripts/map_diel_RNA_to_genes.sh
```
#Calculate coverage

```
qsub ~/batch_scripts/coverage_per_contig.sh
```

#Calculate coverage of different annotations; this is the total coverage, and average per base coverage of each gene (from coverage_per_contig.pl) then summed by annotations, as well as the number of annotations; this has been performed for phylodist, KEGG (KO), and COG;

```
qsub ~/batch_scripts/sum_coverage_by_annotation.sh
```

#total base coverage, per length (base) average coverage, and number of genes within an annotation are joined across samples and then reported wihtin three different files (these are the files on Eric's laptop)

```
perl ~/bin/join_ann_tables.pl /bio/morrise1/martiny_diel_seqs/ coverage.COG
```

### Rerunning mapping for read counts instead of coverage as this will be more appropriate metric

#rerun mapping of RNA reads to get read counts instead of coverage. This script outputs mapped RNA reads per contig with samtools idxstats

```
qsub ~/batch_scripts/map_diel_RNA_to_genes_read_counts.sh
```

#Calculate mapped reads by annotations. As above except with read counts per contig not by per base coverage. This will be input to statistical analyses

```
qsub ~/batch_scripts/sum_read_counts_by_annotations.sh
```

#total read count mapped, total length of target reference contigs that were mapped to, and number of genes within an annotation are joined across samples and then reported wihtin three different files (these are the files on Eric's laptop)

```
perl ~/bin/join_count_ann_tables.pl /dfs3/bio/morrise1/martiny_diel_seqs/ read_counts.KO.phylodist
```

#### This same workflow or similar applies for both coverage based and read count based mapping

#Get read mapping totals and total reads

#number of input reads

```
for i in *metaT;do grep "Output" $i/*report.txt | cut -f 2 -d "|" | sed 's/ //g' | sed 's/,//g' > $i/num_input_reads.txt; done
```
#number of mapped reads comes from coverage_per_contig.sh. Or this can be summed from the files `read_per_contig.txt`
```
perl ~/bin/sum_mapped_read_by_sample.pl /dfs3/bio/morrise1/martiny_diel_seqs/ reads_per_contig.txt > mapped_reads_read_counts_112119.txt
```
#Number of input bases

```
for i in *metaT;do grep "Output" $i/*report.txt | cut -f 4 -d "|" | sed 's/ //g' | sed 's/,//g' > $i/num_bases_mapped.txt; done
```

#get metagenome wide coverage totals (length, total bases mapped, avg. cov.)

```
for i in *metaT;do grep "genome" $i/coverage_by_sequence.txt > $i/total_cov.txt; done
```
#combine to a single file

```
for i in ./*;do echo $i >> totals_mapped_input.txt; cat $i/num_input_reads.txt >> totals_mapped_input.txt; cat $i/num_mapped_reads.txt >> totals_mapped_input.txt; cat $i/num_bases_mapped.txt >> totals_mapped_input.txt; cat $i/total_cov.txt >> totals_mapped_input.txt; done
```

#total counts file can be reformatted with ```reformat_total_counts.pl```

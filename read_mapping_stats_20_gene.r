require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)

source("../ggplot_theme.txt")

fancy_scientific <- function(l) {
	# turn in to character string in scientific notation
	l <- format(l, scientific = TRUE)
	#reformat zeros
	l <- gsub("0e\\+00","0",l)
	# quote the part before the exponent to keep all the digits
	l <- gsub("^(.*)e", "'\\1'e", l)
	# turn the 'e+' into plotmath format
	l <- gsub("e", "%*%10^", l)
	# return this as an expression
	parse(text=l)
}

#total read amounts table
total_mapped = read.csv("totals_mapped_input.reformat.txt", sep = "\t", header = F)
colnames(total_mapped) = c("Sample", "totalReads", "readsMappedCov", "totalBases", "targetsLength", "basesMapped", "avgCov")
total_mapped = total_mapped %>% as.tbl
total_mapped = total_mapped %>% filter(Sample != "P2T31_metaT")

#metadata table
metadata = read.table("metadata.txt", header = T, sep = "\t")
metadata = metadata %>% as.tbl
metadata = metadata %>% filter(Sample != "P2T31_metaT")

mapped.metadata = full_join(total_mapped, metadata, by = "Sample")

#RNA seq data, read_counts is the number of reads mapped per category (e.g. taxonomic string)

#import
read_count.20_gene_taxonomy = read.table("read_counts.20_gene_taxonomy.readCount.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
read_count.20_gene_taxonomy$P2T31_metaT = NULL
read_count.20_gene_taxonomy$X = NULL
colSums(read_count.20_gene_taxonomy[,5:111]) %>% sort

colSums(read_count.20_gene_taxonomy[,5:111])/total_mapped$totalReads



#Join annotated counts to count table & plot
total_mapped = left_join(total_mapped, data.frame(Sample = names(colSums(read_count.20_gene_taxonomy[,5:111])), twenty_gene_taxonomy.annotated = colSums(read_count.20_gene_taxonomy[,5:111]) ))

plot_read_mapping.df = gather(total_mapped, key = "readCat", value = "reads", -Sample) %>%
filter(readCat == "totalReads" | readCat == "twenty_gene_taxonomy.annotated")

plot_read_mapping.df$readCat = factor(as.factor(plot_read_mapping.df$readCat), levels = c("totalReads", "twenty_gene_taxonomy.annotated"))
plot_read_mapping.df = plot_read_mapping.df[order(plot_read_mapping.df$reads), ]
names_for_facets = list("totalReads" = "total reads", "twenty_gene_taxonomy.annotated" = "20 gene taxonomy")

pdf("reads_counts_total_mapped_annotated_20_gene.pdf", width = 8, height = 6)
ggplot(plot_read_mapping.df, aes(Sample, reads) ) +
geom_point() +
facet_wrap(~readCat, labeller = labeller(readCat = names_for_facets)) +
scale_y_log10(labels = fancy_scientific) +
labs(x = "sample") +
my_gg_theme +
theme(axis.text.x = element_blank())
dev.off()

#RNA seq data, ref_len is the total length of reference contigs per category
ref_len.20_gene_taxonomy = read.table("read_counts.20_gene_taxonomy.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len.20_gene_taxonomy$P2T31_metaT = NULL
ref_len.20_gene_taxonomy$X = NULL

num_genes.20_gene_taxonomy = read.table("read_counts.20_gene_taxonomy.readCountNumGenes.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
num_genes.20_gene_taxonomy$P2T31_metaT = NULL
num_genes.20_gene_taxonomy$X = NULL


#use min read count for rarefaction depth
minDepth = (colSums(read_count.20_gene_taxonomy[,5:111]) %>% sort)[1]
read_count.rarefied = data.frame(read_count.20_gene_taxonomy[,1:4], t(rrarefy(t(read_count.20_gene_taxonomy[,5:111]), sample = minDepth)) )

#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len.20_gene_taxonomy = data.frame(read_count.rarefied[,1:4], read_count.rarefied[,5:111]/ref_len.20_gene_taxonomy[,5:111] ) %>% as.tbl
reads_per_len.20_gene_taxonomy[is.na(reads_per_len.20_gene_taxonomy)] = 0

#look at totals by Tax_name
reads_per_len.Tax_name = reads_per_len.20_gene_taxonomy %>% group_by(Tax_name) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
reads_per_len.Tax_name.rowSum = data.frame(Tax_name = reads_per_len.Tax_name$Tax_name, sum = reads_per_len.Tax_name[2:108] %>% rowSums)
reads_per_len.Tax_name.rowSum[order(reads_per_len.Tax_name.rowSum$sum),]
reads_per_len.Tax_name.rowSum %>% filter(Tax_name == "Curtobacterium")

(reads_per_len.Tax_name.rowSum[order(reads_per_len.Tax_name.rowSum$sum, decreasing = T),])[1:10,]


#raw read counts by Tax_name (not len normalized)
read_count.Tax_name = read_count.rarefied %>% group_by(Tax_name) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
read_count.Tax_name.rowSum = data.frame(Tax_name = read_count.Tax_name$Tax_name, sum = read_count.Tax_name[2:108] %>% rowSums)
read_count.Tax_name.rowSum[order(read_count.Tax_name.rowSum$sum),]
sum(read_count.Tax_name.rowSum$sum)

read_count.Tax_name.rowSum %>% filter(Tax_name == "Curtobacterium")

ref_len.Tax_name = ref_len.20_gene_taxonomy %>% group_by(Tax_name) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
ref_len.Tax_name.rowSum = data.frame(Tax_name = ref_len.Tax_name$Tax_name, sum = ref_len.Tax_name[2:108] %>% rowSums)
ref_len.Tax_name.rowSum[order(ref_len.Tax_name.rowSum$sum),]


#Plot top 10 genera by different metrics

top_ten_genera.reads_per_target =
    data.frame(Category = "20 gene taxonomy", (reads_per_len.Tax_name.rowSum[order(reads_per_len.Tax_name.rowSum$sum, decreasing = T),])[1:10,])


top_ten_genera.total_reads =
data.frame(Category = "20 gene taxonomy", (read_count.Tax_name.rowSum[order(read_count.Tax_name.rowSum$sum, decreasing = T),])[1:10,])

top_ten_genera.reads_per_target = top_ten_genera.reads_per_target  %>%
arrange(Category, -sum) %>%
mutate(order = row_number())

p1 = ggplot(top_ten_genera.reads_per_target,
    aes(x = order, y = sum)) +
geom_col() +
facet_wrap(~Category, scales = "free_x") +
labs(y = "reads mapped\nper target length", x = "") +
scale_x_continuous(breaks = top_ten_genera.reads_per_target$order, labels = as.character(top_ten_genera.reads_per_target$Tax_name)) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))


top_ten_genera.total_reads = top_ten_genera.total_reads  %>%
arrange(Category, -sum) %>%
mutate(order = row_number())

p2 = ggplot(top_ten_genera.total_reads, aes(x = order, y = sum)) +
geom_col() +
facet_wrap(~Category, scales = "free_x") +
scale_y_continuous(labels = fancy_scientific) +
scale_x_continuous(breaks = top_ten_genera.total_reads$order, labels = as.character(top_ten_genera.total_reads$Tax_name)) +
labs(y = "reads mapped", x = "") +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("top_ten_genera.20_gene_tax.pdf", width = 12, height = 10)
grid.arrange(p1,p2,nrow = 2)
dev.off()

#RAREFACTION SHOULD BE PERFORMED ON READ_COUNTS TABLE BEFORE BEING NORMALIZED BY LEN. THIS WILL ALLOW RAREFACTION BY READ COUNTS INSTEAD OF ON A 'FUDGED' METRIC (e.g., rounding of length normalized counts). FOR ANY DOWNSTREAM FILTERING (E.G., BY Tax_name) CAN PERFORM FILTERING ON COUNT AND LEN TABLES, DO RAREFACTION ON FILTERED COUNTS, AND THEN NORMALIZE.

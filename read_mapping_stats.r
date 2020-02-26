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

reads_mapped.count_based = read.table("mapped_reads_read_counts_112119.txt", header = F)
colnames(reads_mapped.count_based) = c("Sample", "readsMappedCount") #This is second mapping for counts
total_mapped = left_join(total_mapped, reads_mapped.count_based, by = "Sample")

total_mapped$readsMappedCov/total_mapped$readsMappedCount
#read mapping peforme the same

#metadata table
metadata = read.table("metadata.txt", header = T, sep = "\t")
metadata = metadata %>% as.tbl
metadata = metadata %>% filter(Sample != "P2T31_metaT")

mapped.metadata = full_join(total_mapped, metadata, by = "Sample")

#RNA seq data, read_counts is the number of reads mapped per category (e.g. taxonomic string)

#KO and phylodist
read_count.KO.phylodist = read.table("read_counts.KO.phylodist.readCount.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
read_count.KO.phylodist$P2T31_metaT = NULL
read_count.KO.phylodist$X = NULL
colSums(read_count.KO.phylodist[,9:115]) %>% sort

colSums(read_count.KO.phylodist[,9:115])/total_mapped$totalReads
colSums(read_count.KO.phylodist[,9:115])/total_mapped$readsMappedCount

#KO
read_count.KO = read.table("read_counts.KO.readCount.join", header = T, sep = "\t") %>% as.tbl
read_count.KO$P2T31_metaT = NULL
read_count.KO$X = NULL
colSums(read_count.KO[,2:108]) %>% sort

(colSums(read_count.KO[,2:108]) %>% sort)/(colSums(read_count.KO.phylodist[,9:115]) %>% sort)
#KO and KO by phylodist has same number of mapped reads annotated

colSums(read_count.KO[,2:108])/total_mapped$totalReads
colSums(read_count.KO[,2:108])/total_mapped$readsMappedCount

#phylodist

read_count.phylodist = read.table("read_counts.phylodist.readCount.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
read_count.phylodist$P2T31_metaT = NULL
read_count.phylodist$X = NULL
colSums(read_count.phylodist[,8:114]) %>% sort
(colSums(read_count.KO[,2:108]) %>% sort)/(colSums(read_count.phylodist[,8:114]) %>% sort)
#KO has 40-50% of reads mapped compared to phylodist

colSums(read_count.phylodist[,8:114])/total_mapped$totalReads
colSums(read_count.phylodist[,8:114])/total_mapped$readsMappedCount


#Join annotated counts to count table & plot
total_mapped = left_join(total_mapped, data.frame(Sample = names(colSums(read_count.phylodist[,8:114])), phylodist.annotated = colSums(read_count.phylodist[,8:114]) ))
total_mapped = left_join(total_mapped, data.frame(Sample = names(colSums(read_count.KO.phylodist[,9:115])), KO.phylodist.annotated = colSums(read_count.KO.phylodist[,9:115]) ))

plot_read_mapping.df = gather(total_mapped, key = "readCat", value = "reads", -Sample) %>%
filter(readCat == "totalReads" | readCat == "readsMappedCount" | readCat == "phylodist.annotated" | readCat == "KO.phylodist.annotated")

plot_read_mapping.df$readCat = factor(as.factor(plot_read_mapping.df$readCat), levels = c("totalReads", "readsMappedCount", "phylodist.annotated", "KO.phylodist.annotated"))
plot_read_mapping.df = plot_read_mapping.df[order(plot_read_mapping.df$reads), ]
names_for_facets = list("totalReads" = "total reads", "readsMappedCount" = "reads mapped", "phylodist.annotated" = "phylodist", "KO.phylodist.annotated" = "KO & phylodist")

pdf("reads_counts_total_mapped_annotated.pdf", width = 8, height = 6)
ggplot(plot_read_mapping.df, aes(Sample, reads) ) +
geom_point() +
facet_wrap(~readCat, labeller = labeller(readCat = names_for_facets)) +
scale_y_log10(labels = fancy_scientific) +
labs(x = "sample") +
my_gg_theme +
theme(axis.text.x = element_blank())
dev.off()

#RNA seq data, ref_len is the total length of reference contigs per category
ref_len.KO.phylodist = read.table("read_counts.KO.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len.KO.phylodist$P2T31_metaT = NULL
ref_len.KO.phylodist$X = NULL

num_genes.KO.phylodist = read.table("read_counts.KO.phylodist.readCountNumGenes.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
num_genes.KO.phylodist$P2T31_metaT = NULL
num_genes.KO.phylodist$X = NULL


ref_len.phylodist = read.table("read_counts.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len.phylodist$P2T31_metaT = NULL
ref_len.phylodist$X = NULL

minDepth = (colSums(read_count.KO.phylodist[,9:115]) %>% sort)[1]
read_count.KO.phylodist.rarefied = data.frame(read_count.KO.phylodist[,1:8], t(rrarefy(t(read_count.KO.phylodist[,9:115]), sample = minDepth)) )

minDepth = (colSums(read_count.phylodist[,8:114]) %>% sort)[1]
read_count.phylodist.rarefied = data.frame(read_count.phylodist[,1:7], t(rrarefy(t(read_count.phylodist[,8:114]), sample = minDepth)) )

#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len.KO.phylodist = data.frame(read_count.KO.phylodist.rarefied[,1:8], read_count.KO.phylodist.rarefied[,9:115]/ref_len.KO.phylodist[,9:115] ) %>% as.tbl
reads_per_len.KO.phylodist[is.na(reads_per_len.KO.phylodist)] = 0

reads_per_len.phylodist = data.frame(read_count.phylodist.rarefied[,1:7], read_count.phylodist.rarefied[,8:114]/ref_len.phylodist[,8:114] ) %>% as.tbl
reads_per_len.phylodist[is.na(reads_per_len.phylodist)] = 0

#look at totals by Genus
reads_per_len.Genus = reads_per_len.KO.phylodist %>% group_by(Genus) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
reads_per_len.Genus.rowSum = data.frame(Genus = reads_per_len.Genus$Genus, sum = reads_per_len.Genus[2:108] %>% rowSums)
reads_per_len.Genus.rowSum[order(reads_per_len.Genus.rowSum$sum),]
reads_per_len.Genus.rowSum %>% filter(Genus == "Curtobacterium")

(reads_per_len.Genus.rowSum[order(reads_per_len.Genus.rowSum$sum, decreasing = T),])[1:10,]

reads_per_len.phylodist.Genus = reads_per_len.phylodist %>% group_by(Genus) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
reads_per_len.phylodist.Genus.rowSum = data.frame(Genus = reads_per_len.phylodist.Genus$Genus, sum = reads_per_len.phylodist.Genus[2:108] %>% rowSums)
reads_per_len.phylodist.Genus.rowSum[order(reads_per_len.phylodist.Genus.rowSum$sum),]
reads_per_len.phylodist.Genus.rowSum %>% filter(Genus == "Curtobacterium")
(reads_per_len.phylodist.Genus.rowSum[order(reads_per_len.phylodist.Genus.rowSum$sum, decreasing = T),])[1:10,]
#Curtobacterium is top in KO.phylodist, Pseudomonas is top in phylodist


#raw read counts by genus (not len normalized)
read_count.Genus = read_count.KO.phylodist.rarefied %>% group_by(Genus) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
read_count.Genus.rowSum = data.frame(Genus = read_count.Genus$Genus, sum = read_count.Genus[2:108] %>% rowSums)
read_count.Genus.rowSum[order(read_count.Genus.rowSum$sum),]
sum(read_count.Genus.rowSum$sum)

read_count.Genus.rowSum %>% filter(Genus == "Curtobacterium")

ref_len.Genus = ref_len.KO.phylodist %>% group_by(Genus) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
ref_len.Genus.rowSum = data.frame(Genus = ref_len.Genus$Genus, sum = ref_len.Genus[2:108] %>% rowSums)
ref_len.Genus.rowSum[order(ref_len.Genus.rowSum$sum),]


read_count.phylodist.Genus = read_count.phylodist.rarefied %>% group_by(Genus) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
read_count.phylodist.Genus.rowSum = data.frame(Genus = read_count.phylodist.Genus$Genus, sum = read_count.phylodist.Genus[2:108] %>% rowSums)
read_count.phylodist.Genus.rowSum[order(read_count.phylodist.Genus.rowSum$sum),]
sum(read_count.phylodist.Genus.rowSum$sum)

ref_len.phylodist.Genus = ref_len.phylodist %>% group_by(Genus) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
ref_len.phylodist.Genus.rowSum = data.frame(Genus = ref_len.phylodist.Genus$Genus, sum = ref_len.phylodist.Genus[2:108] %>% rowSums)
ref_len.phylodist.Genus.rowSum[order(ref_len.phylodist.Genus.rowSum$sum),]


#Plot top 10 genera by different metrics

top_ten_genera.reads_per_target = rbind(
    data.frame(Category = "KO & phylodist", (reads_per_len.Genus.rowSum[order(reads_per_len.Genus.rowSum$sum, decreasing = T),])[1:10,]),
    data.frame(Category = "phylodist", (reads_per_len.phylodist.Genus.rowSum[order(reads_per_len.phylodist.Genus.rowSum$sum, decreasing = T),])[1:10,])
)

top_ten_genera.total_reads = rbind(
data.frame(Category = "KO & phylodist", (read_count.Genus.rowSum[order(read_count.Genus.rowSum$sum, decreasing = T),])[1:10,]),
data.frame(Category = "phylodist", (read_count.phylodist.Genus.rowSum[order(read_count.phylodist.Genus.rowSum$sum, decreasing = T),])[1:10,])
)

top_ten_genera.reads_per_target = top_ten_genera.reads_per_target  %>%
arrange(Category, -sum) %>%
mutate(order = row_number())

p1 = ggplot(top_ten_genera.reads_per_target,
    aes(x = order, y = sum)) +
geom_col() +
facet_wrap(~Category, scales = "free_x") +
labs(y = "reads mapped\nper target length", x = "") +
scale_x_continuous(breaks = top_ten_genera.reads_per_target$order, labels = as.character(top_ten_genera.reads_per_target$Genus)) +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))


top_ten_genera.total_reads = top_ten_genera.total_reads  %>%
arrange(Category, -sum) %>%
mutate(order = row_number())

p2 = ggplot(top_ten_genera.total_reads, aes(x = order, y = sum)) +
geom_col() +
facet_wrap(~Category, scales = "free_x") +
scale_y_continuous(labels = fancy_scientific) +
scale_x_continuous(breaks = top_ten_genera.total_reads$order, labels = as.character(top_ten_genera.total_reads$Genus)) +
labs(y = "reads mapped", x = "") +
my_gg_theme +
theme(axis.text.x = element_text(angle = 45, hjust = 1))


pdf("top_ten_genera.pdf", width = 12, height = 10)
grid.arrange(p1,p2,nrow = 2)
dev.off()

#RAREFACTION SHOULD BE PERFORMED ON READ_COUNTS TABLE BEFORE BEING NORMALIZED BY LEN. THIS WILL ALLOW RAREFACTION BY READ COUNTS INSTEAD OF ON A 'FUDGED' METRIC (e.g., rounding of length normalized counts). FOR ANY DOWNSTREAM FILTERING (E.G., BY GENUS) CAN PERFORM FILTERING ON COUNT AND LEN TABLES, DO RAREFACTION ON FILTERED COUNTS, AND THEN NORMALIZE.

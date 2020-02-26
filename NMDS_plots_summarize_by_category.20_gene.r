require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)
require(reshape2)
require(rlang)

source("~/ggplot_theme.txt")

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

#arguments for summarizing and output
args = commandArgs(T)
catToSum = args[1]
outputFile = args[2]

#total read amounts table
total_mapped = read.csv("totals_mapped_input.reformat.txt", sep = "\t", header = F)
colnames(total_mapped) = c("Sample", "totalReads", "readsMappedCov", "totalBases", "targetsLength", "basesMapped", "avgCov")
total_mapped = total_mapped %>% as.tbl
total_mapped = total_mapped %>% filter(Sample != "P2T31_metaT")

reads_mapped.count_based = read.table("mapped_reads_read_counts_112119.txt", header = F)
colnames(reads_mapped.count_based) = c("Sample", "readsMappedCount") #This is second mapping for counts
total_mapped = left_join(total_mapped, reads_mapped.count_based, by = "Sample")

#metadata table
metadata = read.table("metadata.txt", header = T, sep = "\t")
metadata = metadata %>% as.tbl
metadata = metadata %>% filter(Sample != "P2T31_metaT")

mapped.metadata = full_join(total_mapped, metadata, by = "Sample")

#RNA seq data, read_counts is the number of reads mapped per category (e.g. taxonomic string)

#KO and phylodist
read_count = read.table("read_counts.20_gene_taxonomy.readCount.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
read_count$P2T31_metaT = NULL
read_count$X = NULL
read_count$Edge_num = as.factor(read_count$Edge_num)

#RNA seq data, ref_len is the total length of reference contigs per category
ref_len = read.table("read_counts.20_gene_taxonomy.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL
ref_len$Edge_num = as.factor(ref_len$Edge_num)

#RAREFACTION SHOULD BE PERFORMED ON READ_COUNTS TABLE BEFORE BEING NORMALIZED BY LEN. THIS WILL ALLOW RAREFACTION BY READ COUNTS INSTEAD OF ON A 'FUDGED' METRIC (e.g., rounding of length normalized counts). FOR ANY DOWNSTREAM FILTERING (E.G., BY GENUS) CAN PERFORM FILTERING ON COUNT AND LEN TABLES, DO RAREFACTION ON FILTERED COUNTS, AND THEN NORMALIZE.

#use min read count for rarefaction depth
minDepth = (colSums(read_count[,5:111]) %>% sort)[1]
print(c("minimum sequence depth: ", minDepth))
#Rarefy then normalize by target len and summarize by KO for KO delta-delta analysis

read_count.rarefied = data.frame(read_count[,1:4], t(rrarefy(t(read_count[,5:111]), sample = minDepth)) )

#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
read_count.rarefied.sumCat = read_count.rarefied %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
ref_len.sumCat = ref_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

print("summarized count tbl:")
print(read_count.rarefied.sumCat)
#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len = data.frame(read_count.rarefied.sumCat[,1], read_count.rarefied.sumCat[,2:108]/ref_len.sumCat[,2:108] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

#Full dataset NMDS

reads_per_len.nmds = metaMDS(
sqrt(
t(reads_per_len[
#The indexing term filters outs that only appear in one sample
rowSums(reads_per_len[,2:length(colnames(reads_per_len))] > 0) >= 2,
2:length(colnames(reads_per_len))
])
),
distance = "bray", k = 2
)
while(reads_per_len.nmds$converged != TRUE){
    reads_per_len.nmds = metaMDS(sqrt(
    t(reads_per_len[
    #The indexing term filters outs that only appear in one sample
    rowSums(reads_per_len[,2:length(colnames(reads_per_len))] > 0) >= 2,
    2:length(colnames(reads_per_len))
    ])
    ),
    distance = "bray", k = 2, previous.best = reads_per_len.nmds)
}

reads_per_len.nmds.metadata = left_join(
    data.frame(Sample = rownames(reads_per_len.nmds$points), MDS1 = reads_per_len.nmds$points[,1], MDS2 = reads_per_len.nmds$points[,2]),
    mapped.metadata
)

require(RColorBrewer)

p1 = ggplot(reads_per_len.nmds.metadata, aes(MDS1, MDS2, shape = Plot, color = hours.cumulative.RNA)) +
geom_point(size = 3) +
my_gg_theme +
labs(color = "Time (hrs)") +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
theme(legend.title = element_text(size = 22))

p2 = ggplot(filter(reads_per_len.nmds.metadata, Moisture >= 0), aes(MDS1, MDS2, shape = Plot, color = Moisture)) +
geom_point(size = 3) +
my_gg_theme +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "sqrt") +
labs(color = "Moisture") +
theme(legend.title = element_text(size = 22))

p3 = ggplot(reads_per_len.nmds.metadata, aes(MDS1, MDS2, shape = Plot, color = Temperature)) +
geom_point(size = 3) +
my_gg_theme +
labs(color = "Temperature") +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
theme(legend.title = element_text(size = 22))

pdf("20_gene_NMDS.pdf", width = 24, height = 6)
grid.arrange(p1,p2,p3, ncol = 3)
dev.off()

#transform DF

reads_per_len.long.meta =
#long form to add metadata
reads_per_len %>% gather(key = Sample, value = count, -sumCat) %>%
#metadata to filter
left_join(., mapped.metadata, by = "Sample")

#filter by plots and spread for distance calculation

reads_per_len.wide.P1 = reads_per_len.long.meta %>%
filter(Plot == "P1") %>%
select(c("sumCat", "Sample", "count")) %>%
spread(., Sample, count)

reads_per_len.wide.P2 = reads_per_len.long.meta %>%
filter(Plot == "P2") %>%
select(c("sumCat", "Sample", "count")) %>%
spread(., Sample, count)

reads_per_len.wide.P3 = reads_per_len.long.meta %>%
filter(Plot == "P3") %>%
select(c("sumCat", "Sample", "count")) %>%
spread(., Sample, count)

#NMDS by plot

reads_per_len.wide.P1.nmds = metaMDS(sqrt(
    t(reads_per_len.wide.P1[
    #The indexing term filters outs that only appear in one sample
    rowSums(reads_per_len[,2:length(colnames(reads_per_len.wide.P1))] > 0) >= 2,
    2:length(colnames(reads_per_len.wide.P1))
    ])
),
distance = "bray", k = 2
)
while(reads_per_len.wide.P1.nmds$converged != TRUE){
    reads_per_len.wide.P1.nmds = metaMDS(sqrt(
    t(reads_per_len.wide.P1[
    #The indexing term filters outs that only appear in one sample
    rowSums(reads_per_len.wide.P1[,2:length(colnames(reads_per_len.wide.P1))] > 0) >= 2,
    2:length(colnames(reads_per_len.wide.P1))
    ])
    ),
    distance = "bray", k = 2, previous.best = reads_per_len.wide.P1.nmds)
}

reads_per_len.wide.P2.nmds = metaMDS(sqrt(
t(reads_per_len.wide.P2[
#The indexing term filters outs that only appear in one sample
rowSums(reads_per_len[,2:length(colnames(reads_per_len.wide.P2))] > 0) >= 2,
2:length(colnames(reads_per_len.wide.P2))
])
),
distance = "bray", k = 2
)
while(reads_per_len.wide.P2.nmds$converged != TRUE){
    reads_per_len.wide.P2.nmds = metaMDS(sqrt(
    t(reads_per_len.wide.P2[
    #The indexing term filters outs that only appear in one sample
    rowSums(reads_per_len.wide.P2[,2:length(colnames(reads_per_len.wide.P2))] > 0) >= 2,
    2:length(colnames(reads_per_len.wide.P2))
    ])
    ),
    distance = "bray", k = 2, previous.best = reads_per_len.wide.P2.nmds)
}


reads_per_len.wide.P3.nmds = metaMDS(sqrt(
t(reads_per_len.wide.P3[
#The indexing term filters outs that only appear in one sample
rowSums(reads_per_len[,2:length(colnames(reads_per_len.wide.P3))] > 0) >= 2,
2:length(colnames(reads_per_len.wide.P3))
])
),
distance = "bray", k = 2
)
while(reads_per_len.wide.P3.nmds$converged != TRUE){
    reads_per_len.wide.P3.nmds = metaMDS(sqrt(
    t(reads_per_len.wide.P3[
    #The indexing term filters outs that only appear in one sample
    rowSums(reads_per_len.wide.P3[,2:length(colnames(reads_per_len.wide.P3))] > 0) >= 2,
    2:length(colnames(reads_per_len.wide.P3))
    ])
    ),
    distance = "bray", k = 2, previous.best = reads_per_len.wide.P3.nmds)
}

#Add metadata for plotting

reads_per_len.nmds.P1.metadata = left_join(
data.frame(Sample = rownames(reads_per_len.wide.P1.nmds$points), MDS1 = reads_per_len.wide.P1.nmds$points[,1], MDS2 = reads_per_len.wide.P1.nmds$points[,2]),
mapped.metadata
)

reads_per_len.nmds.P2.metadata = left_join(
data.frame(Sample = rownames(reads_per_len.wide.P2.nmds$points), MDS1 = reads_per_len.wide.P2.nmds$points[,1], MDS2 = reads_per_len.wide.P2.nmds$points[,2]),
mapped.metadata
)

reads_per_len.nmds.P3.metadata = left_join(
data.frame(Sample = rownames(reads_per_len.wide.P3.nmds$points), MDS1 = reads_per_len.wide.P3.nmds$points[,1], MDS2 = reads_per_len.wide.P3.nmds$points[,2]),
mapped.metadata
)

reads_per_len.nmds.P1.P2.P3.metadata = rbind(reads_per_len.nmds.P1.metadata, reads_per_len.nmds.P2.metadata, reads_per_len.nmds.P3.metadata)

p1 = ggplot(reads_per_len.nmds.P1.P2.P3.metadata, aes(MDS1, MDS2, color = hours.cumulative.RNA)) +
geom_point(size = 2) +
facet_wrap(~Plot, ncol = 3, scales = "free") +
my_gg_theme +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
labs(color = "Time (hrs)") +
theme(legend.title = element_text(size = 22))

p2 = ggplot(filter(reads_per_len.nmds.P1.P2.P3.metadata, Moisture >= 0), aes(MDS1, MDS2, color = Moisture)) +
geom_point(size = 2) +
facet_wrap(~Plot, ncol = 3, scales = "free") +
my_gg_theme +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "sqrt") +
labs(color = "Moisture") +
theme(legend.title = element_text(size = 22))

p3 = ggplot(reads_per_len.nmds.P1.P2.P3.metadata, aes(MDS1, MDS2, color = Temperature)) +
geom_point(size = 2) +
facet_wrap(~Plot, ncol = 3, scales = "free") +
my_gg_theme +
labs(color = "Temperature") +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
theme(legend.title = element_text(size = 22))

pdf("20_gene_NMDS_by_plot.pdf", width = 12, height = 8)
grid.arrange(p1,p2,p3, ncol = 1)
dev.off()



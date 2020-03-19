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
reads_per_len.nmds = readRDS(file = "intermediate_RDS/20_gene_nmds_full.rds")

reads_per_len.nmds.metadata = left_join(
    data.frame(Sample = rownames(reads_per_len.nmds$points), MDS1 = reads_per_len.nmds$points[,1], MDS2 = reads_per_len.nmds$points[,2]),
    mapped.metadata
)

p1 = ggplot(reads_per_len.nmds.metadata, aes(hours.cumulative.RNA, MDS1, shape = Plot, color = hours.cumulative.RNA)) +
geom_point(size = 3) +
my_gg_theme +
labs(x = "Time (hrs)", color = "Time (hrs)") +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
theme(legend.title = element_text(size = 22))

p2 = ggplot(reads_per_len.nmds.metadata, aes(hours.cumulative.RNA, MDS2, shape = Plot, color = hours.cumulative.RNA)) +
geom_point(size = 3) +
my_gg_theme +
labs(x = "Time (hrs)", color = "Time (hrs)") +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
theme(legend.title = element_text(size = 22))

p3 = ggplot(reads_per_len.nmds.metadata, aes(hours.cumulative.RNA, MDS1, shape = Plot, color = Moisture)) +
geom_point(size = 3) +
my_gg_theme +
labs(x = "Time (hrs)", color = "Moisture") +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "sqrt") +
theme(legend.title = element_text(size = 22))

p4 = ggplot(reads_per_len.nmds.metadata, aes(hours.cumulative.RNA, MDS2, shape = Plot, color = Moisture)) +
geom_point(size = 3) +
my_gg_theme +
labs(x = "Time (hrs)", color = "Moisture") +
scale_color_gradient(low = "#0571b0", high = "#ca0020", trans = "sqrt") +
theme(legend.title = element_text(size = 22))

p5 = ggplot(reads_per_len.nmds.metadata, aes(hours.cumulative.RNA, MDS1, shape = Plot, color = Temperature)) +
geom_point(size = 3) +
my_gg_theme +
labs(x = "Time (hrs)", color = "Temperature") +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
theme(legend.title = element_text(size = 22))

p6 = ggplot(reads_per_len.nmds.metadata, aes(hours.cumulative.RNA, MDS2, shape = Plot, color = Temperature)) +
geom_point(size = 3) +
my_gg_theme +
labs(x = "Time (hrs)", color = "Temperature") +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
theme(legend.title = element_text(size = 22))

pdf("20_gene_NMDS_axes_over_time.pdf", width = 24, height = 10)
grid.arrange(p1,p3,p5,p2,p4,p6, ncol = 3)
dev.off()


reads_per_len.wide.P1.nmds = readRDS(file = "intermediate_RDS/20_gene_nmds.P1.rds")
reads_per_len.wide.P2.nmds = readRDS(file = "intermediate_RDS/20_gene_nmds.P2.rds")
reads_per_len.wide.P3.nmds = readRDS(file = "intermediate_RDS/20_gene_nmds.P3.rds")

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

p1 = ggplot(reads_per_len.nmds.P1.P2.P3.metadata, aes(hours.cumulative.RNA, MDS1)) +
geom_point(size = 2) +
facet_wrap(~Plot, ncol = 3, scales = "free") +
my_gg_theme +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
labs(x = "Time (hrs)") +
theme(legend.title = element_text(size = 22))

p2 = ggplot(reads_per_len.nmds.P1.P2.P3.metadata, aes(hours.cumulative.RNA, MDS2)) +
geom_point(size = 2) +
facet_wrap(~Plot, ncol = 3, scales = "free") +
my_gg_theme +
scale_color_gradient(low = "#0571b0", high = "#ca0020") +
labs(x = "Time (hrs)") +
theme(legend.title = element_text(size = 22))

pdf("20_gene_NMDS_axes_over_time_by_plot.pdf", width = 12, height = 6)
grid.arrange(p1,p2, ncol = 1)
dev.off()



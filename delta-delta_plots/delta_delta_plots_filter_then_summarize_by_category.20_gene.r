require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)
require(reshape2)
require(rlang)

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

#arguments for summarizing and output
args = commandArgs(T)
catToFilter = args[1]
filterValue = args[2]
catToSum = args[3]
outputFile = args[4]

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

#20 gene tax
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

#filter tables

read_count.filtered = read_count %>% filter(!!parse_expr(catToFilter) == filterValue)
ref_len.filtered = ref_len %>% filter(!!parse_expr(catToFilter) == filterValue)

print("filtered count tbl:")
read_count.filtered

#use min read count for rarefaction depth
minDepth = (colSums(read_count.filtered[,5:111]) %>% sort)[1]
print(c("minimum sequence depth: ", minDepth))
#Rarefy then normalize by target len and summarize by KO for KO delta-delta analysis

read_count.rarefied = data.frame(read_count.filtered[,1:4], t(rrarefy(t(read_count.filtered[,5:111]), sample = minDepth)) )

#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
read_count.rarefied.sumCat = read_count.rarefied %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
ref_len.sumCat = ref_len.filtered %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

print("summarized count tbl:")
print(read_count.rarefied.sumCat)
#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len = data.frame(read_count.rarefied.sumCat[,1], read_count.rarefied.sumCat[,2:108]/ref_len.sumCat[,2:108] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0




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

#distance calculations

reads_per_len.P1.distBC = vegdist(t(sqrt(reads_per_len.wide.P1[
    #The indexing term filters outs that only appear in one sample
    rowSums(reads_per_len.wide.P1[,2:length(colnames(reads_per_len.wide.P1))] > 0) >= 2,
    2:length(colnames(reads_per_len.wide.P1))])),
    "bray")

reads_per_len.P2.distBC = vegdist(t(sqrt(reads_per_len.wide.P2[
#The indexing term filters outs that only appear in one sample
rowSums(reads_per_len.wide.P2[,2:length(colnames(reads_per_len.wide.P2))] > 0) >= 2,
2:length(colnames(reads_per_len.wide.P2))])),
"bray")

reads_per_len.P3.distBC = vegdist(t(sqrt(reads_per_len.wide.P3[
#The indexing term filters outs that only appear in one sample
rowSums(reads_per_len.wide.P3[,2:length(colnames(reads_per_len.wide.P3))] > 0) >= 2,
2:length(colnames(reads_per_len.wide.P3))])),
"bray")

#melt trinagular matrix. melt will add zeros for blank values upper triangle so need to remove
reads_per_len.P1.distBC.long = subset(melt(reads_per_len.P1.distBC %>% as.matrix), value!=0)
reads_per_len.P2.distBC.long = subset(melt(reads_per_len.P2.distBC %>% as.matrix), value!=0)
reads_per_len.P3.distBC.long = subset(melt(reads_per_len.P3.distBC %>% as.matrix), value!=0)

###################################################################
#Calculate Euclidean distances for environmental variables and time

P1.time = data.frame(Time = metadata %>% filter(Plot == "P1") %>% select(hours.cumulative.RNA))
rownames(P1.time) = filter(metadata, Plot == "P1")$Sample
P2.time = data.frame(Time = metadata %>% filter(Plot == "P2") %>% select(hours.cumulative.RNA))
rownames(P2.time) = filter(metadata, Plot == "P2")$Sample
P3.time = data.frame(Time = metadata %>% filter(Plot == "P3") %>% select(hours.cumulative.RNA))
rownames(P3.time) = filter(metadata, Plot == "P3")$Sample

P1.time.euc = vegdist(P1.time, "euclidean")
P2.time.euc = vegdist(P2.time, "euclidean")
P3.time.euc = vegdist(P3.time, "euclidean")

P1.timeEuc.long = subset(melt(P1.time.euc %>% as.matrix), value!=0)
P2.timeEuc.long = subset(melt(P2.time.euc %>% as.matrix), value!=0)
P3.timeEuc.long = subset(melt(P3.time.euc %>% as.matrix), value!=0)

#Calculate temp distance
P1.temp = data.frame(temp = metadata %>% filter(Plot == "P1") %>% select(Temperature))
rownames(P1.temp) = filter(metadata, Plot == "P1")$Sample
P2.temp = data.frame(temp = metadata %>% filter(Plot == "P2") %>% select(Temperature))
rownames(P2.temp) = filter(metadata, Plot == "P2")$Sample
P3.temp = data.frame(temp = metadata %>% filter(Plot == "P3") %>% select(Temperature))
rownames(P3.temp) = filter(metadata, Plot == "P3")$Sample

P1.temp.euc = vegdist(P1.temp, "euclidean")
P2.temp.euc = vegdist(P2.temp, "euclidean")
P3.temp.euc = vegdist(P3.temp, "euclidean")

P1.tempEuc.long = subset(melt(P1.temp.euc %>% as.matrix), value!=0)
P2.tempEuc.long = subset(melt(P2.temp.euc %>% as.matrix), value!=0)
P3.tempEuc.long = subset(melt(P3.temp.euc %>% as.matrix), value!=0)

#Calculate moisture distance
P1.moist = data.frame(moist = metadata %>% filter(Plot == "P1") %>% select(Moisture))
rownames(P1.moist) = filter(metadata, Plot == "P1")$Sample
P2.moist = data.frame(moist = metadata %>% filter(Plot == "P2") %>% select(Moisture))
rownames(P2.moist) = filter(metadata, Plot == "P2")$Sample
P3.moist = data.frame(moist = metadata %>% filter(Plot == "P3") %>% select(Moisture))
rownames(P3.moist) = filter(metadata, Plot == "P3")$Sample

P1.moist.euc = vegdist(P1.moist, "euclidean")
P2.moist.euc = vegdist(P2.moist, "euclidean")
P3.moist.euc = vegdist(P3.moist, "euclidean")

P1.moistEuc.long = subset(melt(P1.moist.euc %>% as.matrix), value!=0)
P2.moistEuc.long = subset(melt(P2.moist.euc %>% as.matrix), value!=0)
P3.moistEuc.long = subset(melt(P3.moist.euc %>% as.matrix), value!=0)



#make dfs
#time
P1_time_dist.long = full_join(reads_per_len.P1.distBC.long, P1.timeEuc.long, by = c("Var1", "Var2"))
P2_time_dist.long = full_join(reads_per_len.P2.distBC.long, P2.timeEuc.long, by = c("Var1", "Var2"))
P3_time_dist.long = full_join(reads_per_len.P3.distBC.long, P3.timeEuc.long, by = c("Var1", "Var2"))

P1_time_dist.long$Plot = "P1"
P2_time_dist.long$Plot = "P2"
P3_time_dist.long$Plot = "P3"

P1P2P3_time_dist.long = rbind(P1_time_dist.long, P2_time_dist.long, P3_time_dist.long)
colnames(P1P2P3_time_dist.long) = c("Var1", "Var2", "bray.curtis", "delta.time", "Plot")

P1P2P3_time_dist.long$time.breaks = cut(P1P2P3_time_dist.long$delta.time, breaks = seq(1,73,2), labels = seq(2, 72, 2))

#temp
P1_temp_dist.long = full_join(reads_per_len.P1.distBC.long, P1.tempEuc.long, by = c("Var1", "Var2"))
P2_temp_dist.long = full_join(reads_per_len.P2.distBC.long, P2.tempEuc.long, by = c("Var1", "Var2"))
P3_temp_dist.long = full_join(reads_per_len.P3.distBC.long, P3.tempEuc.long, by = c("Var1", "Var2"))


P1_temp_dist.long$Plot = "P1"
P2_temp_dist.long$Plot = "P2"
P3_temp_dist.long$Plot = "P3"

P1P2P3_temp_dist.long = rbind(P1_temp_dist.long, P2_temp_dist.long, P3_temp_dist.long)
colnames(P1P2P3_temp_dist.long) = c("Var1", "Var2", "bray.curtis", "delta.temp", "Plot")

#max(P1P2P3_temp_dist.long$delta.temp)
P1P2P3_temp_dist.long$temp.breaks = cut(P1P2P3_temp_dist.long$delta.temp, breaks = seq(0,40,2), labels = seq(2, 40, 2))

#moisture
P1_moist_dist.long = full_join(reads_per_len.P1.distBC.long, P1.moistEuc.long, by = c("Var1", "Var2"))
P2_moist_dist.long = full_join(reads_per_len.P2.distBC.long, P2.moistEuc.long, by = c("Var1", "Var2"))
P3_moist_dist.long = full_join(reads_per_len.P3.distBC.long, P3.moistEuc.long, by = c("Var1", "Var2"))

P1_moist_dist.long$Plot = "P1"
P2_moist_dist.long$Plot = "P2"
P3_moist_dist.long$Plot = "P3"

P1P2P3_moist_dist.long = rbind(P1_moist_dist.long, P2_moist_dist.long, P3_moist_dist.long)
colnames(P1P2P3_moist_dist.long) = c("Var1", "Var2", "bray.curtis", "delta.moist", "Plot")

#range(as.numeric(P1P2P3_moist_dist.long$delta.moist))
P1P2P3_moist_dist.long$moist.breaks = cut(P1P2P3_moist_dist.long$delta.moist, breaks = seq(0,1.8,.1), labels = seq(10, 180, 10))



#PLOTS
#TIME
pdf(file = outputFile, width = 6, height = 6)
#pdf("delta_time_bray.KO.pdf", width = 8, height = 5)
print(
ggplot(P1P2P3_time_dist.long, aes(delta.time, bray.curtis, color = Plot)) +
facet_wrap(~Plot, ncol = 1, scales = "free_y") +
geom_point(size= 0.2) +
scale_color_manual(values = c("P1" = "#E69F00", "P2" = "#56B4E9", "P3" = "#009E73"), guide = F) +
labs(x = "Delta time (hrs)", y = "Bray-Curtis distance") +
scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
#scale_y_continuous(limits = c(0.1, 0.4)) +
#geom_smooth(method = "loess") +
my_gg_theme
)
#dev.off()

#pdf("delta_time_bray_boxplot.KO.pdf", width = 12, height = 4)
print(
ggplot(P1P2P3_time_dist.long, aes(time.breaks, bray.curtis, color = Plot)) +
facet_wrap(~Plot, ncol = 1) +
geom_boxplot(outlier.size = 0.5) +
scale_color_manual(values = c("P1" = "#E69F00", "P2" = "#56B4E9", "P3" = "#009E73"), guide = F) +
labs(x = "Delta time (hrs)", y = "Bray-Curtis distance") +
scale_x_discrete(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
#scale_y_continuous(limits = c(0.1, 0.4)) +
#geom_smooth() +
my_gg_theme
)
#dev.off()

#TEMPERATURE
#pdf("delta_temp_bray.KO.pdf", width = 8, height = 5)
print(
ggplot(P1P2P3_temp_dist.long, aes(delta.temp, bray.curtis, color = Plot)) +
facet_wrap(~Plot, ncol = 1) +
geom_point(size= 0.2) +
scale_color_manual(values = c("P1" = "#E69F00", "P2" = "#56B4E9", "P3" = "#009E73"), guide = F) +
labs(x = expression("Delta temperature ("*~degree*C*")"), y = "Bray-Curtis distance") +
#scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
#scale_y_continuous(limits = c(0.1, 0.4)) +
#geom_smooth(method = "loess") +
my_gg_theme
)
#dev.off()

#pdf("delta_temp_bray_boxplot.KO.pdf", width = 12, height = 4)
print(
ggplot(filter(P1P2P3_temp_dist.long, !is.na(temp.breaks)), aes(temp.breaks, bray.curtis, color = Plot)) +
facet_wrap(~Plot, ncol = 1) +
geom_boxplot(outlier.size = 0.5, notch = F) +
scale_color_manual(values = c("P1" = "#E69F00", "P2" = "#56B4E9", "P3" = "#009E73"), guide = F) +
labs(x = expression("Delta temperature ("*~degree*C*")"), y = "Bray-Curtis distance") +
scale_x_discrete(breaks = c(1, 5, 10, 15, 20, 25, 30)) +
#scale_y_continuous(limits = c(0.1, 0.4)) +
#geom_smooth() +
my_gg_theme
)
#dev.off()

#moisture
#pdf("delta_moist_bray.KO.pdf", width = 8, height = 5)
print(
ggplot(P1P2P3_moist_dist.long, aes(delta.moist, bray.curtis, color = Plot)) +
facet_wrap(~Plot, ncol = 1) +
geom_point(size= 0.2) +
scale_color_manual(values = c("P1" = "#E69F00", "P2" = "#56B4E9", "P3" = "#009E73"), guide = F) +
labs(x = "Delta moisture", y = "Bray-Curtis distance") +
#scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
#scale_y_continuous(limits = c(0.1, 0.4)) +
geom_smooth(method = "loess") +
my_gg_theme
)
#dev.off()

#pdf("delta_moist_bray_boxplot.KO.pdf", width = 12, height = 4)
print(
ggplot(filter(P1P2P3_moist_dist.long, !is.na(moist.breaks)), aes(moist.breaks, bray.curtis, color = Plot)) +
facet_wrap(~Plot, ncol = 1) +
#geom_jitter(size = 0.25, color = "black") +
geom_boxplot(outlier.size = 0.5) +
scale_color_manual(values = c("P1" = "#E69F00", "P2" = "#56B4E9", "P3" = "#009E73"), guide = F) +
labs(x = expression("Delta moisture (%)"), y = "Bray-Curtis distance") +
scale_x_discrete(breaks = c(20, 40, 60, 80, 100, 120, 140, 160, 180)) +
#scale_x_discrete(breaks = c(10,50,100,150)) +
#scale_y_continuous(limits = c(0.1, 0.4)) +
#geom_smooth() +
my_gg_theme
)

dev.off()

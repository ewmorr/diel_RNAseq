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


#make dfs
#time
P1_time_dist.long = full_join(P1.tempEuc.long, P1.timeEuc.long, by = c("Var1", "Var2"))
P2_time_dist.long = full_join(P2.tempEuc.long, P2.timeEuc.long, by = c("Var1", "Var2"))
P3_time_dist.long = full_join(P3.tempEuc.long, P3.timeEuc.long, by = c("Var1", "Var2"))

P1_time_dist.long$Plot = "P1"
P2_time_dist.long$Plot = "P2"
P3_time_dist.long$Plot = "P3"

P1P2P3_time_dist.long = rbind(P1_time_dist.long, P2_time_dist.long, P3_time_dist.long)
colnames(P1P2P3_time_dist.long) = c("Var1", "Var2", "delta.temp", "delta.time", "Plot")

P1P2P3_time_dist.long$time.breaks = cut(P1P2P3_time_dist.long$delta.time, breaks = seq(1,73,2), labels = seq(2, 72, 2))

#PLOTS
#TIME
pdf(file = "delta-delta_temp_time.pdf", width = 6, height = 5)
#pdf("delta_time_bray.KO.pdf", width = 8, height = 5)
print(
ggplot(P1P2P3_time_dist.long, aes(delta.time, delta.temp, color = Plot)) +
facet_wrap(~Plot, ncol = 1, scales = "free_y") +
geom_point(size= 0.2) +
scale_color_manual(values = c("P1" = "#E69F00", "P2" = "#56B4E9", "P3" = "#009E73"), guide = F) +
labs(x = "Delta time (hrs)", y = "Temp. Euclidean distance") +
scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
#scale_y_continuous(limits = c(0.1, 0.4)) +
#geom_smooth(method = "loess") +
my_gg_theme
)


#pdf("delta_time_bray_boxplot.KO.pdf", width = 12, height = 4)
print(
ggplot(P1P2P3_time_dist.long, aes(time.breaks, delta.temp, color = Plot)) +
facet_wrap(~Plot, ncol = 1) +
geom_boxplot(outlier.size = 0.5) +
scale_color_manual(values = c("P1" = "#E69F00", "P2" = "#56B4E9", "P3" = "#009E73"), guide = F) +
labs(x = "Delta time (hrs)", y = "Temp. Euclidean distance") +
scale_x_discrete(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
#scale_y_continuous(limits = c(0.1, 0.4)) +
#geom_smooth() +
my_gg_theme
)
dev.off()

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

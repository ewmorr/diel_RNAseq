require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)
require(lubridate)

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

total_mapped = read.csv("totals_mapped_input.reformat.txt", sep = "\t", header = F)
colnames(total_mapped) = c("Sample", "totalReads", "readsMappedCov", "totalBases", "targetsLength", "basesMapped", "avgCov")
total_mapped = total_mapped %>% as.tbl
total_mapped = total_mapped %>% filter(Sample != "P2T31_metaT")

#metadata table
metadata = read.table("metadata.txt", header = T, sep = "\t")
metadata = metadata %>% as.tbl
metadata = metadata %>% filter(Sample != "P2T31_metaT")

metadata$timeOfDay.RNA = (as.numeric(as.POSIXlt(metadata$timeOfDay.RNA, format = "%H:%M") %>% hour) + as.numeric((as.POSIXlt(metadata$timeOfDay.RNA, format = "%H:%M") %>% minute)/60) ) %>% signif(2)

mapped.metadata = full_join(total_mapped, metadata, by = "Sample")


#RNA seq data, ref_len is the total length of reference contigs per category
ref_len = read.table("read_counts.KO.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/KO_rarefied_count.rds")

#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
read_count.rarefied.sumCat = read_count.rarefied %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
ref_len.sumCat = ref_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

print("summarized count tbl:")
print(read_count.rarefied.sumCat)
#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len = data.frame(read_count.rarefied.sumCat[,1], read_count.rarefied.sumCat[,2:108]/ref_len.sumCat[,2:108] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

mapped.metadata = mapped.metadata[match(colnames(reads_per_len[2:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

###################################
#Read KO names

KO_names = read.table("KO_term_names_rarefied_table.txt", sep = "\t")
colnames(KO_names) = c("sumCat", "name")

#look at totals by sumCat
reads_per_len.sumCat = reads_per_len %>% group_by(sumCat) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
reads_per_len.sumCat.rowSum = data.frame(sumCat = reads_per_len.sumCat$sumCat, sum = reads_per_len.sumCat[2:length(colnames(reads_per_len.sumCat))] %>% rowSums)
reads_per_len.sumCat.rowSum[order(reads_per_len.sumCat.rowSum$sum),]

top_ten_KO = (reads_per_len.sumCat.rowSum[order(reads_per_len.sumCat.rowSum$sum, decreasing = T),])[1:12,]

write.table(left_join(top_ten_KO, KO_names), file = "top_12_KO.txt", sep = "\t", quote = F, row.names = F, col.names = F)

write.table(left_join(filter(reads_per_len.sumCat.rowSum[order(reads_per_len.sumCat.rowSum$sum, decreasing = T),], sum > 0), KO_names), file = "all_KO_with_names.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#Relative abundance
reads_per_len.sumCat.relAbd = data.frame(
reads_per_len.sumCat[,1],
t(t(reads_per_len.sumCat[,2:length(reads_per_len.sumCat)])/colSums(reads_per_len.sumCat[,2:length(reads_per_len.sumCat)]))
)

saveRDS(reads_per_len.sumCat.relAbd, "KO_relative_abd.metaG.rds")
write.table(reads_per_len.sumCat.relAbd, file = "KO_rel_abd.metaG.txt", sep = "\t", row.names = F, quote = F)

reads_per_len.sumCat.relAbd.long = pivot_longer(reads_per_len.sumCat.relAbd,
cols = -sumCat,
names_to = "Sample",
values_to = "rel_abd"
)

reads_per_len.sumCat.relAbd.long.metadata = left_join(reads_per_len.sumCat.relAbd.long, mapped.metadata, by = "Sample")
require(RColorBrewer)
reads_per_len.sumCat.relAbd.long.metadata.topTen = reads_per_len.sumCat.relAbd.long.metadata %>% filter(sumCat %in% top_ten_KO$sumCat)

reads_per_len.sumCat.relAbd.long.metadata.topTen$sumCat = factor(reads_per_len.sumCat.relAbd.long.metadata.topTen$sumCat, levels = top_ten_KO$sumCat)

p1 = ggplot(reads_per_len.sumCat.relAbd.long.metadata.topTen, aes(TimePoint, rel_abd, fill = sumCat)) +
geom_col() +
facet_wrap(~Plot) +
scale_fill_brewer(palette = "Paired") +
labs(x = "Time point", y = "Rel. abd.") +
my_gg_theme

pdf("KO_rel_abd_top_ten.metaG.pdf", width = 14, height = 6)
p1
dev.off()



p1 = ggplot(reads_per_len.sumCat.relAbd.long.metadata.topTen, aes(TimePoint, rel_abd)) +
geom_point() +
#geom_smooth() +
facet_grid(sumCat~Plot, scale = "free_y") +
scale_fill_brewer(palette = "Paired") +
labs(x = "Time point", y = "Rel. abd.") +
my_gg_theme

p2 = ggplot(reads_per_len.sumCat.relAbd.long.metadata.topTen %>% filter(TimePoint > 5), aes(timeOfDay.RNA, rel_abd)) +
geom_point() +
geom_smooth() +
facet_grid(Plot~sumCat) +
scale_fill_brewer(palette = "Paired") +
labs(x = "Time of Day", y = "Rel. abd.") +
my_gg_theme

pdf("KO_rel_abd_top_ten_scatterPlots.metaG.pdf", width = 16, height = 20)
print(p1)
print(p2)
dev.off()

pdf("KO_rel_abd_top_ten_TOD_scatterPlots.metaG.pdf", width = 28, height = 10)
#print(p1)
print(p2)
dev.off()

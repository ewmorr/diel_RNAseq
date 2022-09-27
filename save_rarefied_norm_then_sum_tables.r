require(tidyverse)
require(lubridate)
require(rlang)

#METADATA
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

metadata$timeOfDay.RNA = (as.numeric(as.POSIXlt(metadata$timeOfDay.RNA, format = "%H:%M") %>% hour) + as.numeric((as.POSIXlt(metadata$timeOfDay.RNA, format = "%H:%M") %>% minute)/60) ) %>% signif(2)

mapped.metadata = full_join(total_mapped, metadata, by = "Sample")

write.table(mapped.metadata, "tables_for_adam/metadata.txt", sep = "\t", quote = F, row.names = F)

#KO terms
catToSum = "KO"

ref_len = read.table("read_counts.KO.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/KO_rarefied_count.rds")

reads_per_len = data.frame(read_count.rarefied[,1:8], read_count.rarefied[,9:115]/ref_len[,9:115] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

reads_per_len.sumCat = reads_per_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

print("summarized count tbl:")
print(reads_per_len.sumCat)

reads_per_len = reads_per_len.sumCat
colnames(reads_per_len)[1] = "KO"

write.table(reads_per_len, "tables_for_adam/KEGG.rarefied_normalized.txt", sep = "\t", quote = F, row.names = F)

#
#20GENE
#

catToSum = "Edge_num"

ref_len = read.table("read_counts.20_gene_taxonomy.metaG.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL
ref_len$Edge_num = as.factor(ref_len$Edge_num)

read_count.rarefied = readRDS(file = "intermediate_RDS/20_gene_rarefied_counts.metaG.rds")
read_count.rarefied$Edge_num = as.factor(read_count.rarefied$Edge_num)

reads_per_len = data.frame(read_count.rarefied[,1:4], read_count.rarefied[,5:111]/ref_len[,5:111] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
reads_per_len.sumCat = reads_per_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

reads_per_len %>% group_by(Edge_num) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

print("summarized count tbl:")
print(reads_per_len.sumCat)
#Divide read count by refLen to normalize reads recruited by contig length

reads_per_len = reads_per_len.sumCat
colnames(reads_per_len)[1] = "Edge_num"

reads_per_len.tax = right_join(
    data.frame(read_count.rarefied[,1:2]) %>% group_by(Edge_num) %>% unique(),
    reads_per_len
)


#investigate groups
reads_per_len.tax %>% group_by(Tax_name) %>% mutate(n = n()) %>% filter(n > 1)
reads_per_len.tax %>% group_by(Edge_num) %>% mutate(n = n()) %>% filter(n > 1)

write.table(reads_per_len.tax, "tables_for_adam/20_gene_tax.rarefied_normalized.txt", sep = "\t", quote = F, row.names = F)


#
#GH terms
#

catToSum = "Category"

ref_len = read.table("read_counts.GH.readCountRefLen.join", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/GH_rarefied_count.rds")

reads_per_len = data.frame(Category = read_count.rarefied[,1], read_count.rarefied[,2:108]/ref_len[,2:108] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
reads_per_len.sumCat = reads_per_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

print("summarized count tbl:")
print(reads_per_len.sumCat)
#Divide read count by refLen to normalize reads recruited by contig length

reads_per_len = reads_per_len.sumCat

colnames(reads_per_len)[1] = "CAZyme"

write.table(reads_per_len, "tables_for_adam/CAZy.rarefied_normalized.txt", sep = "\t", quote = F, row.names = F)

require(tidyverse)

read_count.rarefied = readRDS(file = "intermediate_RDS/20_gene_rarefied_counts.metaG.rds")

ref_len.20_gene_taxonomy = read.table("read_counts.20_gene_taxonomy.metaG.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len.20_gene_taxonomy$P2T31_metaT = NULL
ref_len.20_gene_taxonomy$X = NULL

#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len.20_gene_taxonomy = data.frame(read_count.rarefied[,1:4], read_count.rarefied[,5:111]/ref_len.20_gene_taxonomy[,5:111] ) %>% as.tbl
reads_per_len.20_gene_taxonomy[is.na(reads_per_len.20_gene_taxonomy)] = 0

#get totals by Tax_name
reads_per_len.Tax_name = reads_per_len.20_gene_taxonomy %>% group_by(Tax_name) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

saveRDS(reads_per_len.Tax_name, file = "intermediate_RDS/reads_per_len.20_gene_tax.rds")
write.table(reads_per_len.Tax_name, file = "reads_per_len.20_gene_tax.txt", sep = "\t", col.names = T, row.names = F, quote = F)

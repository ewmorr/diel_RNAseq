require(vegan)
require(tidyverse)


read_count = read.table("read_counts.KO.phylodist.readCount.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
read_count$P2T31_metaT = NULL
read_count$X = NULL

minDepth = (colSums(read_count[,9:115]) %>% sort)[1]
print(c("minimum sequence depth: ", minDepth))

read_count.rarefied = data.frame(read_count[,1:8], t(rrarefy(t(read_count[,9:115]), sample = minDepth)) )
saveRDS(read_count.rarefied, file = "intermediate_RDS/KO_rarefied_count.rds")

read_count = read.table("read_counts.20_gene_taxonomy.readCount.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
read_count$P2T31_metaT = NULL
read_count$X = NULL
read_count$Edge_num = as.factor(read_count$Edge_num)

minDepth = (colSums(read_count[,5:111]) %>% sort)[1]
print(c("minimum sequence depth: ", minDepth))
#Rarefy then normalize by target len and summarize by KO for KO delta-delta analysis

read_count.rarefied = data.frame(read_count[,1:4], t(rrarefy(t(read_count[,5:111]), sample = minDepth)) )
saveRDS(read_count.rarefied, file = "intermediate_RDS/20_gene_rarefied_count.rds")


read_count = read.table("read_counts.GH.readCount.join", header = T, sep = "\t") %>% as.tbl
read_count$P2T31_metaT = NULL
read_count$X = NULL

minDepth = (colSums(read_count[,2:108]) %>% sort)[1]
print(c("minimum sequence depth: ", minDepth))

read_count.rarefied = data.frame(read_count[,1], t(rrarefy(t(read_count[,2:108]), sample = minDepth)) )
saveRDS(read_count.rarefied, file = "intermediate_RDS/GH_rarefied_count.rds")

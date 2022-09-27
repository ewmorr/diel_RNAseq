require(tidyverse)
require(RobPer)
require(rlang)
require(lubridate)

source("~/ggplot_theme.txt")

####################################
#READ AND PROCESS TABLES
#########################

#arguments for summarizing and output
catToSum = "KO"
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

#RNA seq data, read_counts is the number of reads mapped per category (e.g. taxonomic string)


#RNA seq data, ref_len is the total length of reference contigs per category
ref_len = read.table("read_counts.KO.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/KO_rarefied_count.rds")

reads_per_len = data.frame(read_count.rarefied[,1:8], read_count.rarefied[,9:115]/ref_len[,9:115] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0


#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
reads_per_len.sumCat = reads_per_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

KO_cat_sum = data.frame(
"name" = reads_per_len.sumCat$sumCat,
"RPK.sum" = rowSums(reads_per_len.sumCat[,2:ncol(reads_per_len.sumCat)])*1000
)
KO_cat_sum$name = as.character(KO_cat_sum$name)

reads_per_len = as.data.frame(reads_per_len.sumCat)


KO_names = read.table("KO_term_names_rarefied_table.txt", sep = "\t")
colnames(KO_names) = c("cat", "name")

#extract plots
P1.samples = mapped.metadata %>%
    filter(Plot == "P1") %>%
    select(Sample)
P2.samples = mapped.metadata %>%
    filter(Plot == "P2") %>%
    select(Sample)
P3.samples = mapped.metadata %>%
    filter(Plot == "P3") %>%
    select(Sample)

reads_per_len.P1 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P1.samples$Sample)]
reads_per_len.P2 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P2.samples$Sample)]
reads_per_len.P3 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P3.samples$Sample)]

#filter by frequency
reads_per_len = reads_per_len[rowSums(reads_per_len[,2:ncol(reads_per_len)] > 0) > 107 *.33,]
mapped.metadata = mapped.metadata[match(colnames(reads_per_len[2:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

reads_per_len.P1 = reads_per_len.P1[rowSums(reads_per_len.P1[,2:ncol(reads_per_len.P1)] > 0) > ncol(reads_per_len.P1) *.33,]
reads_per_len.P2 = reads_per_len.P2[rowSums(reads_per_len.P2[,2:ncol(reads_per_len.P2)] > 0) > ncol(reads_per_len.P2) *.33,]
reads_per_len.P3 = reads_per_len.P3[rowSums(reads_per_len.P3[,2:ncol(reads_per_len.P3)] > 0) > ncol(reads_per_len.P3) *.33,]


mapped.metadata.P1 = mapped.metadata %>% filter(Plot == "P1")
mapped.metadata.P2 = mapped.metadata %>% filter(Plot == "P2")
mapped.metadata.P3 = mapped.metadata %>% filter(Plot == "P3")

###############################
#RobPer by row
##############

robper.overall = data.frame(
    cat = vector(mode = "character", length = nrow(reads_per_len)),
    R2 = vector(mode = "numeric", length = nrow(reads_per_len)),
    stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len)){
    temp.dat = data.frame(
        KO = t(reads_per_len[i,2:ncol(reads_per_len)])*1000,
        hrs = mapped.metadata$hours.cumulative.RNA
    )
    
    temp.RobPer = RobPer(
        ts = temp.dat,
        weighting = F,
        periods = 24,
        regression = "L2",
        model= "sine"
    )

    robper.overall$cat[i] = as.character(reads_per_len[i,1])
    robper.overall$R2[i] = temp.RobPer
}

robper.overall.ordered = robper.overall[order(-robper.overall$R2),]
robper.overall.ordered = left_join(robper.overall.ordered, KO_names)

#P1
robper.P1 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P1)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P1)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P1)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P1[i,2:ncol(reads_per_len.P1)])*1000,
    hrs = mapped.metadata.P1$hours.cumulative.RNA
    )
    
    temp.RobPer = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.P1$cat[i] = as.character(reads_per_len.P1[i,1])
    robper.P1$R2[i] = temp.RobPer
}

robper.P1.ordered = robper.P1[order(-robper.P1$R2),]
robper.P1.ordered = left_join(robper.P1.ordered, KO_names)

plot(robper.P1.ordered$R2)
#P2
robper.P2 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P2)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P2)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P2)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P2[i,2:ncol(reads_per_len.P2)])*1000,
    hrs = mapped.metadata.P2$hours.cumulative.RNA
    )
    
    temp.RobPer = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.P2$cat[i] = as.character(reads_per_len.P2[i,1])
    robper.P2$R2[i] = temp.RobPer
}

robper.P2.ordered = robper.P2[order(-robper.P2$R2),]
robper.P2.ordered = left_join(robper.P2.ordered, KO_names)

plot(robper.P2.ordered$R2)

#P3
robper.P3 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P3)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P3)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P3)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P3[i,2:ncol(reads_per_len.P3)])*1000,
    hrs = mapped.metadata.P3$hours.cumulative.RNA
    )
    
    temp.RobPer = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.P3$cat[i] = as.character(reads_per_len.P3[i,1])
    robper.P3$R2[i] = temp.RobPer
}

robper.P3.ordered = robper.P3[order(-robper.P3$R2),]
robper.P3.ordered = left_join(robper.P3.ordered, KO_names)

plot(robper.P1.ordered$R2)

####End KO

##############################
#TAXONOMIC
catToSum = "Edge_num"


#RNA seq data, ref_len is the total length of reference contigs per category
ref_len = read.table("read_counts.20_gene_taxonomy.metaG.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL
ref_len$Edge_num = as.character(ref_len$Edge_num)

read_count.rarefied = readRDS(file = "intermediate_RDS/20_gene_rarefied_counts.metaG.rds")
read_count.rarefied$Edge_num = as.character(read_count.rarefied$Edge_num)

reads_per_len = data.frame(read_count.rarefied[,1:4], read_count.rarefied[,5:111]/ref_len[,5:111] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

taxa_names = read_count.rarefied[1:2]
taxa_names = aggregate(Tax_name ~ Edge_num, taxa_names, unique)
colnames(taxa_names) = c("cat", "name")
#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
reads_per_len.sumCat = reads_per_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

taxa_cat_sum = data.frame(
"name" = reads_per_len.sumCat$sumCat,
"RPK.sum" = rowSums(reads_per_len.sumCat[,2:ncol(reads_per_len.sumCat)])*1000
)
taxa_cat_sum$name = as.character(taxa_cat_sum$name)

reads_per_len = reads_per_len.sumCat

#extract plots
P1.samples = mapped.metadata %>%
filter(Plot == "P1") %>%
select(Sample)
P2.samples = mapped.metadata %>%
filter(Plot == "P2") %>%
select(Sample)
P3.samples = mapped.metadata %>%
filter(Plot == "P3") %>%
select(Sample)

reads_per_len.P1 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P1.samples$Sample)]
reads_per_len.P2 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P2.samples$Sample)]
reads_per_len.P3 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P3.samples$Sample)]

#filter by frequency
reads_per_len = reads_per_len[rowSums(reads_per_len[,2:ncol(reads_per_len)] > 0) > 107 *.33,]
mapped.metadata = mapped.metadata[match(colnames(reads_per_len[2:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

reads_per_len.P1 = reads_per_len.P1[rowSums(reads_per_len.P1[,2:ncol(reads_per_len.P1)] > 0) > ncol(reads_per_len.P1) *.33,]
reads_per_len.P2 = reads_per_len.P2[rowSums(reads_per_len.P2[,2:ncol(reads_per_len.P2)] > 0) > ncol(reads_per_len.P2) *.33,]
reads_per_len.P3 = reads_per_len.P3[rowSums(reads_per_len.P3[,2:ncol(reads_per_len.P3)] > 0) > ncol(reads_per_len.P3) *.33,]

mapped.metadata.P1 = mapped.metadata %>% filter(Plot == "P1")
mapped.metadata.P2 = mapped.metadata %>% filter(Plot == "P2")
mapped.metadata.P3 = mapped.metadata %>% filter(Plot == "P3")


###############################
#robper.20gene by row
##############

robper.20gene.overall = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len)){
    temp.dat = data.frame(
    KO = t(reads_per_len[i,2:ncol(reads_per_len)])*1000,
    hrs = mapped.metadata$hours.cumulative.RNA
    )
    
    temp.robper.20gene = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.20gene.overall$cat[i] = reads_per_len[i,1]
    robper.20gene.overall$R2[i] = temp.robper.20gene
}

robper.20gene.overall$cat = as.character(robper.20gene.overall$cat)

robper.20gene.overall.ordered = robper.20gene.overall[order(-robper.20gene.overall$R2),]
robper.20gene.overall.ordered = left_join(robper.20gene.overall.ordered, taxa_names)

#P1
robper.20gene.P1 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P1)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P1)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P1)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P1[i,2:ncol(reads_per_len.P1)])*1000,
    hrs = mapped.metadata.P1$hours.cumulative.RNA
    )
    
    temp.robper.20gene = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.20gene.P1$cat[i] = reads_per_len.P1[i,1]
    robper.20gene.P1$R2[i] = temp.robper.20gene
}
robper.20gene.P1$cat = as.character(robper.20gene.P1$cat)

robper.20gene.P1.ordered = robper.20gene.P1[order(-robper.20gene.P1$R2),]
robper.20gene.P1.ordered = left_join(robper.20gene.P1.ordered, taxa_names)

plot(robper.20gene.P1.ordered$R2)
#P2
robper.20gene.P2 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P2)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P2)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P2)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P2[i,2:ncol(reads_per_len.P2)])*1000,
    hrs = mapped.metadata.P2$hours.cumulative.RNA
    )
    
    temp.robper.20gene = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.20gene.P2$cat[i] = reads_per_len.P2[i,1]
    robper.20gene.P2$R2[i] = temp.robper.20gene
}
robper.20gene.P2$cat = as.character(robper.20gene.P2$cat)

robper.20gene.P2.ordered = robper.20gene.P2[order(-robper.20gene.P2$R2),]
robper.20gene.P2.ordered = left_join(robper.20gene.P2.ordered, taxa_names)

plot(robper.20gene.P2.ordered$R2)

#P3
robper.20gene.P3 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P3)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P3)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P3)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P3[i,2:ncol(reads_per_len.P3)])*1000,
    hrs = mapped.metadata.P3$hours.cumulative.RNA
    )
    
    temp.robper.20gene = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.20gene.P3$cat[i] = reads_per_len.P3[i,1]
    robper.20gene.P3$R2[i] = temp.robper.20gene
}
robper.20gene.P3$cat = as.character(robper.20gene.P3$cat)

robper.20gene.P3.ordered = robper.20gene.P3[order(-robper.20gene.P3$R2),]
robper.20gene.P3.ordered = left_join(robper.20gene.P3.ordered, taxa_names)

plot(robper.20gene.P3.ordered$R2)

####End taxonomic

##############################
#GH
catToSum = "Category"

ref_len = read.table("read_counts.GH.readCountRefLen.join", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/GH_rarefied_count.rds")


reads_per_len = data.frame(Category = read_count.rarefied[,1], read_count.rarefied[,2:108]/ref_len[,2:108] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0
reads_per_len$Category = reads_per_len$Category %>% as.character
#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
reads_per_len.sumCat = reads_per_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
reads_per_len = reads_per_len.sumCat

GH_cat_sum = data.frame(
    "name" = reads_per_len.sumCat$sumCat,
    "RPK.sum" = rowSums(reads_per_len.sumCat[,2:ncol(reads_per_len.sumCat)])*1000,
    stringsAsFactors = F
)
GH_cat_sum$name = as.character(GH_cat_sum$name)
#extract plots
P1.samples = mapped.metadata %>%
filter(Plot == "P1") %>%
select(Sample)
P2.samples = mapped.metadata %>%
filter(Plot == "P2") %>%
select(Sample)
P3.samples = mapped.metadata %>%
filter(Plot == "P3") %>%
select(Sample)

reads_per_len.P1 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P1.samples$Sample)]
reads_per_len.P2 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P2.samples$Sample)]
reads_per_len.P3 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", P3.samples$Sample)]

#filter by frequency
reads_per_len = reads_per_len[rowSums(reads_per_len[,2:ncol(reads_per_len)] > 0) > 107 *.33,]
mapped.metadata = mapped.metadata[match(colnames(reads_per_len[2:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

reads_per_len.P1 = reads_per_len.P1[rowSums(reads_per_len.P1[,2:ncol(reads_per_len.P1)] > 0) > ncol(reads_per_len.P1) *.33,]
reads_per_len.P2 = reads_per_len.P2[rowSums(reads_per_len.P2[,2:ncol(reads_per_len.P2)] > 0) > ncol(reads_per_len.P2) *.33,]
reads_per_len.P3 = reads_per_len.P3[rowSums(reads_per_len.P3[,2:ncol(reads_per_len.P3)] > 0) > ncol(reads_per_len.P3) *.33,]

mapped.metadata.P1 = mapped.metadata %>% filter(Plot == "P1")
mapped.metadata.P2 = mapped.metadata %>% filter(Plot == "P2")
mapped.metadata.P3 = mapped.metadata %>% filter(Plot == "P3")


###############################
#robper.GH by row
##############

robper.GH.overall = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len)){
    temp.dat = data.frame(
    KO = t(reads_per_len[i,2:ncol(reads_per_len)])*1000,
    hrs = mapped.metadata$hours.cumulative.RNA
    )
    
    temp.robper.GH = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.GH.overall$cat[i] = as.character(reads_per_len[i,1])
    robper.GH.overall$R2[i] = temp.robper.GH
}

robper.GH.overall.ordered = robper.GH.overall[order(-robper.GH.overall$R2),]
robper.GH.overall.ordered$name = robper.GH.overall.ordered$cat

#P1
robper.GH.P1 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P1)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P1)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P1)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P1[i,2:ncol(reads_per_len.P1)])*1000,
    hrs = mapped.metadata.P1$hours.cumulative.RNA
    )
    
    temp.robper.GH = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.GH.P1$cat[i] = reads_per_len.P1[i,1]
    robper.GH.P1$R2[i] = temp.robper.GH
}

robper.GH.P1.ordered = robper.GH.P1[order(-robper.GH.P1$R2),]
robper.GH.P1.ordered$name = robper.GH.P1.ordered$cat

plot(robper.GH.P1.ordered$R2)
#P2
robper.GH.P2 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P2)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P2)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P2)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P2[i,2:ncol(reads_per_len.P2)])*1000,
    hrs = mapped.metadata.P2$hours.cumulative.RNA
    )
    
    temp.robper.GH = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.GH.P2$cat[i] = reads_per_len.P2[i,1]
    robper.GH.P2$R2[i] = temp.robper.GH
}

robper.GH.P2.ordered = robper.GH.P2[order(-robper.GH.P2$R2),]
robper.GH.P2.ordered$name = robper.GH.P2.ordered$cat

plot(robper.GH.P2.ordered$R2)

#P3
robper.GH.P3 = data.frame(
cat = vector(mode = "character", length = nrow(reads_per_len.P3)),
R2 = vector(mode = "numeric", length = nrow(reads_per_len.P3)),
stringsAsFactors = F
)

for( i in 1:nrow(reads_per_len.P3)){
    temp.dat = data.frame(
    KO = t(reads_per_len.P3[i,2:ncol(reads_per_len.P3)])*1000,
    hrs = mapped.metadata.P3$hours.cumulative.RNA
    )
    
    temp.robper.GH = RobPer(
    ts = temp.dat,
    weighting = F,
    periods = 24,
    regression = "L2",
    model= "sine"
    )
    
    robper.GH.P3$cat[i] = reads_per_len.P3[i,1]
    robper.GH.P3$R2[i] = temp.robper.GH
}

robper.GH.P3.ordered = robper.GH.P3[order(-robper.GH.P3$R2),]
robper.GH.P3.ordered$name = robper.GH.P3.ordered$cat

plot(robper.GH.P3.ordered$R2)



##########################
#Combine df for plotting
########################
robper.overall.ordered$category = "KO"
robper.overall.ordered$plot = "Overall"
robper.overall.ordered$rank = seq(1,nrow(robper.overall.ordered), 1)
robper.P1.ordered$category = "KO"
robper.P1.ordered$plot = "P1"
robper.P1.ordered$rank = seq(1,nrow(robper.P1.ordered), 1)
robper.P2.ordered$category = "KO"
robper.P2.ordered$plot = "P2"
robper.P2.ordered$rank = seq(1,nrow(robper.P2.ordered), 1)
robper.P3.ordered$category = "KO"
robper.P3.ordered$plot = "P3"
robper.P3.ordered$rank = seq(1,nrow(robper.P3.ordered), 1)

robper.20gene.overall.ordered$category = "20GT"
robper.20gene.overall.ordered$plot = "Overall"
robper.20gene.overall.ordered$rank = seq(1,nrow(robper.20gene.overall.ordered), 1)
robper.20gene.P1.ordered$category = "20GT"
robper.20gene.P1.ordered$plot = "P1"
robper.20gene.P1.ordered$rank = seq(1,nrow(robper.20gene.P1.ordered), 1)
robper.20gene.P2.ordered$category = "20GT"
robper.20gene.P2.ordered$plot = "P2"
robper.20gene.P2.ordered$rank = seq(1,nrow(robper.20gene.P2.ordered), 1)
robper.20gene.P3.ordered$category = "20GT"
robper.20gene.P3.ordered$plot = "P3"
robper.20gene.P3.ordered$rank = seq(1,nrow(robper.20gene.P3.ordered), 1)

robper.GH.overall.ordered$category = "GH"
robper.GH.overall.ordered$plot = "Overall"
robper.GH.overall.ordered$rank = seq(1,nrow(robper.GH.overall.ordered), 1)
robper.GH.P1.ordered$category = "GH"
robper.GH.P1.ordered$plot = "P1"
robper.GH.P1.ordered$rank = seq(1,nrow(robper.GH.P1.ordered), 1)
robper.GH.P2.ordered$category = "GH"
robper.GH.P2.ordered$plot = "P2"
robper.GH.P2.ordered$rank = seq(1,nrow(robper.GH.P2.ordered), 1)
robper.GH.P3.ordered$category = "GH"
robper.GH.P3.ordered$plot = "P3"
robper.GH.P3.ordered$rank = seq(1,nrow(robper.GH.P3.ordered), 1)
robper.GH.P3.ordered$cat = robper.GH.P3.ordered$cat %>% as.character
robper.GH.P3.ordered$name = robper.GH.P3.ordered$name %>% as.character

lombScargle.plot.df = rbind(
    robper.overall.ordered,
    robper.P1.ordered,
    robper.P2.ordered,
    robper.P3.ordered,
    robper.20gene.overall.ordered,
    robper.20gene.P1.ordered,
    robper.20gene.P2.ordered,
    robper.20gene.P3.ordered,
    robper.GH.overall.ordered,
    robper.GH.P1.ordered,
    robper.GH.P2.ordered,
    robper.GH.P3.ordered

)

rpkSum.df = rbind(
    KO_cat_sum,
    taxa_cat_sum,
    GH_cat_sum
)

colnames(rpkSum.df) = c("cat", "RPK.sum")
lombScargle.plot.df$cat = as.character(lombScargle.plot.df$cat)

lombScargle.plot.df.sums = left_join(lombScargle.plot.df, rpkSum.df)


lombScargle.plot.df.sums$category = factor(lombScargle.plot.df.sums$category, levels = c("KO", "20GT", "GH"))

ggplot(lombScargle.plot.df.sums, aes(x = rank, y = R2)) +
facet_wrap(category~plot, scale = "free_x", ncol = 4) +
geom_point(alpha = 0.25, size = 0.75) +
my_gg_theme +
labs(
x = expression(paste("Rank of R"^2)),
y = expression(paste("R"^2, "of 24-hr sine function"))
) +
theme(
axis.text.x = element_text(angle = 35, hjust = 1)
)

ggplot(lombScargle.plot.df.sums, aes(x = RPK.sum, y = R2)) +
facet_wrap(category~plot, scale = "free_x", ncol = 4) +
scale_x_log10(labels = fancy_scientific) +
geom_point(alpha = 0.25, size = 0.75) +
my_gg_theme +
labs(
x = "Global sequence abundance (reads per kilobase)",
y = expression(paste("R"^2, "of 24-hr sine function"))
) +
theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
)

######################
#Get quantiles

p = c(0.25,0.5,0.75,0.95,0.99,1)

p_names <- map_chr(p, ~paste0(.x*100, "%"))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>%
set_names(nm = p)

p_funs

lombScargle.plot.df.percentiles.wide = lombScargle.plot.df.sums %>%
    group_by(category, plot) %>%
    summarize_at(vars(R2), funs(!!!p_funs))

lombScargle.plot.df.percentiles.long = lombScargle.plot.df.percentiles.wide %>%
    pivot_longer(cols = c(-category,-plot), names_to= "percentile", values_to = "R2")

lombScargle.plot.df.percentiles.long$percentile = as.numeric(lombScargle.plot.df.percentiles.long$percentile) * 100

ggplot(lombScargle.plot.df.percentiles.long, aes(percentile, R2)) +
geom_point() +
geom_line() +
labs(y = expression(paste("R"^2, "of 24-hr sine function"))) +
#scale_y_continuous(limits = c(0,0.8), breaks = c(0,0.25,.5,0.75)) +
scale_x_continuous(breaks = c(25,50,75,100)) +
facet_grid(plot ~ category) +
my_gg_theme

#add 95th percentile to df


lombScargle.plot.df.percentiles.wide = lombScargle.plot.df.sums %>%
group_by(category, plot) %>%
summarize(top_1 = quantile(R2, probs = 0.99, na.rm = T))

#lombScargle.plot.df.sums = na.omit(lombScargle.plot.df.sums)

lombScargle.plot.df.sums = drop_na(lombScargle.plot.df.sums)
lombScargle.plot.df.sums$percentile = vector(mode = "character", length = nrow(lombScargle.plot.df.sums))

for(i in 1:nrow(lombScargle.plot.df.percentiles.wide)){
    lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "percentile"][
        lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "category"] == lombScargle.plot.df.percentiles.wide$category[i] &
        lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "plot"] == lombScargle.plot.df.percentiles.wide$plot[i] &
        lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "R2"] >= lombScargle.plot.df.percentiles.wide$top_1[i]
    ] = "99th"
    
    lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "percentile"][
        lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "category"] == lombScargle.plot.df.percentiles.wide$category[i] &
        lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "plot"] == lombScargle.plot.df.percentiles.wide$plot[i] &
        lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "R2"] < lombScargle.plot.df.percentiles.wide$top_1[i]
    ] = "<99th"
}


lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "percentile"][
lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "category"] == "KO" &
lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "plot"] == "Overall" &
lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "R2"] >= 0.2
] = "99th"



p1 = ggplot(lombScargle.plot.df.sums, aes(x = rank, y = R2, color = percentile)) +
facet_wrap(category~plot, scale = "free_x", ncol = 4) +
geom_point(alpha = 0.75, size = 1) +
my_gg_theme +
labs(
x = expression(paste("Rank of R"^2)),
y = expression(paste("R"^2, "of 24-hr sine function"))
) +
guides(colour = guide_legend(override.aes = list(size=3))) +
scale_color_manual(values = c("99th" = "red", "<99th" = "black")) +
theme(
axis.text.x = element_text(angle = 35, hjust = 1),
legend.title = element_text(size = 18)
)

p2 = ggplot(lombScargle.plot.df.sums, aes(x = RPK.sum, y = R2, color = percentile)) +
facet_wrap(category~plot, scale = "free_x", ncol = 4) +
scale_x_log10(labels = fancy_scientific) +
geom_point(alpha = 0.75, size = 1) +
my_gg_theme +
labs(
x = "Global sequence abundance (reads per kilobase)",
y = expression(paste("R"^2, "of 24-hr sine function"))
) +
guides(colour = guide_legend(override.aes = list(size=3))) +
scale_color_manual(values = c("99th" = "red", "<99th" = "black")) +
theme(
axis.text.x = element_text(angle = 35, hjust = 1),
legend.title = element_text(size = 18)
)

pdf("lomb-scargle_KO_taxa_GH.pdf", width = 12, height = 8)
p2
p1
dev.off()


write.table(
    lombScargle.plot.df.sums %>%
        filter(category == "KO" & percentile == "99th") %>%
        select(plot, name, R2),
    file = "top_99th_KO_lomb-scargle.txt",
    quote = F,
    row.names = F,
    sep = "\t"
)


write.table(
lombScargle.plot.df.sums %>%
filter(category == "20GT" & percentile == "99th") %>%
select(plot, name, R2),
file = "top_99th_taxa_lomb-scargle.txt",
quote = F,
row.names = F,
sep = "\t"
)


write.table(
lombScargle.plot.df.sums %>%
filter(category == "GH" & percentile == "99th") %>%
select(plot, name, R2),
file = "top_99th_GH_lomb-scargle.txt",
quote = F,
row.names = F,
sep = "\t"
)


#################
#99th percentile across all annotations categories

lombScargle.plot.df.percentiles.wide = lombScargle.plot.df.sums %>%
group_by(plot) %>%
summarize(top_1 = quantile(R2, probs = 0.99, na.rm = T))

#lombScargle.plot.df.sums = na.omit(lombScargle.plot.df.sums)

lombScargle.plot.df.sums = drop_na(lombScargle.plot.df.sums)
lombScargle.plot.df.sums$percentile = vector(mode = "character", length = nrow(lombScargle.plot.df.sums))

for(i in 1:nrow(lombScargle.plot.df.percentiles.wide)){
    lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "percentile"][
    lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "plot"] == lombScargle.plot.df.percentiles.wide$plot[i] &
    lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "R2"] >= lombScargle.plot.df.percentiles.wide$top_1[i]
    ] = "99th"
    
    lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "percentile"][
    lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "plot"] == lombScargle.plot.df.percentiles.wide$plot[i] &
    lombScargle.plot.df.sums[,colnames(lombScargle.plot.df.sums) == "R2"] < lombScargle.plot.df.percentiles.wide$top_1[i]
    ] = "<99th"
}




p1 = ggplot(lombScargle.plot.df.sums, aes(x = rank, y = R2, color = percentile)) +
facet_wrap(category~plot, scale = "free_x", ncol = 4) +
geom_point(alpha = 0.75, size = 1) +
my_gg_theme +
labs(
x = expression(paste("Rank of R"^2)),
y = expression(paste("R"^2, "of 24-hr sine function"))
) +
guides(colour = guide_legend(override.aes = list(size=3))) +
scale_color_manual(values = c("99th" = "red", "<99th" = "black")) +
theme(
axis.text.x = element_text(angle = 35, hjust = 1),
legend.title = element_text(size = 18)
)

p2 = ggplot(lombScargle.plot.df.sums, aes(x = RPK.sum, y = R2, color = percentile)) +
facet_wrap(category~plot, scale = "free_x", ncol = 4) +
scale_x_log10(labels = fancy_scientific) +
geom_point(alpha = 0.75, size = 1) +
my_gg_theme +
labs(
x = "Global sequence abundance (reads per kilobase)",
y = expression(paste("R"^2, "of 24-hr sine function"))
) +
guides(colour = guide_legend(override.aes = list(size=3))) +
scale_color_manual(values = c("99th" = "red", "<99th" = "black")) +
theme(
axis.text.x = element_text(angle = 35, hjust = 1),
legend.title = element_text(size = 18)
)

pdf("lomb-scargle_KO_taxa_GH.pdf", width = 12, height = 8)
p2
p1
dev.off()


write.table(
lombScargle.plot.df.sums %>%
filter(category == "KO" & percentile == "99th") %>%
select(plot, name, R2),
file = "top_99th_KO_lomb-scargle.txt",
quote = F,
row.names = F,
sep = "\t"
)


write.table(
lombScargle.plot.df.sums %>%
filter(category == "20GT" & percentile == "99th") %>%
select(plot, name, R2),
file = "top_99th_taxa_lomb-scargle.txt",
quote = F,
row.names = F,
sep = "\t"
)


write.table(
lombScargle.plot.df.sums %>%
filter(category == "GH" & percentile == "99th") %>%
select(plot, name, R2),
file = "top_99th_GH_lomb-scargle.txt",
quote = F,
row.names = F,
sep = "\t"
)

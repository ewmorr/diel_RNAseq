require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)
require(lubridate)
require(rlang)

source("~/ggplot_theme.txt")

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
mapped.metadata$Time.bead.beat = strptime(as.character(mapped.metadata$Time.bead.beat), '%m/%d/%Y %H:%M')


#RNA seq data, ref_len is the total length of reference contigs per category
ref_len = read.table("read_counts.KO.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/KO_rarefied_count.rds")

#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len = data.frame(read_count.rarefied[,1:8], read_count.rarefied[,9:115]/ref_len[,9:115] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

mapped.metadata = mapped.metadata[match(colnames(reads_per_len[9:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

fungal_phyla = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mucoromycota", "Blastocladiomycota", "Zoopagomycota")

reads_per_len.fungi = reads_per_len %>% filter(Phylum %in% fungal_phyla)
reads_per_len.bac = reads_per_len %>% filter(Kindom == "Bacteria")

reads_per_len.fungi.sum = data.frame(sumCat = "Fungi", count = colSums(reads_per_len.fungi[,9:115]), sample = names(reads_per_len.fungi[,9:115]))

reads_per_len.FB.KO = rbind(
    data.frame(sumCat = "Fungi", count = colSums(reads_per_len.fungi[,9:115]), Sample = names(reads_per_len.fungi[,9:115])),
    data.frame(sumCat = "Bacteria", count = colSums(reads_per_len.bac[,9:115]), Sample = names(reads_per_len.bac[,9:115]))
) %>%
pivot_wider(., names_from = "sumCat", values_from = "count") %>%
mutate(FB = Fungi/Bacteria)

reads_per_len.FB.KO.metadata = left_join(reads_per_len.FB.KO, mapped.metadata, by = "Sample")


p1 = ggplot(reads_per_len.FB.KO.metadata, aes(as.POSIXct(Time.bead.beat), FB)) +
geom_point() +
geom_smooth(se = F, color = "dark grey") +
facet_wrap(~Plot, nrow = 3) +
labs(x = "Date time", y = "F:B") +
scale_y_continuous(breaks = c(0,2,4,6)) +
geom_hline(yintercept = 1, linetype = 2) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
labs(title = "A") +
my_gg_theme +
theme(plot.title = element_text(hjust = -0.08))#,
#axis.title.x = element_blank(),
#axis.text.x = element_blank())

require(nlme)

KO.lme = lme(fixed = FB ~ Plot, random = ~1|hours.cumulative.RNA, data = reads_per_len.FB.KO.metadata[,colnames(reads_per_len.FB.KO.metadata) != "Time.bead.beat"] )
anova(KO.lme)
KO.aov = aov(FB ~ Plot,data = reads_per_len.FB.KO.metadata[,colnames(reads_per_len.FB.KO.metadata) != "Time.bead.beat"])
model.tables(KO.aov, "means")

#############
#GHs

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
mapped.metadata$Time.bead.beat = strptime(as.character(mapped.metadata$Time.bead.beat), '%m/%d/%Y %H:%M')

#RNA seq data, ref_len is the total length of reference contigs per category
ref_len = read.table("read_counts.GH.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/GH.phylodist_rarefied_count.rds")

#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len = data.frame(read_count.rarefied[,1:8], read_count.rarefied[,9:115]/ref_len[,9:115] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

mapped.metadata = mapped.metadata[match(colnames(reads_per_len[9:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

fungal_phyla = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mucoromycota", "Blastocladiomycota", "Zoopagomycota")

reads_per_len.fungi = reads_per_len %>% filter(Phylum %in% fungal_phyla)
reads_per_len.bac = reads_per_len %>% filter(Kindom == "Bacteria")

reads_per_len.fungi.sum = data.frame(sumCat = "Fungi", count = colSums(reads_per_len.fungi[,9:115]), sample = names(reads_per_len.fungi[,9:115]))

reads_per_len.FB.GH = rbind(
data.frame(sumCat = "Fungi", count = colSums(reads_per_len.fungi[,9:115]), Sample = names(reads_per_len.fungi[,9:115])),
data.frame(sumCat = "Bacteria", count = colSums(reads_per_len.bac[,9:115]), Sample = names(reads_per_len.bac[,9:115]))
) %>%
pivot_wider(., names_from = "sumCat", values_from = "count") %>%
mutate(FB = Fungi/Bacteria)


reads_per_len.FB.GH.metadata = left_join(reads_per_len.FB.GH, mapped.metadata, by = "Sample")




p2 = ggplot(reads_per_len.FB.GH.metadata, aes(as.POSIXct(Time.bead.beat), FB)) +
geom_point() +
geom_smooth(se = F, color = "dark grey") +
facet_wrap(~Plot, nrow = 3) +
labs(x = "Date time", y = "F:B") +
#scale_y_continuous(breaks = c(0,5,10)) +
geom_hline(yintercept = 1, linetype = 2) +
scale_x_datetime( date_breaks = "8 hours", date_labels = "%m/%d\n%H:%M") +
labs(title = "B") +
my_gg_theme +
theme(
plot.title = element_text(hjust = -0.05),
axis.title.y = element_blank()
)


GH.lme = lme(fixed = FB ~ Plot, random = ~1|hours.cumulative.RNA, data = reads_per_len.FB.metadata[,colnames(reads_per_len.FB.metadata) != "Time.bead.beat"] )
anova(GH.lme)
GH.aov = aov(FB ~ Plot,data = reads_per_len.FB.metadata[,colnames(reads_per_len.FB.metadata) != "Time.bead.beat"])
model.tables(GH.aov, "means")

require(gridExtra)

pdf("Fig2_F-B_date_time.pdf", width = 16, height = 6)
grid.arrange(p1,p2,nrow = 1)
dev.off()


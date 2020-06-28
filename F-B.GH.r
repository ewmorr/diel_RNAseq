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


#RNA seq data, ref_len is the total length of reference contigs per category
ref_len = read.table("read_counts.GH.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/GH.phylodist_rarefied_count.rds")

#To have dplyr interpret args correctly need to !! and parse_expr() (bc of NSE)
#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len = data.frame(read_count.rarefied[,1:8], read_count.rarefied[,9:115]/ref_len[,9:115] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

mapped.metadata = mapped.metadata[match(colnames(reads_per_len[2:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

fungal_phyla = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Mucoromycota", "Blastocladiomycota", "Zoopagomycota")

reads_per_len.fungi = reads_per_len %>% filter(Phylum %in% fungal_phyla)
reads_per_len.bac = reads_per_len %>% filter(Kindom == "Bacteria")

reads_per_len.fungi.sum = data.frame(sumCat = "Fungi", count = colSums(reads_per_len.fungi[,9:115]), sample = names(reads_per_len.fungi[,9:115]))

reads_per_len.FB = rbind(
    data.frame(sumCat = "Fungi", count = colSums(reads_per_len.fungi[,9:115]), Sample = names(reads_per_len.fungi[,9:115])),
    data.frame(sumCat = "Bacteria", count = colSums(reads_per_len.bac[,9:115]), Sample = names(reads_per_len.bac[,9:115]))
) %>%
pivot_wider(., names_from = "sumCat", values_from = "count") %>%
mutate(FB = Fungi/Bacteria)

reads_per_len.FB.metadata = left_join(reads_per_len.FB, mapped.metadata, by = "Sample")

mapped.metadata = mapped.metadata[match(colnames(reads_per_len.FB[2:length(colnames(reads_per_len.FB))]), mapped.metadata$Sample),]



p1 = ggplot(reads_per_len.FB.metadata, aes(TimePoint, FB)) +
geom_point() +
geom_smooth() +
facet_wrap(~Plot) +
labs(x = "Time point", y = "F:B") +
#scale_y_continuous(breaks = c(0,2,4,6)) +
geom_hline(yintercept = 1, linetype = 2) +
my_gg_theme

exclude_timepoints = c(0,1,2,3,4,5)

p2 = ggplot(reads_per_len.FB.metadata %>% filter(!TimePoint %in% exclude_timepoints), aes(timeOfDay.RNA, FB)) +
geom_point() +
geom_smooth() +
facet_wrap(~Plot) +
labs(x = "Time of day", y = "F:B") +
geom_hline(yintercept = 1, linetype = 2) +
my_gg_theme

pdf("FB_GH.pdf", width = 14, height = 6)
grid.arrange(p1, p2, nrow = 2)
dev.off()

p1 = ggplot(reads_per_len.FB.metadata %>% filter(Moisture >= 0), aes(Moisture, FB)) +
geom_point(aes(color = Plot)) +
geom_smooth(method = "lm", formula = y ~ poly(x, 3)) +
#geom_smooth(method = "loess") +
#facet_wrap(~Plot) +
labs(x = "Moisture", y = "F:B") +
#scale_y_continuous(breaks = c(0,2,4,6)) +
scale_colour_manual(values = c("P1" = "#E69F00","P2" = "#56B4E9", "P3" = "#009E73")) +
geom_hline(yintercept = 1, linetype = 2) +
my_gg_theme

ggplot(reads_per_len.FB.metadata %>% filter(Moisture >= 0), aes(Moisture, FB, color = Plot)) +
geom_point() +
geom_smooth(method = "lm") +
#geom_smooth(method = "loess") +
#facet_wrap(~Plot) +
labs(x = "Moisture", y = "F:B") +
scale_y_continuous(breaks = c(0,2,4,6)) +
scale_colour_manual(values = c("P1" = "#E69F00","P2" = "#56B4E9", "P3" = "#009E73")) +
geom_hline(yintercept = 1, linetype = 2) +
my_gg_theme


p2 = ggplot(reads_per_len.FB.metadata, aes(Temperature, FB)) +
geom_point(aes(color = Plot)) +
geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
#geom_smooth(method = "loess") +
labs(x = "Temperature", y = "F:B") +
scale_y_continuous(breaks = c(0,2,4,6)) +
scale_colour_manual(values = c("P1" = "#E69F00","P2" = "#56B4E9", "P3" = "#009E73")) +
geom_hline(yintercept = 1, linetype = 2) +
my_gg_theme

ggplot(reads_per_len.FB.metadata %>% filter(Moisture >= 0), aes(Temperature, FB)) +
geom_point(aes(color = Moisture)) +
geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
#geom_smooth(method = "loess") +
labs(x = "Temperature", y = "F:B") +
scale_y_continuous(breaks = c(0,2,4,6)) +
scale_colour_continuous(trans = "sqrt") +
geom_hline(yintercept = 1, linetype = 2) +
my_gg_theme

ggplot(reads_per_len.FB.metadata %>% filter(Moisture >= 0), aes(Temperature, FB, color = Plot)) +
geom_point() +
geom_smooth(method = "lm") +
#geom_smooth(method = "loess") +
#facet_wrap(~Plot) +
labs(x = "Temperature", y = "F:B") +
scale_y_continuous(breaks = c(0,2,4,6)) +
scale_colour_manual(values = c("P1" = "#E69F00","P2" = "#56B4E9", "P3" = "#009E73")) +
geom_hline(yintercept = 1, linetype = 2) +
my_gg_theme


summary(lm(FB ~ poly(Moisture, 3), data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))

summary(lm(FB ~ Plot/Moisture - 1, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))
summary(lm(FB ~ Moisture*Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))
summary(lm(FB ~ Moisture+Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))

summary(lm(FB ~ Plot/Temperature - 1, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))
summary(lm(FB ~ Temperature*Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))
summary(lm(FB ~ Temperature+Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))


summary(lme(fixed = FB ~ poly(Moisture,3), random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))
summary(lme(fixed = FB ~ Temperature, random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))

anova(lme(fixed = FB ~ poly(Moisture,2), random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))
anova(lme(fixed = FB ~ Temperature, random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0)))

mod1 = lme(fixed = FB ~ Moisture, random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0))
mod2 = lme(fixed = FB ~ poly(Moisture,2), random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0))
mod3 = lme(fixed = FB ~ poly(Moisture,3), random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0))

anova(mod1,mod2)
anova(mod2,mod3)
#third order poly is best fit
anova(mod3)


mod1 = lme(fixed = FB ~ Temperature, random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0))
mod2 = lme(fixed = FB ~ poly(Temperature,2), random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0))
mod3 = lme(fixed = FB ~ poly(Temperature,3), random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0))

anova(mod1,mod2)
anova(mod2,mod3)
#third order poly is best fit
anova(mod3)
mod4 = lme(fixed = FB ~ poly(Temperature,3) + poly(Moisture,3), random = ~1|Plot, data = reads_per_len.FB.metadata %>% filter(Moisture >= 0))

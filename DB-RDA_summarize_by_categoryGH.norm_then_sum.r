require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)
require(reshape2)
require(rlang)
require(lubridate)

source("~/ggplot_theme.txt")

#arguments for summarizing and output
args = commandArgs(T)
catToSum = args[1]
outputFile = args[2]

catToSum = "Category"
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

mapped.metadata = mapped.metadata[match(colnames(reads_per_len[2:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

#DBRDA
reads_per_len.DBRDA = capscale(
    sqrt(
        t(reads_per_len[
            #The indexing term filters out things that only appear in one sample
            rowSums(reads_per_len[,2:length(colnames(reads_per_len))] > 0) >= 2,
            2:length(colnames(reads_per_len))
        ]))
    ~ sin(2*pi*mapped.metadata$hours.cumulative.RNA/24) + cos(2*pi*mapped.metadata$hours.cumulative.RNA/24),
distance = "bray",
metaMDSdist = T
)

anova(reads_per_len.DBRDA, by = "terms")


#DB-RDA with ordistep

#backward and forward model selection with upper scope

samples_to_exclude = c("P1T0_metaT", "P2T0_metaT", "P3T0_metaT",
"P1T1_metaT", "P2T1_metaT", "P3T1_metaT",
"P1T2_metaT", "P2T2_metaT", "P3T2_metaT",
"P1T3_metaT", "P2T3_metaT", "P3T3_metaT",
"P1T4_metaT", "P2T4_metaT", "P3T4_metaT",
"P1T5_metaT", "P2T5_metaT", "P3T5_metaT"
)
reads_per_len.no_outliers = reads_per_len[,!colnames(reads_per_len) %in% samples_to_exclude]

reads_per_len.for_ord = data.frame(t(reads_per_len.no_outliers[
#The indexing term filters outs that only appear in one sample
rowSums(reads_per_len.no_outliers[,2:length(colnames(reads_per_len.no_outliers))] > 0) >= 2,
2:length(colnames(reads_per_len.no_outliers))
])
)

mapped.metadata.no_outliers = filter(mapped.metadata, !Sample %in% samples_to_exclude)

harmonics_phase_shift = data.frame(
X0 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 0),
X1 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 1/24),
X2 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 2/24),
X3 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 3/24),
X4 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 4/24),
X5 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 5/24),
X6 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 6/24),
X7 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 7/24),
X8 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 8/24),
X9 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 9/24),
X10 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 10/24),
X11 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 11/24),
X12 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 12/24),
X13 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 13/24),
X14 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 14/24),
X15 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 15/24),
X16 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 16/24),
X17 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 17/24),
X18 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 18/24),
X19 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 19/24),
X20 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 20/24),
X21 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 21/24),
X22 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 22/24),
X23 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 23/24)
)


reads_per_len.DBRDA.harmonics = capscale(
sqrt(reads_per_len.for_ord) ~ 1,
distance = "bray",
metaMDSdist = T,
data = harmonics_phase_shift
)

phase_mod_res = data.frame(phase = seq(1,24,1), p.val = vector(length = 24, mode = "numeric"), sumSq = vector(length = 24, mode = "numeric"))

for(i in 1:length(colnames(harmonics_phase_shift))){
    print(i)
    reads_per_len.DBRDA.harmonics = capscale(
    sqrt(reads_per_len.for_ord) ~ harmonics_phase_shift[,i],
    distance = "bray",
    metaMDSdist = T,
    )
    temp.var = anova(reads_per_len.DBRDA.harmonics)
    phase_mod_res$phase[i] = i
    phase_mod_res$p.val[i] = temp.var$`Pr(>F)`[1]
    phase_mod_res$sumSq[i] = temp.var$SumOfSqs[1]
}



potential_explanatory_scope = data.frame(plot = mapped.metadata.no_outliers$Plot, moisture = scale(mapped.metadata.no_outliers$Moisture), temperature = scale(mapped.metadata.no_outliers$Temperature), time = (mapped.metadata.no_outliers$hours.cumulative.RNA), harmonic = harmonics_phase_shift[,which.max(phase_mod_res$sumSq)], timeOfDay = mapped.metadata.no_outliers$timeOfDay.RNA
)
#plot(potential_explanatory_scope$diel_sin+potential_explanatory_scope$diel_cos ~ potential_explanatory_scope$timeOfDay)
p1 = ggplot(potential_explanatory_scope, aes(y = harmonic, x = timeOfDay)) +
geom_point() +
labs(x = "time of day") +
my_gg_theme

p2 = ggplot(potential_explanatory_scope, aes(y = harmonic, x = time))+
geom_point() +
labs(x = "hours from experiment start") +
my_gg_theme

pdf("harmonics_overOtime.pdf", width = 8, height = 6)
grid.arrange(p1,p2,ncol=1)
dev.off()
#plot(potential_explanatory_scope$diel_sin ~ potential_explanatory_scope$temperature)
#plot(potential_explanatory_scope$diel_cos ~ potential_explanatory_scope$temperature)

#After ordistep try adding harmonics


potential_explanatory_scope = data.frame(plot = mapped.metadata.no_outliers$Plot, moisture = scale(mapped.metadata.no_outliers$Moisture), temperature = scale(mapped.metadata.no_outliers$Temperature), time = (mapped.metadata.no_outliers$hours.cumulative.RNA), harmonic = harmonics_phase_shift[,which.max(phase_mod_res$sumSq)]#, timeOfDay = mapped.metadata.no_outliers$timeOfDay.RNA
)


reads_per_len.DBRDA = capscale(
sqrt(reads_per_len.for_ord) ~ 1,
distance = "bray",
metaMDSdist = T,
data = potential_explanatory_scope
)


reads_per_len.DBRDA.all = capscale(
sqrt(reads_per_len.for_ord) ~ .,
distance = "bray",
metaMDSdist = T,
data = potential_explanatory_scope
)

reads_per_len.DBRDA.ordistep.highStepHighPerm = ordistep(reads_per_len.DBRDA, scope = formula(reads_per_len.DBRDA.all), direction = "both", steps = 100, permutations = how(nperm = 999))
reads_per_len.DBRDA.ordistep.highStepHighPerm
reads_per_len.DBRDA.ordistep.highStepHighPerm$anova
reads_per_len.DBRDA.ordistep.highStepHighPerm.terms = anova(reads_per_len.DBRDA.ordistep.highStepHighPerm, by = "terms")

tot_sum_sq = sum(reads_per_len.DBRDA.ordistep.highStepHighPerm.terms$SumOfSqs)
expln_sum_sq = sum(reads_per_len.DBRDA.ordistep.highStepHighPerm.terms$SumOfSqs[1:2])
expln_sum_sq/tot_sum_sq
var1_expln_var = reads_per_len.DBRDA.ordistep.highStepHighPerm.terms$SumOfSqs[1]/tot_sum_sq
var2_expln_var = reads_per_len.DBRDA.ordistep.highStepHighPerm.terms$SumOfSqs[2]/tot_sum_sq

#Plot == 0.125096 tot. variance
#moisture == 0.02141625 tot. variance


plot(reads_per_len.DBRDA.ordistep.highStepHighPerm)

potential_explanatory_scope = data.frame(plot = mapped.metadata.no_outliers$Plot, moisture = scale(mapped.metadata.no_outliers$Moisture), temperature = scale(mapped.metadata.no_outliers$Temperature), time = scale(mapped.metadata.no_outliers$hours.cumulative.RNA), harmonic = harmonics_phase_shift[,which.max(phase_mod_res$sumSq)], timeOfDay = mapped.metadata.no_outliers$timeOfDay.RNA, timepoint = mapped.metadata.no_outliers$TimePoint
)

plot.df = data.frame(reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$u[,1:2], potential_explanatory_scope)
plot.df.biplot = data.frame(label = rownames(reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$biplot)[3], CAP1 = reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$biplot[3,1], CAP2 = reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$biplot[3,2])

#eigen vals
(reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$eig / sum(reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$eig))*(expln_sum_sq/tot_sum_sq)

pdf("dbRDA_best_fit_GH.norm_then_sum.pdf", width = 8, height = 6)
ggplot() +
geom_point(data = plot.df, aes(CAP1, CAP2, shape = plot, fill = timeOfDay), size = 3.5, position = position_jitter(height = 0, width = 0.005)) +
geom_segment(data = plot.df.biplot, aes(x = 0, y = 0, xend = CAP1/3, yend = CAP2/3), color = "black", arrow = arrow(length = unit(0.015, "npc"))) +
#scale_color_gradient(low = "#0571b0", high = "#ca0020") +
scale_shape_manual(values = c("P1" = 21,"P2" = 22, "P3" = 23)) +
geom_text(data = plot.df.biplot, aes(CAP1/2.75, CAP2/2.75, label = label), size = 5) +
scale_fill_gradient2(low = "black", high = "black", mid = "white", midpoint = 12, breaks = c(0,12,24), limits = c(0,24)) +
labs(x = "dbRDA axis 1 (9.7% tot. variance)",
y = "dbRDA axis 2 (2.0% tot. variance)",
title = "Best fit dbRDA terms: plot, moisture\nConstrained variance = 0.16, P = 0.001",
shape = "Plot",
fill = "Time of day"
) +
my_gg_theme +
theme(
legend.title = element_text(size = 20, hjust = 0),
plot.title = element_text(size = 18, hjust = 0)
)
dev.off()



#####################
#33% frequency

num_samples = (reads_per_len.no_outliers %>% colnames %>% length)-1
min_samples = floor(num_samples*0.33)

reads_per_len.for_ord = data.frame(t(reads_per_len.no_outliers[
#The indexing term filters outs that appear in at least 29 samples which is 33% frequency
rowSums(reads_per_len.no_outliers[,2:length(colnames(reads_per_len.no_outliers))] > 0) >= min_samples,
2:length(colnames(reads_per_len.no_outliers))
])
)


mapped.metadata.no_outliers = filter(mapped.metadata, !Sample %in% samples_to_exclude)

harmonics_phase_shift = data.frame(
X0 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 0),
X1 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 1/24),
X2 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 2/24),
X3 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 3/24),
X4 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 4/24),
X5 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 5/24),
X6 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 6/24),
X7 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 7/24),
X8 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 8/24),
X9 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 9/24),
X10 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 10/24),
X11 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 11/24),
X12 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 12/24),
X13 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 13/24),
X14 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 14/24),
X15 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 15/24),
X16 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 16/24),
X17 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 17/24),
X18 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 18/24),
X19 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 19/24),
X20 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 20/24),
X21 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 21/24),
X22 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 22/24),
X23 = cos(2*pi*mapped.metadata.no_outliers$hours.cumulative.RNA/24 + 23/24)
)


reads_per_len.DBRDA.harmonics = capscale(
sqrt(reads_per_len.for_ord) ~ 1,
distance = "bray",
metaMDSdist = T,
data = harmonics_phase_shift
)

phase_mod_res = data.frame(phase = seq(1,24,1), p.val = vector(length = 24, mode = "numeric"), sumSq = vector(length = 24, mode = "numeric"))

for(i in 1:length(colnames(harmonics_phase_shift))){
    print(i)
    reads_per_len.DBRDA.harmonics = capscale(
    sqrt(reads_per_len.for_ord) ~ harmonics_phase_shift[,i],
    distance = "bray",
    metaMDSdist = T,
    )
    temp.var = anova(reads_per_len.DBRDA.harmonics)
    phase_mod_res$phase[i] = i
    phase_mod_res$p.val[i] = temp.var$`Pr(>F)`[1]
    phase_mod_res$sumSq[i] = temp.var$SumOfSqs[1]
}



potential_explanatory_scope = data.frame(plot = mapped.metadata.no_outliers$Plot, moisture = scale(mapped.metadata.no_outliers$Moisture), temperature = scale(mapped.metadata.no_outliers$Temperature), time = (mapped.metadata.no_outliers$hours.cumulative.RNA), harmonic = harmonics_phase_shift[,which.max(phase_mod_res$sumSq)], timeOfDay = mapped.metadata.no_outliers$timeOfDay.RNA
)
#plot(potential_explanatory_scope$diel_sin+potential_explanatory_scope$diel_cos ~ potential_explanatory_scope$timeOfDay)
p1 = ggplot(potential_explanatory_scope, aes(y = harmonic, x = timeOfDay)) +
geom_point() +
labs(x = "time of day") +
my_gg_theme

p2 = ggplot(potential_explanatory_scope, aes(y = harmonic, x = time))+
geom_point() +
labs(x = "hours from experiment start") +
my_gg_theme

pdf("harmonics_overOtime.pdf", width = 8, height = 6)
grid.arrange(p1,p2,ncol=1)
dev.off()
#plot(potential_explanatory_scope$diel_sin ~ potential_explanatory_scope$temperature)
#plot(potential_explanatory_scope$diel_cos ~ potential_explanatory_scope$temperature)

#After ordistep try adding harmonics


potential_explanatory_scope = data.frame(plot = mapped.metadata.no_outliers$Plot, moisture = scale(mapped.metadata.no_outliers$Moisture), temperature = scale(mapped.metadata.no_outliers$Temperature), time = (mapped.metadata.no_outliers$hours.cumulative.RNA), harmonic = harmonics_phase_shift[,which.max(phase_mod_res$sumSq)]#, timeOfDay = mapped.metadata.no_outliers$timeOfDay.RNA
)


reads_per_len.DBRDA = capscale(
sqrt(reads_per_len.for_ord) ~ 1,
distance = "bray",
metaMDSdist = T,
data = potential_explanatory_scope
)


reads_per_len.DBRDA.all = capscale(
sqrt(reads_per_len.for_ord) ~ .,
distance = "bray",
metaMDSdist = T,
data = potential_explanatory_scope
)

reads_per_len.DBRDA.ordistep.highStepHighPerm = ordistep(reads_per_len.DBRDA, scope = formula(reads_per_len.DBRDA.all), direction = "both", steps = 100, permutations = how(nperm = 999))
reads_per_len.DBRDA.ordistep.highStepHighPerm
reads_per_len.DBRDA.ordistep.highStepHighPerm$anova
reads_per_len.DBRDA.ordistep.highStepHighPerm.terms = anova(reads_per_len.DBRDA.ordistep.highStepHighPerm, by = "terms")

tot_sum_sq = sum(reads_per_len.DBRDA.ordistep.highStepHighPerm.terms$SumOfSqs)
expln_sum_sq = sum(reads_per_len.DBRDA.ordistep.highStepHighPerm.terms$SumOfSqs[1:2])
expln_sum_sq/tot_sum_sq
var1_expln_var = reads_per_len.DBRDA.ordistep.highStepHighPerm.terms$SumOfSqs[1]/tot_sum_sq
var2_expln_var = reads_per_len.DBRDA.ordistep.highStepHighPerm.terms$SumOfSqs[2]/tot_sum_sq

#Plot == 0.125096 tot. variance
#moisture == 0.02141625 tot. variance


plot(reads_per_len.DBRDA.ordistep.highStepHighPerm)

potential_explanatory_scope = data.frame(plot = mapped.metadata.no_outliers$Plot, moisture = scale(mapped.metadata.no_outliers$Moisture), temperature = scale(mapped.metadata.no_outliers$Temperature), time = scale(mapped.metadata.no_outliers$hours.cumulative.RNA), harmonic = harmonics_phase_shift[,which.max(phase_mod_res$sumSq)], timeOfDay = mapped.metadata.no_outliers$timeOfDay.RNA, timepoint = mapped.metadata.no_outliers$TimePoint
)

plot.df = data.frame(reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$u[,1:2], potential_explanatory_scope)
plot.df.biplot = data.frame(label = rownames(reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$biplot)[3], CAP1 = reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$biplot[3,1], CAP2 = reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$biplot[3,2])

#eigen vals
(reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$eig / sum(reads_per_len.DBRDA.ordistep.highStepHighPerm$CCA$eig))*(expln_sum_sq/tot_sum_sq)

pdf("dbRDA_best_fit_GH.norm_then_sum.33_perc_freq.pdf", width = 8, height = 6)
ggplot() +
geom_point(data = plot.df, aes(CAP1, CAP2, shape = plot, fill = timeOfDay), size = 3.5, position = position_jitter(height = 0, width = 0.005)) +
geom_segment(data = plot.df.biplot, aes(x = 0, y = 0, xend = CAP1/3, yend = CAP2/3), color = "black", arrow = arrow(length = unit(0.015, "npc"))) +
#scale_color_gradient(low = "#0571b0", high = "#ca0020") +
scale_shape_manual(values = c("P1" = 21,"P2" = 22, "P3" = 23)) +
geom_text(data = plot.df.biplot, aes(CAP1/2.75, CAP2/2.75, label = label), size = 5) +
scale_fill_gradient2(low = "black", high = "black", mid = "white", midpoint = 12, breaks = c(0,12,24), limits = c(0,24)) +
labs(x = "dbRDA axis 1 (12.8% tot. variance)",
y = "dbRDA axis 2 (1.9% tot. variance)",
title = "Best fit dbRDA terms: plot, moisture\nConstrained variance = 0.17, P = 0.001",
shape = "Plot",
fill = "Time of day"
) +
my_gg_theme +
theme(
legend.title = element_text(size = 20, hjust = 0),
plot.title = element_text(size = 18, hjust = 0)
)
dev.off()

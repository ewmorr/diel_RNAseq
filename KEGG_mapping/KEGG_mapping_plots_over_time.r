require(tidyverse)
require(RColorBrewer)
source("~/ggplot_theme.txt")

KO_names = read.table("KO_term_names_rarefied_table.txt", sep = "\t")

P1_sig = read.table("harmonic_nls.KO_P1.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P2_sig = read.table("harmonic_nls.KO_P2.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P3_sig = read.table("harmonic_nls.KO_P3.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")

P1_peaks = read.table("harmonic_nls.KO_P1.top_600_R2.peaks.txt", header = T, sep = "\t")
P2_peaks = read.table("harmonic_nls.KO_P2.top_600_R2.peaks.txt", header = T, sep = "\t")
P3_peaks = read.table("harmonic_nls.KO_P3.top_600_R2.peaks.txt", header = T, sep = "\t")

P1_peaks$KO = sub("KO:", "", P1_peaks$KO)
P2_peaks$KO = sub("KO:", "", P2_peaks$KO)
P3_peaks$KO = sub("KO:", "", P3_peaks$KO)

colnames(P1_sig) = c("A", "B", "C", "KO")
colnames(P2_sig) = c("A", "B", "C", "KO")
colnames(P3_sig) = c("A", "B", "C", "KO")

#Filter pathways
P1_sig.path = P1_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")
P2_sig.path = P2_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")
P3_sig.path = P3_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")

P1_sig.peaks = left_join(P1_sig.path, P1_peaks)
P2_sig.peaks = left_join(P2_sig.path, P2_peaks)
P3_sig.peaks = left_join(P3_sig.path, P3_peaks)

P1_sig.peaks$Plot = "Plot 1"
P2_sig.peaks$Plot = "Plot 2"
P3_sig.peaks$Plot = "Plot 3"

P1_P2_P3_sig.peaks = rbind(P1_sig.peaks, P2_sig.peaks) %>%
    rbind(., P3_sig.peaks)

P1_P2_P3_sig.peaks$A = sub("^A[0-9]+ ", "", P1_P2_P3_sig.peaks$A, perl = T)
P1_P2_P3_sig.peaks$B = sub("^B [0-9]+ ", "", P1_P2_P3_sig.peaks$B, perl = T)

#Round peaks to nearest hour and replace 24 with 0 for plot binning
P1_P2_P3_sig.peaks$peaks.hour = round(P1_P2_P3_sig.peaks$peak)
P1_P2_P3_sig.peaks$peaks.hour[P1_P2_P3_sig.peaks$peaks.hour == 24] = 0

#Remove C and unique to get rid of the egregious multimappers in single B cats
P1_P2_P3_sig.peaks.uniques = P1_P2_P3_sig.peaks %>% select(-C)
P1_P2_P3_sig.peaks.uniques = unique(P1_P2_P3_sig.peaks.uniques)


#Longer legend
#If binning half hourly 49 bins splits the data properly. If binning hourly without rounding 25
#better to use data rounded to the nearest hour and geom_bar instead of geom_histogram to avoid binning altogether
#Hourly probably makes the most sense for measurement/sigDig reasons given we measured every two hours
p1 = ggplot(P1_P2_P3_sig.peaks, aes(x = peaks.hour, group = interaction(peaks.hour, A), fill = A)) +
geom_bar() +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(
palette = "Paired",
labels = c("Cellular Processes" = "Cellular\nProcesses", "Environmental Information Processing" = "Environmental\nInformation\nProcessing", "Genetic Information Processing" = "Genetic\nInformation\nProcessing", "Metabolism" = "Metabolism")
) +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24)) +
labs(x = "Time of day", y = "KO count", fill = "A level KEGG") +
my_gg_theme +
theme(
legend.text = element_text(size = 15),
legend.title = element_text(size = 16),
legend.key.height=unit(2, "cm")
) +
guides(fill=guide_legend(override.aes=list(shape=22, size = 3)))

p1 = ggplot(P1_P2_P3_sig.peaks, aes(x = peaks.hour, group = interaction(peaks.hour, A), fill = A)) +
geom_bar() +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(
palette = "Paired",
labels = c("Cellular Processes" = "Cellular\nProcesses", "Environmental Information Processing" = "Environmental\nInformation\nProcessing", "Genetic Information Processing" = "Genetic\nInformation\nProcessing", "Metabolism" = "Metabolism")
) +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24)) +
labs(x = "Time of day", y = "KO count", fill = "A level KEGG", title = "A") +
my_gg_theme +
theme(
plot.title = element_text(hjust = -0.1),
legend.text = element_text(size = 15),
legend.title = element_text(size = 16),
legend.key.height=unit(2, "cm"),
legend.key = element_rect(color = "white")
)

###################
#Panel B plot
#B level matabolism
p2 = ggplot(P1_P2_P3_sig.peaks %>% filter(A == "Metabolism"),
    aes(x = peaks.hour, group = interaction(peaks.hour, B), fill = B)) +
geom_bar() +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(palette = "Paired") +
scale_x_continuous(limits = c(0,24), breaks = c(0, 4, 8, 12, 16, 20)) +
labs(x = "Time of day", y = "KO count", fill = "B level KEGG ––\nwithin Metabolism", title = "B") +
my_gg_theme +
theme(
plot.title = element_text(hjust = -0.125),
legend.text = element_text(size = 15),
legend.title = element_text(size = 16),
legend.key.height=unit(0.75, "cm")
)

#B level Environmental Information Processing
gplot(P1_P2_P3_sig.peaks.uniques %>% filter(A == "Environmental Information Processing"),
aes(x = peaks.hour, group = interaction(peaks.hour, B), fill = B)) +
geom_bar() +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(palette = "Paired") +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24)) +
labs(x = "Time of day", y = "KO count") +
my_gg_theme+
theme(legend.text = element_text(size = 15))

#B level Genetic Information Processing
ggplot(P1_P2_P3_sig.peaks.uniques %>% filter(A == "Genetic Information Processing"),
aes(x = peaks.hour, group = interaction(peaks.hour, B), fill = B)) +
geom_bar() +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(palette = "Paired") +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24)) +
labs(x = "Time of day", y = "KO count") +
my_gg_theme+
theme(legend.text = element_text(size = 15))

#B level Cellular Processes
ggplot(P1_P2_P3_sig.peaks.uniques %>% filter(A == "Cellular Processes"),
aes(x = peaks.hour, group = interaction(peaks.hour, B), fill = B)) +
geom_bar() +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(palette = "Paired") +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24)) +
labs(x = "Time of day", y = "KO count") +
my_gg_theme+
theme(legend.text = element_text(size = 15))

#For adjacent colored lines indicating A level grouping try plotting geom_bar twice and for the second overwrite fill and/or make the bar narrow and line width (size) thick

#All together

#B level Environmental Information Processing

P1_P2_P3_sig.peaks.uniques %>% group_by(A,B) %>% summarize(n()) %>% print(n = Inf)

P1_P2_P3_sig.peaks.uniques$A = factor(P1_P2_P3_sig.peaks.uniques$A, levels = c("Cellular Processes", "Environmental Information Processing", "Genetic Information Processing", "Metabolism"))
Bs_ordered = P1_P2_P3_sig.peaks.uniques %>% group_by(A,B) %>% summarize(n())
P1_P2_P3_sig.peaks.uniques$B = factor(P1_P2_P3_sig.peaks.uniques$B, levels = Bs_ordered$B)

#Changing ordering so the groups of four at A level go with like colors
P1_P2_P3_sig.peaks.uniques.p3_dat = P1_P2_P3_sig.peaks.uniques
P1_P2_P3_sig.peaks.uniques.p3_dat$A = factor(P1_P2_P3_sig.peaks.uniques.p3_dat$A, levels = c("Cellular Processes", "Genetic Information Processing", "Environmental Information Processing"))
P1_P2_P3_sig.peaks.uniques.p3_dat$B = factor(P1_P2_P3_sig.peaks.uniques.p3_dat$B, levels = c(Bs_ordered$B[1:4], Bs_ordered$B[8:11], Bs_ordered$B[7:5]))

P1_P2_P3_sig.peaks.p3_dat = P1_P2_P3_sig.peaks %>% filter(A != "Metabolism")
P1_P2_P3_sig.peaks.p3_dat$A = factor(P1_P2_P3_sig.peaks.p3_dat$A, levels = c("Cellular Processes", "Genetic Information Processing", "Environmental Information Processing"))
P1_P2_P3_sig.peaks.p3_dat$B = factor(P1_P2_P3_sig.peaks.p3_dat$B, levels = c(Bs_ordered$B[1:4], Bs_ordered$B[8:11], Bs_ordered$B[7:5]))


###################
#Panel C
p3 = ggplot(P1_P2_P3_sig.peaks.p3_dat %>%
    filter(
        A == "Environmental Information Processing" |
        A == "Genetic Information Processing" |
        A == "Cellular Processes"
    ),
    aes(x = peaks.hour, group = interaction(peaks.hour, B), fill = B)) +
geom_bar() +
geom_bar(aes(color = A, x = peaks.hour - 0.4), width = 0.08, size = 0.6) +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(palette = "Paired") +
scale_color_manual(values = c("black", "grey60", "grey40")) +#, guide = guide_legend(keywidth = ) +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24)) +
labs(x = "Time of day", y = "KO count", fill = "B level KEGG", color = "A level KEGG", title = "C") +
my_gg_theme+
theme(
plot.title = element_text(hjust = -0.125, ),
legend.text = element_text(size = 15),
legend.title = element_text(size = 16)
) +
guides(
    color = guide_legend(
        order = 1,
        keywidth = unit(0.2, "cm"),
        override.aes = list(fill = c("black", "grey60", "grey40"))
    )
)


#PLOTS FOR PRINTING
require(gridExtra)

p1 = ggplot(P1_P2_P3_sig.peaks, aes(x = peaks.hour, group = interaction(peaks.hour, A), fill = A)) +
geom_bar() +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(palette = "Paired") +
#scale_fill_brewer(
#palette = "Paired",
#labels = c("Cellular Processes" = "Cellular\nProcesses", "Environmental Information Processing" = "Environmental\nInformation\nProcessing", "Genetic Information Processing" = "Genetic\nInformation\nProcessing", "Metabolism" = "Metabolism")
#) +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24)) +
#labs(x = "Time of day", y = "KO count", fill = "A level KEGG", title = "A") +
labs(x = "", y = "KO count", fill = "A level KEGG", title = "A") +
my_gg_theme +
theme(
plot.title = element_text(hjust = -0.135, margin = margin(t = 0, b = -15)),
legend.text = element_text(size = 15),
legend.title = element_text(size = 16),
#legend.key.height=unit(2, "cm"),
legend.key.height=unit(0.75, "cm"),
legend.key.width=unit(0.5, "cm"),
#keywidth = unit(0.75, "cm"),
legend.key = element_rect(color = "white"),
axis.title.x = element_blank()
)

(P1_P2_P3_sig.peaks %>% filter(A == "Metabolism"))$B %>% unique

p2 = ggplot(P1_P2_P3_sig.peaks %>% filter(A == "Metabolism"),
aes(x = peaks.hour, group = interaction(peaks.hour, B), fill = B)) +
geom_bar() +
facet_wrap(~Plot, nrow = 3) +
#scale_fill_brewer(palette = "Paired") +
scale_fill_brewer(palette = "Paired",
    labels = c(
        "Carbohydrate metabolism" = "Carbohydrate metabolism",
        "Biosynthesis of other secondary metabolites" = "Biosynthesis of other\nsecondary metabolites",
        "Lipid metabolism" = "Lipid metabolism",
        "Metabolism of cofactors and vitamins" = "Metabolism of cofactors and vitamins",
        "Amino acid metabolism" = "Amino acid metabolism",
        "Nucleotide metabolism" = "Nucleotide metabolism",
        "Glycan biosynthesis and metabolism" = "Glycan biosynthesis and metabolism",
        "Xenobiotics biodegradation and metabolism" = "Xenobiotics biodegradation\nand metabolism",
        "Metabolism of other amino acids" = "Metabolism of other amino acids",
        "Metabolism of terpenoids and polyketides" = "Metabolism of terpenoids\nand polyketides",
        "Energy metabolism" = "Energy metabolism"
    )
) +
scale_x_continuous(limits = c(0,24), breaks = c(0, 4, 8, 12, 16, 20)) +
scale_y_continuous(breaks = c(0,15,30,45)) +
#labs(x = "Time of day", y = "KO count", fill = "B level KEGG\nwithin Metabolism", title = "B") +
labs(x = "", y = "KO count", fill = "B level KEGG\nwithin Metabolism", title = "B") +
my_gg_theme +
theme(
plot.title = element_text(hjust = -0.135, margin = margin(t = 0, b = -15)),
legend.text = element_text(size = 15),
legend.title = element_text(size = 16),
legend.key.height=unit(1, "cm"),
axis.title.x = element_blank()
)

P1_P2_P3_sig.peaks.p3_dat$fill_color = vector(mode = "character", length = nrow(P1_P2_P3_sig.peaks.p3_dat))

P1_P2_P3_sig.peaks.p3_dat[P1_P2_P3_sig.peaks.p3_dat$A == "Environmental Information Processing", "fill_color"] = "black"
P1_P2_P3_sig.peaks.p3_dat[P1_P2_P3_sig.peaks.p3_dat$A == "Genetic Information Processing", "fill_color"] = "grey60"
P1_P2_P3_sig.peaks.p3_dat[P1_P2_P3_sig.peaks.p3_dat$A == "Cellular Processes", "fill_color"] = "grey40"

p3 = ggplot(P1_P2_P3_sig.peaks.p3_dat %>%
filter(
A == "Environmental Information Processing" |
A == "Genetic Information Processing" |
A == "Cellular Processes"
),
aes(x = peaks.hour, group = interaction(peaks.hour, B), fill = B)) +
geom_bar() +
geom_bar(aes(color = A, x = peaks.hour - 0.41), width = 0.08, size = 0.7) +
facet_wrap(~Plot, nrow = 3) +
scale_fill_brewer(palette = "Paired") +
scale_color_manual(values = c("black", "grey60", "grey40")) +#, guide = guide_legend(keywidth = ) +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24)) +
labs(x = "Time of day", y = "KO count", fill = "B level KEGG", color = "A level KEGG", title = "C") +
my_gg_theme+
theme(
plot.title = element_text(hjust = -0.135, margin = margin(t = 10, b = -15) ),
legend.text = element_text(size = 15),
legend.title = element_text(size = 16)
) +
guides(
color = guide_legend(
order = 1,
keywidth = unit(0.2, "cm"),
override.aes = list(fill = c("black", "grey60", "grey40"))
)
)

pdf("KEGG_TOD_peaks.pdf", width = 30, height = 6)
grid.arrange(p1,p2,p3,ncol = 3, widths = c(0.3,0.415,0.385))
dev.off()


pdf("KEGG_TOD_peaks.smaller.pdf", width = 24, height = 6)
grid.arrange(p1,p2,p3,ncol = 3, widths = c(0.3,0.415,0.385))
dev.off()

pdf("KEGG_TOD_peaks.long.pdf", width = 10, height = 15)
grid.arrange(p1,p2,p3,ncol = 1, heights = c(0.3, 0.35, 0.35))
dev.off()


pdf("KEGG_TOD_peaks.A.pdf", width = 6, height = 4)
p1
dev.off()


pdf("KEGG_TOD_peaks.B.pdf", width = 8, height = 6)
p2
dev.off()


pdf("KEGG_TOD_peaks.C.pdf", width = 8, height = 6)
p3
dev.off()





require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)
require(lubridate)

source("~/ggplot_theme.txt")


#total read amounts table
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
ref_len = read.table("read_counts.GH.readCountRefLen.join", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/GH_rarefied_count.rds")


#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len = data.frame(Category = read_count.rarefied[,1], read_count.rarefied[,2:108]/ref_len[,2:108])
reads_per_len[is.na(reads_per_len)] = 0
#filter by frequency
reads_per_len = reads_per_len[rowSums(reads_per_len[,2:ncol(reads_per_len)] > 0) > 107 *.33,]
mapped.metadata = mapped.metadata[match(colnames(reads_per_len[2:length(colnames(reads_per_len))]), mapped.metadata$Sample),]

reads_per_len.sumCat = reads_per_len %>% group_by(Category) %>% summarize_if(is.numeric,sum,na.rm = TRUE)


###################################
#Read GH names

#look at totals by sumCat
reads_per_len.sumCat.rowSum = data.frame(sumCat = reads_per_len.sumCat$Category, sum = reads_per_len.sumCat[2:length(colnames(reads_per_len.sumCat))] %>% rowSums)
plot(reads_per_len.sumCat.rowSum[order(reads_per_len.sumCat.rowSum$sum),"sum"])

top_ten_GH = (reads_per_len.sumCat.rowSum[order(reads_per_len.sumCat.rowSum$sum, decreasing = T),])[1:71,]
top_ten_GH$global_RA = top_ten_GH$sum/sum(reads_per_len.sumCat.rowSum$sum)

#Use top ten genera list to filter table for NLS fitting
#var names in mapped meta Plot, hours.cumulative.RNA


harmonic_nls.coefs.list = list()
harmonic_nls.peaks.list = list()

for(i in 1:length(top_ten_GH[,1])){

    harmonic_nls.coefs = data.frame(
        plot = c(rep("all", 3), rep("P1", 3), rep("P2", 3), rep("P3", 3)),
        coef = rep(c("amp", "phase", "C"), 4),
        est = vector(mode = "numeric" , length = 12),
        p.val = vector(mode = "numeric" , length = 12),
        stringsAsFactors = F
    )
    harmonic.nls.peaks = data.frame(
    #GH name will be the list index
    peak = vector(mode = "numeric", 3),
    plot = c("P1", "P2", "P3"),
    stringsAsFactors = F
    )

    temp.table = data.frame(
        counts = reads_per_len.sumCat[reads_per_len.sumCat$Category == top_ten_GH[i,1],2:108] %>% unlist,
        Sample = colnames(reads_per_len.sumCat)[2:108]
    )
    temp.table.meta = full_join(temp.table, mapped.metadata, by = "Sample")
    #Calculate time of day column
    temp.table.meta$TOD.temp = temp.table.meta$hours.cumulative.temp + 10.07
    temp.table.meta$TOD.temp[ temp.table.meta$TOD.temp > 24 &  temp.table.meta$TOD.temp <= 48] =  temp.table.meta$TOD.temp[ temp.table.meta$TOD.temp > 24 &  temp.table.meta$TOD.temp <= 48] - 24
    temp.table.meta$TOD.temp[ temp.table.meta$TOD.temp > 48 &  temp.table.meta$TOD.temp <= 72] =  temp.table.meta$TOD.temp[ temp.table.meta$TOD.temp > 48 &  temp.table.meta$TOD.temp <= 72] - 48
    temp.table.meta$TOD.temp[ temp.table.meta$TOD.temp > 72] =  temp.table.meta$TOD.temp[ temp.table.meta$TOD.temp > 72] - 72
    
    #Initial LME
    harmonic_lm = lm(counts ~ sin(2*pi*hours.cumulative.RNA/24) + cos(2*pi*hours.cumulative.RNA/24), data = temp.table.meta)

    #NLS starting params
    #coefs from estimated linear model
    #These will be used for full model plus plotwise models
    b0 <- coef(harmonic_lm)[[1]]
    alpha <- coef(harmonic_lm)[[2]]
    beta <- coef(harmonic_lm)[[3]]
    #r is amplitude and "offset" (or phase) is phi
    r <- sqrt(alpha^2 + beta^2)
    phi <- atan2(beta, alpha)

    #x and y for full model
    y = temp.table.meta$counts
    x.time = temp.table.meta$hours.cumulative.RNA
    full_dat.nls = nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))
    #extract coefs
    #amp
    harmonic_nls.coefs$est[1] = summary(full_dat.nls)$coefficients[1,1]
    harmonic_nls.coefs$p.val[1] = summary(full_dat.nls)$coefficients[1,4]
    #phase
    harmonic_nls.coefs$est[2] = summary(full_dat.nls)$coefficients[2,1]
    harmonic_nls.coefs$p.val[2] = summary(full_dat.nls)$coefficients[2,4]
    #C
    harmonic_nls.coefs$est[3] = summary(full_dat.nls)$coefficients[3,1]
    harmonic_nls.coefs$p.val[3] = summary(full_dat.nls)$coefficients[3,4]

    #P1
    y = (temp.table.meta %>% filter(Plot == "P1"))$counts
    x.time = (temp.table.meta %>% filter(Plot == "P1"))$hours.cumulative.RNA
    P1.nls = nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))
    #extract coefs
    #amp
    harmonic_nls.coefs$est[4] = summary(P1.nls)$coefficients[1,1]
    harmonic_nls.coefs$p.val[4] = summary(P1.nls)$coefficients[1,4]
    #phase
    harmonic_nls.coefs$est[5] = summary(P1.nls)$coefficients[2,1]
    harmonic_nls.coefs$p.val[5] = summary(P1.nls)$coefficients[2,4]
    #C
    harmonic_nls.coefs$est[6] = summary(P1.nls)$coefficients[3,1]
    harmonic_nls.coefs$p.val[6] = summary(P1.nls)$coefficients[3,4]

    #TOD extraction
    TOD.df = data.frame(
        TOD = (temp.table.meta %>% filter(Plot == "P1"))$TOD.temp,
        prediction = predict(P1.nls)
    )
    TOD.df = TOD.df[order(TOD.df$prediction, decreasing = T),]
    top_three = TOD.df[1:3, "TOD"]
    #Need to check in case values span between 11 and 1 pm because average will not work
    if(min(top_three) <= 1 & max(top_three) >= 23){
        top_three[top_three < 1] = top_three[top_three < 1] + 24
    }
    peak_time = sum(top_three)/3
    harmonic.nls.peaks$peak[1] = peak_time

    #P2
    y = (temp.table.meta %>% filter(Plot == "P2"))$counts
    x.time = (temp.table.meta %>% filter(Plot == "P2"))$hours.cumulative.RNA
    P2.nls = nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))
    #amp
    harmonic_nls.coefs$est[7] = summary(P2.nls)$coefficients[1,1]
    harmonic_nls.coefs$p.val[7] = summary(P2.nls)$coefficients[1,4]
    #phase
    harmonic_nls.coefs$est[8] = summary(P2.nls)$coefficients[2,1]
    harmonic_nls.coefs$p.val[8] = summary(P2.nls)$coefficients[2,4]
    #C
    harmonic_nls.coefs$est[9] = summary(P2.nls)$coefficients[3,1]
    harmonic_nls.coefs$p.val[9] = summary(P2.nls)$coefficients[3,4]
    #TOD extraction
    TOD.df = data.frame(
    TOD = (temp.table.meta %>% filter(Plot == "P2"))$TOD.temp,
    prediction = predict(P2.nls)
    )
    TOD.df = TOD.df[order(TOD.df$prediction, decreasing = T),]
    top_three = TOD.df[1:3, "TOD"]
    #Need to check in case values span between 11 and 1 pm because average will not work
    if(min(top_three) <= 1 & max(top_three) >= 23){
        top_three[top_three < 1] = top_three[top_three < 1] + 24
    }
    peak_time = sum(top_three)/3
    harmonic.nls.peaks$peak[2] = peak_time

    #P3
    y = (temp.table.meta %>% filter(Plot == "P3"))$counts
    x.time = (temp.table.meta %>% filter(Plot == "P3"))$hours.cumulative.RNA
    P3.nls = nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))
    #amp
    harmonic_nls.coefs$est[10] = summary(P3.nls)$coefficients[1,1]
    harmonic_nls.coefs$p.val[10] = summary(P3.nls)$coefficients[1,4]
    #phase
    harmonic_nls.coefs$est[11] = summary(P3.nls)$coefficients[2,1]
    harmonic_nls.coefs$p.val[11] = summary(P3.nls)$coefficients[2,4]
    #C
    harmonic_nls.coefs$est[12] = summary(P3.nls)$coefficients[3,1]
    harmonic_nls.coefs$p.val[12] = summary(P3.nls)$coefficients[3,4]
    #TOD extraction
    TOD.df = data.frame(
    TOD = (temp.table.meta %>% filter(Plot == "P3"))$TOD.temp,
    prediction = predict(P3.nls)
    )
    TOD.df = TOD.df[order(TOD.df$prediction, decreasing = T),]
    top_three = TOD.df[1:3, "TOD"]
    #Need to check in case values span between 11 and 1 pm because average will not work
    if(min(top_three) <= 1 & max(top_three) >= 23){
        top_three[top_three < 1] = top_three[top_three < 1] + 24
    }
    peak_time = sum(top_three)/3
    harmonic.nls.peaks$peak[3] = peak_time

    harmonic_nls.coefs.list[[as.character(top_ten_GH[i,1])]] = harmonic_nls.coefs
    harmonic_nls.peaks.list[[as.character(top_ten_GH[i,1])]] = harmonic.nls.peaks
}

df <- do.call("rbind", harmonic_nls.coefs.list)
df$cat = gsub('\\.[[:digit:]]+$', '', rownames(df))
df$cat = gsub('(CBM[[:digit:]]+);NA', '\\1;CBM unknown/multiple targets', df$cat)
df$cat = gsub('cCellulose', 'CBM cellulose', df$cat)
df$cat = gsub('cXylan', 'CBM xylan', df$cat)
df$cat = gsub('cChitin', 'CBM chitin', df$cat)
df$cat = gsub('cStarch', 'CBM starch', df$cat)
df$cat = gsub('Fibronectin;NA', 'Fibronectin;Fibronectin', df$cat)
df$cat = gsub('GH7;NA', 'GH7;Cellulose', df$cat)
df$cat = gsub('GH17;NA', 'GH17;Cellulose', df$cat)
df$cat = gsub('cMixed.Polysaccharides', 'CBM mixed polysaccharides', df$cat)
df$cat = gsub(';NA', ';Unknown/multiple targets', df$cat)
df$cat = gsub('Mixed.Polysaccharides', 'Mixed polysaccharides', df$cat)
df$cat = gsub('Other.Plant.Polysaccharides', 'Mixed polysaccharides', df$cat)
df$cat = gsub('Other.Animal.Polysaccharides', 'Mixed polysaccharides', df$cat)
df$substrate = gsub('.+;(.+)', '\\1', df$cat)

df$fdr.p = 0
df[df$plot == "P1", "fdr.p"] = p.adjust(df[df$plot == "P1", "p.val"], method = "fdr")
df[df$plot == "P2", "fdr.p"] = p.adjust(df[df$plot == "P2", "p.val"], method = "fdr")
df[df$plot == "P3", "fdr.p"] = p.adjust(df[df$plot == "P3", "p.val"], method = "fdr")
df[df$plot == "all", "fdr.p"] = p.adjust(df[df$plot == "all", "p.val"], method = "fdr")

write.table(df, "harmonic_nls.GH.freq_33.txt", row.names = F, sep = "\t", quote = F)

df.2 <- do.call("rbind", harmonic_nls.peaks.list)
df.2$cat = gsub('\\.[[:digit:]]+$', '', rownames(df.2))
df.2$cat = gsub('(CBM[[:digit:]]+);NA', '\\1;CBM unknown/multiple targets', df.2$cat)
df.2$cat = gsub('cCellulose', 'CBM cellulose', df.2$cat)
df.2$cat = gsub('cXylan', 'CBM xylan', df.2$cat)
df.2$cat = gsub('cChitin', 'CBM chitin', df.2$cat)
df.2$cat = gsub('cStarch', 'CBM starch', df.2$cat)
df.2$cat = gsub('Fibronectin;NA', 'Fibronectin;Fibronectin', df.2$cat)
df.2$cat = gsub('GH7;NA', 'GH7;Cellulose', df.2$cat)
df.2$cat = gsub('GH17;NA', 'GH17;Cellulose', df.2$cat)
df.2$cat = gsub('cMixed.Polysaccharides', 'CBM mixed polysaccharides', df.2$cat)
df.2$cat = gsub(';NA', ';Unknown/multiple targets', df.2$cat)
df.2$cat = gsub('Mixed.Polysaccharides', 'Mixed polysaccharides', df.2$cat)
df.2$cat = gsub('Other.Plant.Polysaccharides', 'Mixed polysaccharides', df.2$cat)
df.2$cat = gsub('Other.Animal.Polysaccharides', 'Mixed polysaccharides', df.2$cat)
df.2$substrate = gsub('.+;(.+)', '\\1', df.2$cat)

write.table(df.2, "harmonic_nls.GH.freq_33.peaks.txt", row.names = F, sep = "\t", quote = F)

df.sig = df %>% filter(fdr.p < 0.05 & coef == "amp" & plot != "all")

df.sig %>% nrow
df.sig %>% filter(plot == "P1") %>% nrow
df.sig %>% filter(plot == "P2") %>% nrow
df.sig %>% filter(plot == "P3") %>% nrow

df.sig.peak = left_join(df.sig, df.2, by = c("cat", "plot", "substrate"))

df.sig.peak$peak.hour = round(df.sig.peak$peak)
df.sig.peak$peak.hour[df.sig.peak$peak.hour == 24] = 0

df.sig.peak[df.sig.peak$plot == "P1", "plot"] = "Plot 1"
df.sig.peak[df.sig.peak$plot == "P2", "plot"] = "Plot 2"
df.sig.peak[df.sig.peak$plot == "P3", "plot"] = "Plot 3"

require(RColorBrewer)

p1 = ggplot(df.sig.peak, aes(x = peak.hour, group = interaction(peak.hour, substrate), fill = substrate)) +
geom_bar() +
facet_wrap(~plot, nrow = 3) +
scale_fill_manual(values = c(brewer.pal(12,"Paired"), "dark grey")) +
scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20, 24), limits = c(0,23.5)) +
scale_y_continuous(breaks = c(0,5,10)) +
labs(x = "Time of day", y = "Number of GH genes") +
my_gg_theme +
theme(
legend.text = element_text(size = 15)
)

pdf("GH_sig_diel_TOD_by_plot.pdf", width = 8, height = 5)
p1
dev.off()


#colnames(top_ten_GH) = c("cat", "sum", "global_RA")

#df.sig.peak.taxaRA = left_join(df.sig.peak, top_ten_GH)
#write.table(df.sig.peak.taxaRA, "GH_sig_diel_peaks.txt", row.names = F, quote = F, sep = "\t")



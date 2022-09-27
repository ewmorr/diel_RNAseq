require(tidyverse)
require(gridExtra)
require(scales)
require(RobPer)
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


#READ DATA
#RNA seq data, ref_len is the total length of reference contigs per category
ref_len = read.table("read_counts.KO.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/KO_rarefied_count.rds")

reads_per_len = data.frame(read_count.rarefied[,1:8], read_count.rarefied[,9:115]/ref_len[,9:115] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0


#summarize by KO
reads_per_len.sumCat = reads_per_len %>% group_by(sumCat = KO) %>%
    summarize_if(is.numeric,sum,na.rm = TRUE)
#sum counts for most abundant
KO_cat_sum = data.frame(
"name" = reads_per_len.sumCat$sumCat,
"RPK.sum" = rowSums(reads_per_len.sumCat[,2:ncol(reads_per_len.sumCat)])*1000
)
KO_cat_sum$name = as.character(KO_cat_sum$name)

reads_per_len = as.data.frame(reads_per_len.sumCat)

#KO names
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

reads_per_len.P1 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", as.character(P1.samples$Sample))]
reads_per_len.P2 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat", as.character(P2.samples$Sample))]
reads_per_len.P3 = reads_per_len[,colnames(reads_per_len) %in% c("sumCat",  as.character(P3.samples$Sample))]

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

#overall
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

plot(robper.P3.ordered$R2)

####End Lomb-scargle

#Filter KO tables by top ten percent of robPer R2
quantile(robper.overall.ordered$R2, probs = 0.9)
quantile(robper.P1.ordered$R2, probs = 0.9)
quantile(robper.P2.ordered$R2, probs = 0.9)
quantile(robper.P3.ordered$R2, probs = 0.9)

robper.overall.ordered %>% filter(R2 > 0.2) %>% nrow
robper.P1.ordered %>% filter(R2 > 0.2) %>% nrow
robper.P2.ordered %>% filter(R2 > 0.2) %>% nrow
robper.P3.ordered %>% filter(R2 > 0.2) %>% nrow

#Take top 600 KO in terms of R2. This is ~ top 10%
top_KOs.overall = robper.overall.ordered$cat[1:600]
top_KOs.P1 = robper.P1.ordered$cat[1:600]
top_KOs.P2 = robper.P2.ordered$cat[1:600]
top_KOs.P3 = robper.P3.ordered$cat[1:600]

reads_per_len.overall.top_600 = reads_per_len %>% filter(sumCat %in% top_KOs.overall)
reads_per_len.P1.top_600 = reads_per_len.P1 %>% filter(sumCat %in% top_KOs.P1)
reads_per_len.P2.top_600 = reads_per_len.P2 %>% filter(sumCat %in% top_KOs.P2)
reads_per_len.P3.top_600 = reads_per_len.P3 %>% filter(sumCat %in% top_KOs.P3)

reads_per_len.overall.top_600$sumCat %>% head
reads_per_len.P1.top_600$sumCat %>% head
reads_per_len.P2.top_600$sumCat %>% head
reads_per_len.P3.top_600$sumCat %>% head

#Will run NLS separately for overaall and then each plot

#overall
harmonic_nls.coefs.list.overall = list()
harmonic_nls.peaks.overall = data.frame(
    KO = vector(mode = "character", length(top_KOs.overall)),
    peak = vector(mode = "numeric", length(top_KOs.overall)),
    stringsAsFactors = F
)

for(i in 1:length(top_KOs.overall)){

    harmonic_nls.coefs = data.frame(
        coef = c("amp", "phase", "C"),
        est = vector(mode = "numeric" , length = 3),
        p.val = vector(mode = "numeric" , length = 3),
        stringsAsFactors = F
    )

    temp.table = data.frame(
        counts = reads_per_len[reads_per_len$sumCat == top_KOs.overall[i],2:108] %>% unlist,
        Sample = colnames(reads_per_len)[2:108]
    )
    temp.table.meta = left_join(temp.table, mapped.metadata, by = "Sample")

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
    #Calculate peak
    TOD.df = data.frame(TOD = temp.table.meta$hours.cumulative.temp + 10.07, prediction = predict(full_dat.nls))
    TOD.df$TOD[TOD.df$TOD > 24 & TOD.df$TOD <= 48] = TOD.df$TOD[TOD.df$TOD > 24 & TOD.df$TOD <= 48] - 24
    TOD.df$TOD[TOD.df$TOD > 48 & TOD.df$TOD <= 72] = TOD.df$TOD[TOD.df$TOD > 48 & TOD.df$TOD <= 72] - 48
    TOD.df$TOD[TOD.df$TOD > 72] = TOD.df$TOD[TOD.df$TOD > 72] - 72

    TOD.df = TOD.df[order(TOD.df$prediction, decreasing = T),]
    top_three = TOD.df[1:3, "TOD"]
    #Need to check in case values span between 11 and 1 pm because average will not work
    if(min(top_three) <= 1 & max(top_three) >= 23){
        top_three[top_three < 1] = top_three[top_three < 1] + 24
    }
    peak_time = sum(top_three)/3
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
    #peak

    harmonic_nls.coefs.list.overall[[as.character(top_KOs.overall[i])]] = harmonic_nls.coefs
    harmonic_nls.peaks.overall$KO[i] = as.character(top_KOs.overall[i])
    harmonic_nls.peaks.overall$peak[i] = peak_time
}

df.overall <- do.call("rbind", harmonic_nls.coefs.list.overall)
df.overall$cat = sub('\\.[[:digit:]]+$', '', rownames(df.overall))

df.overall$fdr.p = p.adjust(df.overall$p.val, method = "fdr")

df.overall %>% filter(coef == "amp" & p.val < 0.05) %>% summarize(n())
df.overall %>% filter(coef == "amp" & fdr.p < 0.05) %>% summarize(n())

sig_KOs.overall = (df.overall %>% filter(coef == "amp" & fdr.p < 0.05))$cat
df.overall.sig = df.overall %>% filter(cat %in% sig_KOs.overall)
df.overall.sig.names = left_join(df.overall.sig, KO_names)

write.table(harmonic_nls.peaks.overall, "harmonic_nls.KO_overall.top_600_R2.peaks.txt", row.names = F, sep = "\t")

write.table(df.overall.sig.names, "harmonic_nls.KO_overall.top_600_R2.sig.txt", row.names = F, sep = "\t", quote = T)

df.overall.sig.names.wide = pivot_wider(df.overall.sig.names, names_from = coef, values_from = c(est, p.val, fdr.p) )
write.table(df.overall.sig.names.wide, "harmonic_nls.KO_overall.top_600_R2.sig.wide.txt", row.names = F, sep = "\t", quote = T)
write.table(data.frame(KO = df.overall.sig.names.wide$cat), "harmonic_nls.KO_overall.top_600_R2.sig.KOs.txt", row.names = F, col.names = F, quote = F, sep = "\t")
#end overall


#P1
harmonic_nls.coefs.list.P1 = list()
harmonic_nls.peaks.P1 = data.frame(
KO = vector(mode = "character", length(top_KOs.overall)),
peak = vector(mode = "numeric", length(top_KOs.overall)),
stringsAsFactors = F
)

for(i in 1:length(top_KOs.P1)){
    #print(i)
    harmonic_nls.coefs = data.frame(
    coef = c("amp", "phase", "C"),
    est = vector(mode = "numeric" , length = 3),
    p.val = vector(mode = "numeric" , length = 3),
    stringsAsFactors = F
    )
    
    temp.table = data.frame(
    counts = reads_per_len.P1[reads_per_len.P1$sumCat == top_KOs.P1[i],2:ncol(reads_per_len.P1)] %>% unlist,
    Sample = colnames(reads_per_len.P1)[2:ncol(reads_per_len.P1)]
    )
    temp.table.meta = left_join(temp.table, mapped.metadata, by = "Sample")
    
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
    possibleError = tryCatch(
        nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1)),
        error=function(e) e
    )
    if(inherits(possibleError,"error")){
        harmonic_nls.coefs$est[1] = "no convergence"
        harmonic_nls.coefs$p.val[1] = 1
        #phase
        harmonic_nls.coefs$est[2] = "no convergence"
        harmonic_nls.coefs$p.val[2] = 1
        #C
        harmonic_nls.coefs$est[3] = "no convergence"
        harmonic_nls.coefs$p.val[3] = 1

    }else{
        full_dat.nls = nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))
        #Calculate peak
        TOD.df = data.frame(TOD = temp.table.meta$hours.cumulative.temp + 10.07, prediction = predict(full_dat.nls))
        TOD.df$TOD[TOD.df$TOD > 24 & TOD.df$TOD <= 48] = TOD.df$TOD[TOD.df$TOD > 24 & TOD.df$TOD <= 48] - 24
        TOD.df$TOD[TOD.df$TOD > 48 & TOD.df$TOD <= 72] = TOD.df$TOD[TOD.df$TOD > 48 & TOD.df$TOD <= 72] - 48
        TOD.df$TOD[TOD.df$TOD > 72] = TOD.df$TOD[TOD.df$TOD > 72] - 72
        
        TOD.df = TOD.df[order(TOD.df$prediction, decreasing = T),]
        top_three = TOD.df[1:3, "TOD"]
        #Need to check in case values span between 11 and 1 pm because average will not work
        if(min(top_three) <= 1 & max(top_three) >= 23){
            top_three[top_three < 1] = top_three[top_three < 1] + 24
        }
        peak_time = sum(top_three)/3
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
 
        harmonic_nls.coefs.list.P1[[as.character(top_KOs.P1[i])]] = harmonic_nls.coefs
        harmonic_nls.peaks.P1$KO[i] = as.character(top_KOs.P1[i])
        harmonic_nls.peaks.P1$peak[i] = peak_time
    }
}
write.table(harmonic_nls.peaks.P1, "harmonic_nls.KO_P1.top_600_R2.peaks.txt", row.names = F, sep = "\t", quote = F)

df.P1 <- do.call("rbind", harmonic_nls.coefs.list.P1)
df.P1$cat = sub('\\.[[:digit:]]+$', '', rownames(df.P1))
df.P1$fdr.p = p.adjust(df.P1$p.val, method = "fdr")

df.P1 %>% filter(coef == "amp" & p.val < 0.05) %>% summarize(n())
df.P1 %>% filter(coef == "amp" & fdr.p < 0.05) %>% summarize(n())

sig_KOs.P1 = (df.P1 %>% filter(coef == "amp" & fdr.p < 0.05))$cat
df.P1.sig = df.P1 %>% filter(cat %in% sig_KOs.P1)
df.P1.sig.names = left_join(df.P1.sig, KO_names)

write.table(df.P1.sig.names, "harmonic_nls.KO_P1.top_600_R2.sig.txt", row.names = F, sep = "\t", quote = T)

df.P1.sig.names.wide = pivot_wider(df.P1.sig.names, names_from = coef, values_from = c(est, p.val, fdr.p) )
write.table(df.P1.sig.names.wide, "harmonic_nls.KO_P1.top_600_R2.sig.wide.txt", row.names = F, sep = "\t", quote = T)
write.table(data.frame(KO = df.P1.sig.names.wide$cat), "harmonic_nls.KO_P1.top_600_R2.sig.KOs.txt", row.names = F, col.names = F, quote = F, sep = "\t")

non_sig_KOs.P1 = data.frame(KO = top_KOs.P1) %>% filter(!KO %in% sig_KOs.P1)
write.table(non_sig_KOs.P1, "harmonic_nls.KO_P1.top_600_R2.NOT_sig.KOs.txt", row.names = F, col.names = F, quote = F, sep = "\t")
#end P1


#P2
harmonic_nls.coefs.list.P2 = list()
harmonic_nls.peaks.P2 = data.frame(
KO = vector(mode = "character", length(top_KOs.overall)),
peak = vector(mode = "numeric", length(top_KOs.overall)),
stringsAsFactors = F
)

for(i in 1:length(top_KOs.P2)){
    #print(i)
    harmonic_nls.coefs = data.frame(
    coef = c("amp", "phase", "C"),
    est = vector(mode = "numeric" , length = 3),
    p.val = vector(mode = "numeric" , length = 3),
    stringsAsFactors = F
    )
    
    temp.table = data.frame(
    counts = reads_per_len.P2[reads_per_len.P2$sumCat == top_KOs.P2[i],2:ncol(reads_per_len.P2)] %>% unlist,
    Sample = colnames(reads_per_len.P2)[2:ncol(reads_per_len.P2)]
    )
    temp.table.meta = left_join(temp.table, mapped.metadata, by = "Sample")
    
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
    possibleError = tryCatch(
    nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1)),
    error=function(e) e
    )
    if(inherits(possibleError,"error")){
        harmonic_nls.coefs$est[1] = "no convergence"
        harmonic_nls.coefs$p.val[1] = 1
        #phase
        harmonic_nls.coefs$est[2] = "no convergence"
        harmonic_nls.coefs$p.val[2] = 1
        #C
        harmonic_nls.coefs$est[3] = "no convergence"
        harmonic_nls.coefs$p.val[3] = 1
        
    }else{
        full_dat.nls = nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))
        #Calculate peak
        TOD.df = data.frame(TOD = temp.table.meta$hours.cumulative.temp + 10.07, prediction = predict(full_dat.nls))
        TOD.df$TOD[TOD.df$TOD > 24 & TOD.df$TOD <= 48] = TOD.df$TOD[TOD.df$TOD > 24 & TOD.df$TOD <= 48] - 24
        TOD.df$TOD[TOD.df$TOD > 48 & TOD.df$TOD <= 72] = TOD.df$TOD[TOD.df$TOD > 48 & TOD.df$TOD <= 72] - 48
        TOD.df$TOD[TOD.df$TOD > 72] = TOD.df$TOD[TOD.df$TOD > 72] - 72
        
        TOD.df = TOD.df[order(TOD.df$prediction, decreasing = T),]
        top_three = TOD.df[1:3, "TOD"]
        #Need to check in case values span between 11 and 1 pm because average will not work
        if(min(top_three) <= 1 & max(top_three) >= 23){
            top_three[top_three < 1] = top_three[top_three < 1] + 24
        }
        peak_time = sum(top_three)/3
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

        harmonic_nls.coefs.list.P2[[as.character(top_KOs.P2[i])]] = harmonic_nls.coefs
        harmonic_nls.peaks.P2$KO[i] = as.character(top_KOs.P2[i])
        harmonic_nls.peaks.P2$peak[i] = peak_time
    }
}
write.table(harmonic_nls.peaks.P2, "harmonic_nls.KO_P2.top_600_R2.peaks.txt", row.names = F, sep = "\t", quote = F)

df.P2 <- do.call("rbind", harmonic_nls.coefs.list.P2)
df.P2$cat = sub('\\.[[:digit:]]+$', '', rownames(df.P2))
df.P2$fdr.p = p.adjust(df.P2$p.val, method = "fdr")

df.P2 %>% filter(coef == "amp" & p.val < 0.05) %>% summarize(n())
df.P2 %>% filter(coef == "amp" & fdr.p < 0.05) %>% summarize(n())

sig_KOs.P2 = (df.P2 %>% filter(coef == "amp" & fdr.p < 0.05))$cat
df.P2.sig = df.P2 %>% filter(cat %in% sig_KOs.P2)
df.P2.sig.names = left_join(df.P2.sig, KO_names)

write.table(df.P2.sig.names, "harmonic_nls.KO_P2.top_600_R2.sig.txt", row.names = F, sep = "\t", quote = T)

df.P2.sig.names.wide = pivot_wider(df.P2.sig.names, names_from = coef, values_from = c(est, p.val, fdr.p) )
write.table(df.P2.sig.names.wide, "harmonic_nls.KO_P2.top_600_R2.sig.wide.txt", row.names = F, sep = "\t", quote = T)
write.table(data.frame(KO = df.P2.sig.names.wide$cat), "harmonic_nls.KO_P2.top_600_R2.sig.KOs.txt", row.names = F, col.names = F, quote = F, sep = "\t")
non_sig_KOs.P2 = data.frame(KO = top_KOs.P2) %>% filter(!KO %in% sig_KOs.P2)
write.table(non_sig_KOs.P2, "harmonic_nls.KO_P2.top_600_R2.NOT_sig.KOs.txt", row.names = F, col.names = F, quote = F, sep = "\t")

#end P2


#P3
harmonic_nls.coefs.list.P3 = list()
harmonic_nls.peaks.P3 = data.frame(
KO = vector(mode = "character", length(top_KOs.overall)),
peak = vector(mode = "numeric", length(top_KOs.overall)),
stringsAsFactors = F
)

for(i in 1:length(top_KOs.P3)){
    #print(i)
    harmonic_nls.coefs = data.frame(
    coef = c("amp", "phase", "C"),
    est = vector(mode = "numeric" , length = 3),
    p.val = vector(mode = "numeric" , length = 3),
    stringsAsFactors = F
    )
    
    temp.table = data.frame(
    counts = reads_per_len.P3[reads_per_len.P3$sumCat == top_KOs.P3[i],2:ncol(reads_per_len.P3)] %>% unlist,
    Sample = colnames(reads_per_len.P3)[2:ncol(reads_per_len.P3)]
    )
    temp.table.meta = left_join(temp.table, mapped.metadata, by = "Sample")
    
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
    possibleError = tryCatch(
    nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1)),
    error=function(e) e
    )
    if(inherits(possibleError,"error")){
        harmonic_nls.coefs$est[1] = "no convergence"
        harmonic_nls.coefs$p.val[1] = 1
        #phase
        harmonic_nls.coefs$est[2] = "no convergence"
        harmonic_nls.coefs$p.val[2] = 1
        #C
        harmonic_nls.coefs$est[3] = "no convergence"
        harmonic_nls.coefs$p.val[3] = 1
        
    }else{
        full_dat.nls = nls(y ~ amp * sin(2*pi*x.time/24+phase)+C, start = list(amp = r, phase = phi, C = 1))
        #Calculate peak
        TOD.df = data.frame(TOD = temp.table.meta$hours.cumulative.temp + 10.07, prediction = predict(full_dat.nls))
        TOD.df$TOD[TOD.df$TOD > 24 & TOD.df$TOD <= 48] = TOD.df$TOD[TOD.df$TOD > 24 & TOD.df$TOD <= 48] - 24
        TOD.df$TOD[TOD.df$TOD > 48 & TOD.df$TOD <= 72] = TOD.df$TOD[TOD.df$TOD > 48 & TOD.df$TOD <= 72] - 48
        TOD.df$TOD[TOD.df$TOD > 72] = TOD.df$TOD[TOD.df$TOD > 72] - 72
        
        TOD.df = TOD.df[order(TOD.df$prediction, decreasing = T),]
        top_three = TOD.df[1:3, "TOD"]
        #Need to check in case values span between 11 and 1 pm because average will not work
        if(min(top_three) <= 1 & max(top_three) >= 23){
            top_three[top_three < 1] = top_three[top_three < 1] + 24
        }
        peak_time = sum(top_three)/3
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

        harmonic_nls.coefs.list.P3[[as.character(top_KOs.P3[i])]] = harmonic_nls.coefs
        harmonic_nls.peaks.P3$KO[i] = as.character(top_KOs.P3[i])
        harmonic_nls.peaks.P3$peak[i] = peak_time
    }
}
write.table(harmonic_nls.peaks.P3, "harmonic_nls.KO_P3.top_600_R2.peaks.txt", row.names = F, sep = "\t", quote = F)

df.P3 <- do.call("rbind", harmonic_nls.coefs.list.P3)
df.P3$cat = sub('\\.[[:digit:]]+$', '', rownames(df.P3))
df.P3$fdr.p = p.adjust(df.P3$p.val, method = "fdr")

df.P3 %>% filter(coef == "amp" & p.val < 0.05) %>% summarize(n())
df.P3 %>% filter(coef == "amp" & fdr.p < 0.05) %>% summarize(n())

sig_KOs.P3 = (df.P3 %>% filter(coef == "amp" & fdr.p < 0.05))$cat
df.P3.sig = df.P3 %>% filter(cat %in% sig_KOs.P3)
df.P3.sig.names = left_join(df.P3.sig, KO_names)

write.table(df.P3.sig.names, "harmonic_nls.KO_P3.top_600_R2.sig.txt", row.names = F, sep = "\t", quote = T)

df.P3.sig.names.wide = pivot_wider(df.P3.sig.names, names_from = coef, values_from = c(est, p.val, fdr.p) )
write.table(df.P3.sig.names.wide, "harmonic_nls.KO_P3.top_600_R2.sig.wide.txt", row.names = F, sep = "\t", quote = T)
write.table(data.frame(KO = df.P3.sig.names.wide$cat), "harmonic_nls.KO_P3.top_600_R2.sig.KOs.txt", row.names = F, col.names = F, quote = F, sep = "\t")
non_sig_KOs.P3 = data.frame(KO = top_KOs.P3) %>% filter(!KO %in% sig_KOs.P3)
write.table(non_sig_KOs.P3, "harmonic_nls.KO_P3.top_600_R2.NOT_sig.KOs.txt", row.names = F, col.names = F, quote = F, sep = "\t")

#end P3




require(dplyr)
require(tidyr)

top_12_KO = read.table("top_12_KO.norm_then_sum.txt", sep = "\t", header = F)
colnames(top_12_KO) = c("KO", "count", "name")

P1_sig_KO = read.table("harmonic_nls.KO_P1.top_600_R2.sig.txt", header = T)
P2_sig_KO = read.table("harmonic_nls.KO_P2.top_600_R2.sig.txt", header = T)
P3_sig_KO = read.table("harmonic_nls.KO_P3.top_600_R2.sig.txt", header = T)


P1_sig_KO %>% filter(cat %in%  as.character(top_12_KO$KO))

P2_sig_KO %>% filter(cat %in%  as.character(top_12_KO$KO))

P3_sig_KO %>% filter(cat %in%  as.character(top_12_KO$KO))


#calculate amp-C ratio
P1_sig_KO.wide = P1_sig_KO %>% pivot_wider(names_from = coef, values_from = c(est, p.val, fdr.p))
P1_sig_KO.wide$amp.C = P1_sig_KO.wide$est_amp/P1_sig_KO.wide$est_C
P1_sig_KO.wide$amp.C %>% range
sum(P1_sig_KO.wide$amp.C)/length(P1_sig_KO.wide$amp.C)
P1_sig_KO.wide$amp.C %>% sd

P2_sig_KO.wide = P2_sig_KO %>% pivot_wider(names_from = coef, values_from = c(est, p.val, fdr.p))
P2_sig_KO.wide$amp.C = P2_sig_KO.wide$est_amp/P2_sig_KO.wide$est_C
P2_sig_KO.wide$amp.C %>% range
sum(P2_sig_KO.wide$amp.C)/length(P2_sig_KO.wide$amp.C)
P2_sig_KO.wide$amp.C %>% sd

P3_sig_KO.wide = P3_sig_KO %>% pivot_wider(names_from = coef, values_from = c(est, p.val, fdr.p))
P3_sig_KO.wide$amp.C = P3_sig_KO.wide$est_amp/P3_sig_KO.wide$est_C
P3_sig_KO.wide$amp.C %>% range
sum(P3_sig_KO.wide$amp.C)/length(P3_sig_KO.wide$amp.C)
P3_sig_KO.wide$amp.C %>% sd


#GHs

top_10_GH = read.table("top_10_GH.txt", sep = "\t")

GH_sig_diel = read.table("harmonic_nls.GH.freq_33.txt", header = T, sep = "\t")

GH_sig_diel.wide = GH_sig_diel %>% pivot_wider(names_from = coef, values_from = c(est, p.val, fdr.p))
GH_sig_diel.wide.sig = filter(GH_sig_diel.wide, fdr.p_amp < 0.05 & plot != "all")

GH_sig_diel.wide.sig$amp.C = GH_sig_diel.wide.sig$est_amp/GH_sig_diel.wide.sig$est_C
GH_sig_diel.wide.sig %>% nrow

GH_sig_diel.wide.sig %>% filter(cat %in% top_10_GH$V1) %>% print(n = Inf, width = Inf)

GH_sig_diel.wide.sig$amp.C %>% range
sum(GH_sig_diel.wide.sig$amp.C)/length(GH_sig_diel.wide.sig$amp.C)
sd(GH_sig_diel.wide.sig$amp.C)

GH_sig_diel.wide.sig$amp.C %>% range
sum(GH_sig_diel.wide.sig$amp.C)/length(GH_sig_diel.wide.sig$amp.C)
sd(GH_sig_diel.wide.sig$amp.C)

sum((GH_sig_diel.wide.sig %>% filter(plot == "P1"))$amp.C)/length((GH_sig_diel.wide.sig %>% filter(plot == "P1"))$amp.C)
sd((GH_sig_diel.wide.sig %>% filter(plot == "P1"))$amp.C)

sum((GH_sig_diel.wide.sig %>% filter(plot == "P2"))$amp.C)/length((GH_sig_diel.wide.sig %>% filter(plot == "P2"))$amp.C)
sd((GH_sig_diel.wide.sig %>% filter(plot == "P2"))$amp.C)

sum((GH_sig_diel.wide.sig %>% filter(plot == "P3"))$amp.C)/length((GH_sig_diel.wide.sig %>% filter(plot == "P3"))$amp.C)
sd((GH_sig_diel.wide.sig %>% filter(plot == "P3"))$amp.C)








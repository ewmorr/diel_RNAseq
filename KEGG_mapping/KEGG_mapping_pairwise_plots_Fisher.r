require(tidyverse)
require(RColorBrewer)
source("~/ggplot_theme.txt")

overall_mappings = read.table("harmonic_nls.KO_P1.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P1_mappings = read.table("harmonic_nls.KO_P1.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P2_mappings = read.table("harmonic_nls.KO_P2.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P3_mappings = read.table("harmonic_nls.KO_P3.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")

colnames(overall_mappings) = c("A", "B", "C", "D")
colnames(P1_mappings) = c("A", "B", "C", "D")
colnames(P2_mappings) = c("A", "B", "C", "D")
colnames(P3_mappings) = c("A", "B", "C", "D")

#Fisher test for differences in number of significant between plots
P1_sig_tot = 101
P2_sig_tot = 35
P3_sig_tot = 205

P1_v_P2 = matrix(c(101, 600-101, 35, 600-35), nrow = 2, dimnames = list(c("diel", "not diel"), c("P1", "P2")))
P1_v_P3 = matrix(c(101, 600-101, 205, 600-205), nrow = 2, dimnames = list(c("diel", "not diel"), c("P1", "P2")))
P2_v_P3 = matrix(c(35, 600-35, 205, 600-205), nrow = 2, dimnames = list(c("diel", "not diel"), c("P1", "P2")))

fisher.test(P1_v_P2, alternative = "t") #sig odds ratio greater than 1 (P1 > P2)
fisher.test(P1_v_P3, alternative = "t") #sig odds ratio less than 1 (P3 > P1)
fisher.test(P2_v_P3, alternative = "t") #sig odds ratio less than 1 (P3 > P2)

###############
#A level
P1_mappings %>% filter(A == "A09150 Organismal Systems")
P2_mappings %>% filter(A == "A09150 Organismal Systems")
P3_mappings %>% filter(A == "A09150 Organismal Systems")
#Organismal is for higher level organisms (nervous system, digestive system, etc

P1_A = P1_mappings %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems") %>% group_by(A) %>% summarize(KOs = n())
P2_A = P2_mappings %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems") %>% group_by(A) %>% summarize(KOs = n())
P3_A = P3_mappings %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems") %>% group_by(A) %>% summarize(KOs = n())

P1_A$Plot = "P1"
P2_A$Plot = "P2"
P3_A$Plot = "P3"

P1_P2_A = rbind(P1_A, P2_A)
P1_P2_P3_A = rbind(P1_P2_A, P3_A)
P1_P2_P3_A %>% print(n = Inf)

display.brewer.all()

ggplot(P1_P2_P3_A, aes(x = Plot, y = KOs, fill = A)) +
geom_bar(stat = "identity") +
scale_fill_brewer(palette = "Paired") +
my_gg_theme

#For Fisher test at A level use only pathway levels rather than BRITE
P1_A = P1_mappings %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite") %>% group_by(A) %>% summarize(KOs = n())
P2_A = P2_mappings %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite") %>% group_by(A) %>% summarize(KOs = n())
P3_A = P3_mappings %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite") %>% group_by(A) %>% summarize(KOs = n())
P1_A$Plot = "P1"
P2_A$Plot = "P2"
P3_A$Plot = "P3"

#Joins
P1_P2_A.wide = full_join(P1_A, P2_A, by = c("A"))
P1_P3_A.wide = full_join(P1_A, P3_A, by = c("A"))
P2_P3_A.wide = full_join(P2_A, P3_A, by = c("A"))

P1_v_P2.A_level_fisher = data.frame(
    A = vector(mode = "character", length = 4),
    odds_ratio = vector(mode = "numeric", length = 4),
    p.value = vector(mode = "numeric", length = 4),
    stringsAsFactors = F
)

P1_v_P3.A_level_fisher = data.frame(
    A = vector(mode = "character", length = 4),
    odds_ratio = vector(mode = "numeric", length = 4),
    p.value = vector(mode = "numeric", length = 4),
    stringsAsFactors = F
)

P2_v_P3.A_level_fisher = data.frame(
    A = vector(mode = "character", length = 4),
    odds_ratio = vector(mode = "numeric", length = 4),
    p.value = vector(mode = "numeric", length = 4),
    stringsAsFactors = F
)

for(i in 1:4){
    temp.mat = matrix(
    c(
        P1_P2_A.wide$KOs.x[i],
        sum(P1_P2_A.wide$KOs.x) - P1_P2_A.wide$KOs.x[i],
        P1_P2_A.wide$KOs.y[i],
        sum(P1_P2_A.wide$KOs.y) - P1_P2_A.wide$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P1_P2_A.wide$A[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P1_v_P2.A_level_fisher$A[i] = as.character(P2_P3_A.wide$A[i])
    P1_v_P2.A_level_fisher$odds_ratio[i] = fisher_result$estimate
    P1_v_P2.A_level_fisher$p.value[i] = fisher_result$p.value
}


for(i in 1:4){
    temp.mat = matrix(
    c(
    P1_P3_A.wide$KOs.x[i],
    sum(P1_P3_A.wide$KOs.x) - P1_P3_A.wide$KOs.x[i],
    P1_P3_A.wide$KOs.y[i],
    sum(P1_P3_A.wide$KOs.y) - P1_P3_A.wide$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P1_P3_A.wide$A[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P1_v_P3.A_level_fisher$A[i] = as.character(P2_P3_A.wide$A[i])
    P1_v_P3.A_level_fisher$odds_ratio[i] = fisher_result$estimate
    P1_v_P3.A_level_fisher$p.value[i] = fisher_result$p.value
}

for(i in 1:4){
    temp.mat = matrix(
    c(
    P2_P3_A.wide$KOs.x[i],
    sum(P2_P3_A.wide$KOs.x) - P2_P3_A.wide$KOs.x[i],
    P2_P3_A.wide$KOs.y[i],
    sum(P2_P3_A.wide$KOs.y) - P2_P3_A.wide$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P2_P3_A.wide$A[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P2_v_P3.A_level_fisher$A[i] = as.character(P2_P3_A.wide$A[i])
    P2_v_P3.A_level_fisher$odds_ratio[i] = fisher_result$estimate
    P2_v_P3.A_level_fisher$p.value[i] = fisher_result$p.value
}





###############
#B level
P1_B = P1_mappings %>%
    filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems") %>%
    group_by(A,B) %>%
    summarize(KOs = n())
P2_B = P2_mappings %>%
    filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems") %>%
    group_by(A,B) %>%
    summarize(KOs = n())
P3_B = P3_mappings %>%
    filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems") %>%
    group_by(A,B) %>%
    summarize(KOs = n())


P1_B$Plot = "P1"
P2_B$Plot = "P2"
P3_B$Plot = "P3"

P1_P2_B = rbind(P1_B, P2_B)
P1_P2_P3_B = rbind(P1_P2_B, P3_B)
P1_P2_P3_B %>% filter(A == "A09100 Metabolism") %>% print(n = Inf)
P1_P2_P3_B %>% filter(A == "A09180 Brite Hierarchies") %>% print(n = Inf)

P1_P2_P3_B.full = full_join(P1_B, P2_B) %>% full_join(., P3_B)
write.table(P1_P2_P3_B.full, "P1_P2_P3_B.full.txt", quote = T, sep = "\t", row.names = F, col.names = T)

ggplot(P1_P2_P3_B %>% filter(A == "A09100 Metabolism"), aes(x = Plot, y = KOs, fill = B)) +
geom_bar(stat = "identity") +
scale_fill_brewer(palette = "Paired") +
my_gg_theme

#Perform pairwise joins for Fisher's exact
#Run BRITE and pathways separately

#BRITE
#P1 P2
P1_P2_B.wide.brite = full_join(P1_B, P2_B, by = c("A", "B")) %>%
    filter(A == "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite")

#P1 P3
P1_P3_B.wide.brite = full_join(P1_B, P3_B, by = c("A", "B")) %>%
    filter(A == "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite")

#P2 P3
P2_P3_B.wide.brite = full_join(P2_B, P3_B, by = c("A", "B")) %>%
    filter(A == "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite")

#PATH
#P1 P2
P1_P2_B.wide.path = full_join(P1_B, P2_B, by = c("A", "B")) %>%
    filter(A != "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite")
P1_P2_B.wide.path$KOs.x[is.na(P1_P2_B.wide.path$KOs.x)] = 0
P1_P2_B.wide.path$Plot.x[is.na(P1_P2_B.wide.path$Plot.x)] = "P1"
P1_P2_B.wide.path$KOs.y[is.na(P1_P2_B.wide.path$KOs.y)] = 0
P1_P2_B.wide.path$Plot.y[is.na(P1_P2_B.wide.path$Plot.y)] = "P2"

#P1 P3
P1_P3_B.wide.path = full_join(P1_B, P3_B, by = c("A", "B")) %>%
    filter(A != "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite")
P1_P3_B.wide.path$KOs.x[is.na(P1_P3_B.wide.path$KOs.x)] = 0
P1_P3_B.wide.path$Plot.x[is.na(P1_P3_B.wide.path$Plot.x)] = "P1"
P1_P3_B.wide.path$KOs.y[is.na(P1_P3_B.wide.path$KOs.y)] = 0
P1_P3_B.wide.path$Plot.y[is.na(P1_P3_B.wide.path$Plot.y)] = "P2"

#P2 P3
P2_P3_B.wide.path = full_join(P2_B, P3_B, by = c("A", "B")) %>%
    filter(A != "A09180 Brite Hierarchies" & A != "A09190 Not Included in Pathway or Brite")
P2_P3_B.wide.path$KOs.x[is.na(P2_P3_B.wide.path$KOs.x)] = 0
P2_P3_B.wide.path$Plot.x[is.na(P2_P3_B.wide.path$Plot.x)] = "P1"
P2_P3_B.wide.path$KOs.y[is.na(P2_P3_B.wide.path$KOs.y)] = 0
P2_P3_B.wide.path$Plot.y[is.na(P2_P3_B.wide.path$Plot.y)] = "P2"

####################################################
#Example matrix for Fisher
#The form is P1 cat over P1 not cat, P2 cat over P2 not cat, for example
P1_v_P2.carb_met = matrix(
    c(
        P1_P2_B.wide.path$KOs.x[1],
        sum(P1_P2_B.wide.path$KOs.x) - P1_P2_B.wide.path$KOs.x[1],
        P1_P2_B.wide.path$KOs.y[1],
        sum(P1_P2_B.wide.path$KOs.y) - P1_P2_B.wide.path$KOs.y[1]
    ),
    nrow = 2,
    dimnames = list(c(P1_P2_B.wide.path$B[1], "not"), c("P1", "P2"))
)
P1_v_P2.carb_met.fisher = fisher.test(P1_v_P2.carb_met, alternative = "t")
####################################################

#Prepare for loop
#BRITE
P1_v_P2.B_level_fisher.brite = data.frame(
    B = vector(mode = "character", length = nrow(P1_P2_B.wide.brite)),
    odds_ratio = vector(mode = "numeric", length = nrow(P1_P2_B.wide.brite)),
    p.value = vector(mode = "numeric", length = nrow(P1_P2_B.wide.brite)),
    stringsAsFactors = F
)
P1_v_P3.B_level_fisher.brite = data.frame(
B = vector(mode = "character", length = nrow(P1_P3_B.wide.brite)),
odds_ratio = vector(mode = "numeric", length = nrow(P1_P3_B.wide.brite)),
p.value = vector(mode = "numeric", length = nrow(P1_P3_B.wide.brite)),
stringsAsFactors = F
)
P2_v_P3.B_level_fisher.brite = data.frame(
B = vector(mode = "character", length = nrow(P2_P3_B.wide.brite)),
odds_ratio = vector(mode = "numeric", length = nrow(P2_P3_B.wide.brite)),
p.value = vector(mode = "numeric", length = nrow(P2_P3_B.wide.brite)),
stringsAsFactors = F
)
#PATH
P1_v_P2.B_level_fisher.path = data.frame(
B = vector(mode = "character", length = nrow(P1_P2_B.wide.path)),
odds_ratio = vector(mode = "numeric", length = nrow(P1_P2_B.wide.path)),
p.value = vector(mode = "numeric", length = nrow(P1_P2_B.wide.path)),
stringsAsFactors = F
)
P1_v_P3.B_level_fisher.path = data.frame(
B = vector(mode = "character", length = nrow(P1_P3_B.wide.path)),
odds_ratio = vector(mode = "numeric", length = nrow(P1_P3_B.wide.path)),
p.value = vector(mode = "numeric", length = nrow(P1_P3_B.wide.path)),
stringsAsFactors = F
)
P2_v_P3.B_level_fisher.path = data.frame(
B = vector(mode = "character", length = nrow(P2_P3_B.wide.path)),
odds_ratio = vector(mode = "numeric", length = nrow(P2_P3_B.wide.path)),
p.value = vector(mode = "numeric", length = nrow(P2_P3_B.wide.path)),
stringsAsFactors = F
)

#BRITE loops
for(i in 1:nrow(P1_P2_B.wide.brite)){
    temp.mat = matrix(
    c(
        P1_P2_B.wide.brite$KOs.x[i],
        sum(P1_P2_B.wide.brite$KOs.x) - P1_P2_B.wide.brite$KOs.x[i],
        P1_P2_B.wide.brite$KOs.y[i],
        sum(P1_P2_B.wide.brite$KOs.y) - P1_P2_B.wide.brite$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P1_P2_B.wide.brite$B[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P1_v_P2.B_level_fisher.brite$B[i] = P1_P2_B.wide.brite$B[i]
    P1_v_P2.B_level_fisher.brite$odds_ratio[i] = fisher_result$estimate
    P1_v_P2.B_level_fisher.brite$p.value[i] = fisher_result$p.value
}
P1_v_P2.B_level_fisher.brite

for(i in 1:nrow(P1_P3_B.wide.brite)){
    temp.mat = matrix(
    c(
    P1_P3_B.wide.brite$KOs.x[i],
    sum(P1_P3_B.wide.brite$KOs.x) - P1_P3_B.wide.brite$KOs.x[i],
    P1_P3_B.wide.brite$KOs.y[i],
    sum(P1_P3_B.wide.brite$KOs.y) - P1_P3_B.wide.brite$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P1_P3_B.wide.brite$B[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P1_v_P3.B_level_fisher.brite$B[i] = P1_P3_B.wide.brite$B[i]
    P1_v_P3.B_level_fisher.brite$odds_ratio[i] = fisher_result$estimate
    P1_v_P3.B_level_fisher.brite$p.value[i] = fisher_result$p.value
}
P1_v_P3.B_level_fisher.brite

for(i in 1:nrow(P2_P3_B.wide.brite)){
    temp.mat = matrix(
    c(
    P2_P3_B.wide.brite$KOs.x[i],
    sum(P2_P3_B.wide.brite$KOs.x) - P2_P3_B.wide.brite$KOs.x[i],
    P2_P3_B.wide.brite$KOs.y[i],
    sum(P2_P3_B.wide.brite$KOs.y) - P2_P3_B.wide.brite$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P2_P3_B.wide.brite$B[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P2_v_P3.B_level_fisher.brite$B[i] = P2_P3_B.wide.brite$B[i]
    P2_v_P3.B_level_fisher.brite$odds_ratio[i] = fisher_result$estimate
    P2_v_P3.B_level_fisher.brite$p.value[i] = fisher_result$p.value
}
P2_v_P3.B_level_fisher.brite

#PATH loops
for(i in 1:nrow(P1_P2_B.wide.path)){
    temp.mat = matrix(
    c(
    P1_P2_B.wide.path$KOs.x[i],
    sum(P1_P2_B.wide.path$KOs.x) - P1_P2_B.wide.path$KOs.x[i],
    P1_P2_B.wide.path$KOs.y[i],
    sum(P1_P2_B.wide.path$KOs.y) - P1_P2_B.wide.path$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P1_P2_B.wide.path$B[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P1_v_P2.B_level_fisher.path$B[i] = P1_P2_B.wide.path$B[i]
    P1_v_P2.B_level_fisher.path$odds_ratio[i] = fisher_result$estimate
    P1_v_P2.B_level_fisher.path$p.value[i] = fisher_result$p.value
}
P1_v_P2.B_level_fisher.path

for(i in 1:nrow(P1_P3_B.wide.path)){
    temp.mat = matrix(
    c(
    P1_P3_B.wide.path$KOs.x[i],
    sum(P1_P3_B.wide.path$KOs.x) - P1_P3_B.wide.path$KOs.x[i],
    P1_P3_B.wide.path$KOs.y[i],
    sum(P1_P3_B.wide.path$KOs.y) - P1_P3_B.wide.path$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P1_P3_B.wide.path$B[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P1_v_P3.B_level_fisher.path$B[i] = P1_P3_B.wide.path$B[i]
    P1_v_P3.B_level_fisher.path$odds_ratio[i] = fisher_result$estimate
    P1_v_P3.B_level_fisher.path$p.value[i] = fisher_result$p.value
}
P1_v_P3.B_level_fisher.path

for(i in 1:nrow(P2_P3_B.wide.path)){
    temp.mat = matrix(
    c(
    P2_P3_B.wide.path$KOs.x[i],
    sum(P2_P3_B.wide.path$KOs.x) - P2_P3_B.wide.path$KOs.x[i],
    P2_P3_B.wide.path$KOs.y[i],
    sum(P2_P3_B.wide.path$KOs.y) - P2_P3_B.wide.path$KOs.y[i]
    ),
    nrow = 2,
    dimnames = list(c(P2_P3_B.wide.path$B[i], "not"), c("P1", "P2"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "t")
    P2_v_P3.B_level_fisher.path$B[i] = P2_P3_B.wide.path$B[i]
    P2_v_P3.B_level_fisher.path$odds_ratio[i] = fisher_result$estimate
    P2_v_P3.B_level_fisher.path$p.value[i] = fisher_result$p.value
}
P2_v_P3.B_level_fisher.path


P2_mappings %>% filter(B == "B 09143 Cell growth and death") #P2 > P1 & P3
P1_mappings %>% filter(B == "B 09108 Metabolism of cofactors and vitamins") #P1 > P3
P3_mappings %>% filter(B == "B 09131 Membrane transport") #P3 > P1
P2_mappings %>% filter(B == "B 09132 Signal transduction") #P2 > P3


################################
#Have not run Fisher for C level
################################

###############
###############
#C level
P1_C = P1_mappings %>% filter(A != "A09160 Human Diseases") %>% group_by(A,B,C) %>% summarize(KOs = n())
P2_C = P2_mappings %>% filter(A != "A09160 Human Diseases") %>% group_by(A,B,C) %>% summarize(KOs = n())
P3_C = P3_mappings %>% filter(A != "A09160 Human Diseases") %>% group_by(A,B,C) %>% summarize(KOs = n())

P1_C$Plot = "Plot 1"
P2_C$Plot = "Plot 2"
P3_C$Plot = "Plot 3"

P1_P2_C = rbind(P1_C, P2_C)
P1_P2_P3_C = rbind(P1_P2_C, P3_C)

P1_P2_P3_C %>% filter(A == "A09180 Brite Hierarchies" & B == "B 09181 Protein families: metabolism") %>% print(n = Inf,width = Inf)

P1_P2_P3_C %>% filter(A == "A09100 Metabolism" & B == "B 09101 Carbohydrate metabolism") %>% print(n = Inf,width = Inf)

P1_P2_P3_C %>% filter(A == "A09100 Metabolism") %>% print(n = Inf)
P1_P2_P3_C %>% filter(A == "A09180 Brite Hierarchies") %>% print(n = Inf)


display.brewer.all()

ggplot(P1_P2_P3_B %>% filter(A == "A09100 Metabolism"), aes(x = Plot, y = KOs, fill = B)) +
geom_bar(stat = "identity") +
scale_fill_brewer(palette = "Paired") +
my_gg_theme

##################
#For fishers exact, the question is what to use for the denominator possibliities are the total number of mappings, the number of mapping per A level category, or the number of significant diel signal KOs. The latter is probably best, this does not account for KOs with multiple mappings, but if using B level categories there should not be many multiple mapping per category. The test is then for which categories are enriched relative to the total significant

#Summarize at C and full_join for 0 counts

P1_C = P1_mappings %>% filter(A != "A09160 Human Diseases") %>% group_by(A,B,C) %>% summarize(KOs = n())
P2_C = P2_mappings %>% filter(A != "A09160 Human Diseases") %>% group_by(A,B,C) %>% summarize(KOs = n())
P3_C = P3_mappings %>% filter(A != "A09160 Human Diseases") %>% group_by(A,B,C) %>% summarize(KOs = n())

P1_C$Plot = "Plot 1"
P2_C$Plot = "Plot 2"
P3_C$Plot = "Plot 3"

full_join(P1_C, P2_C, by = c("A", "B", "Plot")) %>% print(n = Inf)



















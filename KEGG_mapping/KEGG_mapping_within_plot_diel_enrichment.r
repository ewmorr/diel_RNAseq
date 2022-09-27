require(tidyverse)
require(RColorBrewer)
source("~/ggplot_theme.txt")

P1_sig = read.table("harmonic_nls.KO_P1.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P2_sig = read.table("harmonic_nls.KO_P2.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P3_sig = read.table("harmonic_nls.KO_P3.top_600_R2.sig.KOs.KEGG_map.txt", header = F, sep = "\t")

P1_NS = read.table("harmonic_nls.KO_P1.top_600_R2.NOT_sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P2_NS = read.table("harmonic_nls.KO_P2.top_600_R2.NOT_sig.KOs.KEGG_map.txt", header = F, sep = "\t")
P3_NS = read.table("harmonic_nls.KO_P3.top_600_R2.NOT_sig.KOs.KEGG_map.txt", header = F, sep = "\t")

colnames(P1_sig) = c("A", "B", "C", "D")
colnames(P2_sig) = c("A", "B", "C", "D")
colnames(P3_sig) = c("A", "B", "C", "D")
colnames(P1_NS) = c("A", "B", "C", "D")
colnames(P2_NS) = c("A", "B", "C", "D")
colnames(P3_NS) = c("A", "B", "C", "D")

#Filter for brite vs pathways

P1_sig.brite = P1_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A == "A09180 Brite Hierarchies")
P2_sig.brite = P2_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A == "A09180 Brite Hierarchies")
P3_sig.brite = P3_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A == "A09180 Brite Hierarchies")

P1_NS.brite = P1_NS %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A == "A09180 Brite Hierarchies")
P2_NS.brite = P2_NS %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A == "A09180 Brite Hierarchies")
P3_NS.brite = P3_NS %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A == "A09180 Brite Hierarchies")


P1_sig.path = P1_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")
P2_sig.path = P2_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")
P3_sig.path = P3_sig %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")

P1_NS.path = P1_NS %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")
P2_NS.path = P2_NS %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")
P3_NS.path = P3_NS %>% filter(A != "A09160 Human Diseases" & A != "A09150 Organismal Systems" & A != "A09190 Not Included in Pathway or Brite" & A != "A09180 Brite Hierarchies")

#########################
#Fisher test comparisons
#A level (PATH only

P1_sig.path.A = P1_sig.path %>% group_by(A) %>% summarize(n = n())
P1_NS.path.A = P1_NS.path %>% group_by(A) %>% summarize(n = n())
P2_sig.path.A = P2_sig.path %>% group_by(A) %>% summarize(n = n())
P2_NS.path.A = P2_NS.path %>% group_by(A) %>% summarize(n = n())
P3_sig.path.A = P3_sig.path %>% group_by(A) %>% summarize(n = n())
P3_NS.path.A = P3_NS.path %>% group_by(A) %>% summarize(n = n())

P1_example_mat = matrix(
    c(
        P1_sig.path.A$n[1],
        sum(P1_sig.path.A$n) - P1_sig.path.A$n[1],
        P1_NS.path.A$n[1],
       sum(P1_NS.path.A$n) - P1_NS.path.A$n[1]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P1_P2_B.wide.path$A[1]), "not"), c("diel", "notDiel"))
)
P1_example_mat
fisher.test(P1_example_mat, alternative = "g")

#LOOPS
#P1
for(i in 1:4){
    temp.mat = matrix(
        c(
            P1_sig.path.A$n[i],
            sum(P1_sig.path.A$n) - P1_sig.path.A$n[i],
            P1_NS.path.A$n[i],
            sum(P1_NS.path.A$n) - P1_NS.path.A$n[i]
        ),
    nrow = 2,
    dimnames = list(c(as.character(P1_sig.path.A$A[i]), "not"), c("diel", "notDiel"))
    )
    print(as.character(P1_sig.path.A$A[i]))
    print(fisher.test(temp.mat, alternative = "g"))
}

for(i in 1:4){
    temp.mat = matrix(
    c(
    P2_sig.path.A$n[i],
    sum(P2_sig.path.A$n) - P2_sig.path.A$n[i],
    P2_NS.path.A$n[i],
    sum(P2_NS.path.A$n) - P2_NS.path.A$n[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P2_sig.path.A$A[i]), "not"), c("diel", "notDiel"))
    )
    print(as.character(P2_sig.path.A$A[i]))
    print(fisher.test(temp.mat, alternative = "g"))
}

for(i in 1:4){
    temp.mat = matrix(
    c(
    P3_sig.path.A$n[i],
    sum(P3_sig.path.A$n) - P3_sig.path.A$n[i],
    P3_NS.path.A$n[i],
    sum(P3_NS.path.A$n) - P3_NS.path.A$n[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P3_sig.path.A$A[i]), "not"), c("diel", "notDiel"))
    )
    print(as.character(P3_sig.path.A$A[i]))
    print(fisher.test(temp.mat, alternative = "g"))
}

#P3 metabolism is the only sig
matrix(
c(
P3_sig.path.A$n[1],
sum(P3_sig.path.A$n) - P3_sig.path.A$n[1],
P3_NS.path.A$n[i],
sum(P3_NS.path.A$n) - P3_NS.path.A$n[1]
),
nrow = 2,
dimnames = list(c(as.character(P3_sig.path.A$A[1]), "not"), c("diel", "notDiel"))
)


#############
#For B level need to take care that NS categories are in same order by index
#(i.e., are all cats present in both)
#This is done by left join and NA=0 lines

#PATH
P1_sig.path.B = P1_sig.path %>% group_by(B) %>% summarize(n = n())
P1_NS.path.B = P1_NS.path %>% group_by(B) %>% summarize(n = n())
P2_sig.path.B = P2_sig.path %>% group_by(B) %>% summarize(n = n())
P2_NS.path.B = P2_NS.path %>% group_by(B) %>% summarize(n = n())
P3_sig.path.B = P3_sig.path %>% group_by(B) %>% summarize(n = n())
P3_NS.path.B = P3_NS.path %>% group_by(B) %>% summarize(n = n())

P1_sig.path.B$sig = "sig"
P1_NS.path.B$sig = "NS"
P2_sig.path.B$sig = "sig"
P2_NS.path.B$sig = "NS"
P3_sig.path.B$sig = "sig"
P3_NS.path.B$sig = "NS"

P1.sigNS.path = left_join(P1_sig.path.B, P1_NS.path.B, by = "B")
P2.sigNS.path = left_join(P2_sig.path.B, P2_NS.path.B, by = "B")
P3.sigNS.path = left_join(P3_sig.path.B, P3_NS.path.B, by = "B")

P3.sigNS.path$sig.y[is.na(P3.sigNS.path$sig.y)] = "NS"
P3.sigNS.path$n.y[is.na(P3.sigNS.path$n.y)] = 0


#BRITE
P1_sig.brite.B = P1_sig.brite %>% group_by(B) %>% summarize(n = n())
P1_NS.brite.B = P1_NS.brite %>% group_by(B) %>% summarize(n = n())
P2_sig.brite.B = P2_sig.brite %>% group_by(B) %>% summarize(n = n())
P2_NS.brite.B = P2_NS.brite %>% group_by(B) %>% summarize(n = n())
P3_sig.brite.B = P3_sig.brite %>% group_by(B) %>% summarize(n = n())
P3_NS.brite.B = P3_NS.brite %>% group_by(B) %>% summarize(n = n())

P1_sig.brite.B$sig = "sig"
P1_NS.brite.B$sig = "NS"
P2_sig.brite.B$sig = "sig"
P2_NS.brite.B$sig = "NS"
P3_sig.brite.B$sig = "sig"
P3_NS.brite.B$sig = "NS"

P1.sigNS.brite = left_join(P1_sig.brite.B, P1_NS.brite.B, by = "B")
P2.sigNS.brite = left_join(P2_sig.brite.B, P2_NS.brite.B, by = "B")
P3.sigNS.brite = left_join(P3_sig.brite.B, P3_NS.brite.B, by = "B")

####################
#Example matrix
####################
temp.mat = matrix(
c(
P1.sigNS.path$n.x[1],
sum(P1.sigNS.path$n.x) - P1.sigNS.path$n.x[1],
P1.sigNS.path$n.y[1],
sum(P1.sigNS.path$n.y) - P1.sigNS.path$n.y[1]
),
nrow = 2,
dimnames = list(c(as.character(P1.sigNS.path$B[1]), "not"), c("diel", "notDiel"))
)
####################

#P1 PATH
P1.B_level.fisher_enrichment.path = data.frame(
B = vector(mode = "character", length = nrow(P1.sigNS.path)),
odds_ratio = vector(mode = "numeric", length = nrow(P1.sigNS.path)),
p.value = vector(mode = "numeric", length = nrow(P1.sigNS.path)),
stringsAsFactors = F
)

for(i in 1:nrow(P1.path.sigNS.path)){
    temp.mat = matrix(
    c(
        P1.sigNS.path$n.x[i],
        sum(P1.sigNS.path$n.x) - P1.sigNS.path$n.x[i],
        P1.sigNS.path$n.y[i],
        sum(P1.sigNS.path$n.y) - P1.sigNS.path$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P1.sigNS.path$B[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P1.B_level.fisher_enrichment.path$B[i] = as.character(P1.sigNS.path$B[i])
    P1.B_level.fisher_enrichment.path$odds_ratio[i] = fisher_result$estimate
    P1.B_level.fisher_enrichment.path$p.value[i] = fisher_result$p.value
}
P1.B_level.fisher_enrichment.path

#P2 PATH
P2.B_level.fisher_enrichment.path = data.frame(
B = vector(mode = "character", length = nrow(P2.sigNS.path)),
odds_ratio = vector(mode = "numeric", length = nrow(P2.sigNS.path)),
p.value = vector(mode = "numeric", length = nrow(P2.sigNS.path)),
stringsAsFactors = F
)

for(i in 1:nrow(P2.sigNS.path)){
    temp.mat = matrix(
    c(
    P2.path.sigNS.path$n.x[i],
    sum(P2.sigNS.path$n.x) - P2.sigNS.path$n.x[i],
    P2.sigNS.path$n.y[i],
    sum(P2.sigNS.path$n.y) - P2.sigNS.path$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P2.sigNS.path$B[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P2.B_level.fisher_enrichment.path$B[i] = as.character(P2.sigNS.path$B[i])
    P2.B_level.fisher_enrichment.path$odds_ratio[i] = fisher_result$estimate
    P2.B_level.fisher_enrichment.path$p.value[i] = fisher_result$p.value
}
P2.B_level.fisher_enrichment.path

#P3 PATH
P3.B_level.fisher_enrichment.path = data.frame(
B = vector(mode = "character", length = nrow(P3.sigNS.path)),
odds_ratio = vector(mode = "numeric", length = nrow(P3.sigNS.path)),
p.value = vector(mode = "numeric", length = nrow(P3.sigNS.path)),
stringsAsFactors = F
)

temp.mat = matrix(
c(
P3.sigNS.path$n.x[1],
sum(P3.sigNS.path$n.x) - P3.sigNS.path$n.x[1],
P3.sigNS.path$n.y[1],
sum(P3.sigNS.path$n.y) - P3.sigNS.path$n.y[1]
),
nrow = 2,
dimnames = list(c(as.character(P3.sigNS.path$B[1]), "not"), c("diel", "notDiel"))
)


for(i in 1:nrow(P3.path.sigNS.path)){
    temp.mat = matrix(
        c(
            P3.sigNS.path$n.x[i],
            sum(P3.sigNS.path$n.x) - P3.sigNS.path$n.x[i],
            P3.sigNS.path$n.y[i],
            sum(P3.sigNS.path$n.y) - P3.sigNS.path$n.y[i]
        ),
        nrow = 2,
        dimnames = list(c(as.character(P3.sigNS.path$B[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P3.B_level.fisher_enrichment.path$B[i] = as.character(P3.sigNS.path$B[i])
    P3.B_level.fisher_enrichment.path$odds_ratio[i] = fisher_result$estimate
    P3.B_level.fisher_enrichment.path$p.value[i] = fisher_result$p.value
}
P3.B_level.fisher_enrichment.path



#P1 brite
P1.B_level.fisher_enrichment.brite = data.frame(
B = vector(mode = "character", length = nrow(P1.sigNS.brite)),
odds_ratio = vector(mode = "numeric", length = nrow(P1.sigNS.brite)),
p.value = vector(mode = "numeric", length = nrow(P1.sigNS.brite)),
stringsAsFactors = F
)

for(i in 1:nrow(P1.sigNS.brite)){
    temp.mat = matrix(
    c(
    P1.sigNS.brite$n.x[i],
    sum(P1.sigNS.brite$n.x) - P1.sigNS.brite$n.x[i],
    P1.sigNS.brite$n.y[i],
    sum(P1.sigNS.brite$n.y) - P1.sigNS.brite$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P1.sigNS.brite$B[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P1.B_level.fisher_enrichment.brite$B[i] = as.character(P1.sigNS.brite$B[i])
    P1.B_level.fisher_enrichment.brite$odds_ratio[i] = fisher_result$estimate
    P1.B_level.fisher_enrichment.brite$p.value[i] = fisher_result$p.value
}
P1.B_level.fisher_enrichment.brite

#P2 brite
P2.B_level.fisher_enrichment.brite = data.frame(
B = vector(mode = "character", length = nrow(P2.sigNS.brite)),
odds_ratio = vector(mode = "numeric", length = nrow(P2.sigNS.brite)),
p.value = vector(mode = "numeric", length = nrow(P2.sigNS.brite)),
stringsAsFactors = F
)

for(i in 1:nrow(P2.sigNS.brite)){
    temp.mat = matrix(
    c(
    P2.sigNS.brite$n.x[i],
    sum(P2.sigNS.brite$n.x) - P2.sigNS.brite$n.x[i],
    P2.sigNS.brite$n.y[i],
    sum(P2.sigNS.brite$n.y) - P2.sigNS.brite$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P2.sigNS.brite$B[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P2.B_level.fisher_enrichment.brite$B[i] = as.character(P2.sigNS.brite$B[i])
    P2.B_level.fisher_enrichment.brite$odds_ratio[i] = fisher_result$estimate
    P2.B_level.fisher_enrichment.brite$p.value[i] = fisher_result$p.value
}
P2.B_level.fisher_enrichment.brite

#P3 brite
P3.B_level.fisher_enrichment.brite = data.frame(
B = vector(mode = "character", length = nrow(P3.sigNS.brite)),
odds_ratio = vector(mode = "numeric", length = nrow(P3.sigNS.brite)),
p.value = vector(mode = "numeric", length = nrow(P3.sigNS.brite)),
stringsAsFactors = F
)

for(i in 1:nrow(P3.sigNS.brite)){
    temp.mat = matrix(
    c(
    P3.sigNS.brite$n.x[i],
    sum(P3.sigNS.brite$n.x) - P3.sigNS.brite$n.x[i],
    P3.sigNS.brite$n.y[i],
    sum(P3.sigNS.brite$n.y) - P3.sigNS.brite$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P3.sigNS.brite$B[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P3.B_level.fisher_enrichment.brite$B[i] = as.character(P3.sigNS.brite$B[i])
    P3.B_level.fisher_enrichment.brite$odds_ratio[i] = fisher_result$estimate
    P3.B_level.fisher_enrichment.brite$p.value[i] = fisher_result$p.value
}
P3.B_level.fisher_enrichment.brite




#######################
#C level

#PATH
P1_sig.path.C = P1_sig.path %>% group_by(A, B, C) %>% summarize(n = n())
P1_NS.path.C = P1_NS.path %>% group_by(A, B, C) %>% summarize(n = n())
P2_sig.path.C = P2_sig.path %>% group_by(A, B, C) %>% summarize(n = n())
P2_NS.path.C = P2_NS.path %>% group_by(A, B, C) %>% summarize(n = n())
P3_sig.path.C = P3_sig.path %>% group_by(A, B, C) %>% summarize(n = n())
P3_NS.path.C = P3_NS.path %>% group_by(A, B, C) %>% summarize(n = n())

P1_sig.path.C$sig = "sig"
P1_NS.path.C$sig = "NS"
P2_sig.path.C$sig = "sig"
P2_NS.path.C$sig = "NS"
P3_sig.path.C$sig = "sig"
P3_NS.path.C$sig = "NS"

P1.sigNS.path.C = left_join(P1_sig.path.C, P1_NS.path.C, by = c("A", "B", "C"))
P2.sigNS.path.C = left_join(P2_sig.path.C, P2_NS.path.C, by = c("A", "B", "C"))
P3.sigNS.path.C = left_join(P3_sig.path.C, P3_NS.path.C, by = c("A", "B", "C"))

P1.sigNS.path.C$sig.y[is.na(P1.sigNS.path.C$sig.y)] = "NS"
P1.sigNS.path.C$n.y[is.na(P1.sigNS.path.C$n.y)] = 0
P2.sigNS.path.C$sig.y[is.na(P2.sigNS.path.C$sig.y)] = "NS"
P2.sigNS.path.C$n.y[is.na(P2.sigNS.path.C$n.y)] = 0
P3.sigNS.path.C$sig.y[is.na(P3.sigNS.path.C$sig.y)] = "NS"
P3.sigNS.path.C$n.y[is.na(P3.sigNS.path.C$n.y)] = 0


#BRITE
P1_sig.brite.C = P1_sig.brite %>% group_by(A, B, C) %>% summarize(n = n())
P1_NS.brite.C = P1_NS.brite %>% group_by(A, B, C) %>% summarize(n = n())
P2_sig.brite.C = P2_sig.brite %>% group_by(A, B, C) %>% summarize(n = n())
P2_NS.brite.C = P2_NS.brite %>% group_by(A, B, C) %>% summarize(n = n())
P3_sig.brite.C = P3_sig.brite %>% group_by(A, B, C) %>% summarize(n = n())
P3_NS.brite.C = P3_NS.brite %>% group_by(A, B, C) %>% summarize(n = n())

P1_sig.brite.C$sig = "sig"
P1_NS.brite.C$sig = "NS"
P2_sig.brite.C$sig = "sig"
P2_NS.brite.C$sig = "NS"
P3_sig.brite.C$sig = "sig"
P3_NS.brite.C$sig = "NS"

P1.sigNS.brite.C = left_join(P1_sig.brite.C, P1_NS.brite.C, by = c("A", "B", "C"))
P2.sigNS.brite.C = left_join(P2_sig.brite.C, P2_NS.brite.C, by = c("A", "B", "C"))
P3.sigNS.brite.C = left_join(P3_sig.brite.C, P3_NS.brite.C, by = c("A", "B", "C"))

P1.sigNS.brite.C$sig.y[is.na(P1.sigNS.brite.C$sig.y)] = "NS"
P1.sigNS.brite.C$n.y[is.na(P1.sigNS.brite.C$n.y)] = 0
P2.sigNS.brite.C$sig.y[is.na(P2.sigNS.brite.C$sig.y)] = "NS"
P2.sigNS.brite.C$n.y[is.na(P2.sigNS.brite.C$n.y)] = 0
P3.sigNS.brite.C$sig.y[is.na(P3.sigNS.brite.C$sig.y)] = "NS"
P3.sigNS.brite.C$n.y[is.na(P3.sigNS.brite.C$n.y)] = 0



#P1 PATH
P1.C_level.fisher_enrichment.path = data.frame(
C = vector(mode = "character", length = nrow(P1.sigNS.path.C)),
odds_ratio = vector(mode = "numeric", length = nrow(P1.sigNS.path.C)),
p.value = vector(mode = "numeric", length = nrow(P1.sigNS.path.C)),
stringsAsFactors = F
)

for(i in 1:nrow(P1.sigNS.path.C)){
    temp.mat = matrix(
    c(
    P1.sigNS.path.C$n.x[i],
    sum(P1.sigNS.path.C$n.x) - P1.sigNS.path.C$n.x[i],
    P1.sigNS.path.C$n.y[i],
    sum(P1.sigNS.path.C$n.y) - P1.sigNS.path.C$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P1.sigNS.path.C$C[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P1.C_level.fisher_enrichment.path$C[i] = as.character(P1.sigNS.path.C$C[i])
    P1.C_level.fisher_enrichment.path$odds_ratio[i] = fisher_result$estimate
    P1.C_level.fisher_enrichment.path$p.value[i] = fisher_result$p.value
}
P1.C_level.fisher_enrichment.path

#P2 PATH
P2.C_level.fisher_enrichment.path = data.frame(
C = vector(mode = "character", length = nrow(P2.sigNS.path.C)),
odds_ratio = vector(mode = "numeric", length = nrow(P2.sigNS.path.C)),
p.value = vector(mode = "numeric", length = nrow(P2.sigNS.path.C)),
stringsAsFactors = F
)

for(i in 1:nrow(P2.sigNS.path.C)){
    temp.mat = matrix(
    c(
    P2.sigNS.path.C$n.x[i],
    sum(P2.sigNS.path.C$n.x) - P2.sigNS.path.C$n.x[i],
    P2.sigNS.path.C$n.y[i],
    sum(P2.sigNS.path.C$n.y) - P2.sigNS.path.C$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P2.sigNS.path.C$C[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P2.C_level.fisher_enrichment.path$C[i] = as.character(P2.sigNS.path.C$C[i])
    P2.C_level.fisher_enrichment.path$odds_ratio[i] = fisher_result$estimate
    P2.C_level.fisher_enrichment.path$p.value[i] = fisher_result$p.value
}
P2.C_level.fisher_enrichment.path

#P3 PATH
P3.C_level.fisher_enrichment.path = data.frame(
C = vector(mode = "character", length = nrow(P3.sigNS.path.C)),
odds_ratio = vector(mode = "numeric", length = nrow(P3.sigNS.path.C)),
p.value = vector(mode = "numeric", length = nrow(P3.sigNS.path.C)),
stringsAsFactors = F
)

for(i in 1:nrow(P3.sigNS.path.C)){
    temp.mat = matrix(
    c(
    P3.sigNS.path.C$n.x[i],
    sum(P3.sigNS.path.C$n.x) - P3.sigNS.path.C$n.x[i],
    P3.sigNS.path.C$n.y[i],
    sum(P3.sigNS.path.C$n.y) - P3.sigNS.path.C$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P3.sigNS.path.C$C[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P3.C_level.fisher_enrichment.path$C[i] = as.character(P3.sigNS.path.C$C[i])
    P3.C_level.fisher_enrichment.path$odds_ratio[i] = fisher_result$estimate
    P3.C_level.fisher_enrichment.path$p.value[i] = fisher_result$p.value
}
P3.C_level.fisher_enrichment.path



#P1 BRITE
P1.C_level.fisher_enrichment.brite = data.frame(
C = vector(mode = "character", length = nrow(P1.sigNS.brite.C)),
odds_ratio = vector(mode = "numeric", length = nrow(P1.sigNS.brite.C)),
p.value = vector(mode = "numeric", length = nrow(P1.sigNS.brite.C)),
stringsAsFactors = F
)

for(i in 1:nrow(P1.sigNS.brite.C)){
    temp.mat = matrix(
    c(
    P1.sigNS.brite.C$n.x[i],
    sum(P1.sigNS.brite.C$n.x) - P1.sigNS.brite.C$n.x[i],
    P1.sigNS.brite.C$n.y[i],
    sum(P1.sigNS.brite.C$n.y) - P1.sigNS.brite.C$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P1.sigNS.brite.C$C[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P1.C_level.fisher_enrichment.brite$C[i] = as.character(P1.sigNS.brite.C$C[i])
    P1.C_level.fisher_enrichment.brite$odds_ratio[i] = fisher_result$estimate
    P1.C_level.fisher_enrichment.brite$p.value[i] = fisher_result$p.value
}
P1.C_level.fisher_enrichment.brite

#P2 BRITE
P2.C_level.fisher_enrichment.brite = data.frame(
C = vector(mode = "character", length = nrow(P2.sigNS.brite.C)),
odds_ratio = vector(mode = "numeric", length = nrow(P2.sigNS.brite.C)),
p.value = vector(mode = "numeric", length = nrow(P2.sigNS.brite.C)),
stringsAsFactors = F
)

for(i in 1:nrow(P2.sigNS.brite.C)){
    temp.mat = matrix(
    c(
    P2.sigNS.brite.C$n.x[i],
    sum(P2.sigNS.brite.C$n.x) - P2.sigNS.brite.C$n.x[i],
    P2.sigNS.brite.C$n.y[i],
    sum(P2.sigNS.brite.C$n.y) - P2.sigNS.brite.C$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P2.sigNS.brite.C$C[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P2.C_level.fisher_enrichment.brite$C[i] = as.character(P2.sigNS.brite.C$C[i])
    P2.C_level.fisher_enrichment.brite$odds_ratio[i] = fisher_result$estimate
    P2.C_level.fisher_enrichment.brite$p.value[i] = fisher_result$p.value
}
P2.C_level.fisher_enrichment.brite

#P3 BRITE
P3.C_level.fisher_enrichment.brite = data.frame(
C = vector(mode = "character", length = nrow(P3.sigNS.brite.C)),
odds_ratio = vector(mode = "numeric", length = nrow(P3.sigNS.brite.C)),
p.value = vector(mode = "numeric", length = nrow(P3.sigNS.brite.C)),
stringsAsFactors = F
)

for(i in 1:nrow(P3.sigNS.brite.C)){
    temp.mat = matrix(
    c(
    P3.sigNS.brite.C$n.x[i],
    sum(P3.sigNS.brite.C$n.x) - P3.sigNS.brite.C$n.x[i],
    P3.sigNS.brite.C$n.y[i],
    sum(P3.sigNS.brite.C$n.y) - P3.sigNS.brite.C$n.y[i]
    ),
    nrow = 2,
    dimnames = list(c(as.character(P3.sigNS.brite.C$C[i]), "not"), c("diel", "notDiel"))
    )
    fisher_result = fisher.test(temp.mat, alternative = "g")
    P3.C_level.fisher_enrichment.brite$C[i] = as.character(P3.sigNS.brite.C$C[i])
    P3.C_level.fisher_enrichment.brite$odds_ratio[i] = fisher_result$estimate
    P3.C_level.fisher_enrichment.brite$p.value[i] = fisher_result$p.value
}
P3.C_level.fisher_enrichment.brite






require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)
require(lubridate)
require(gtable)

source("~/ggplot_theme.txt")

setwd("../../martiny_diel_seqs/") 

#########################
#Read count data
#########################

#total read amounts table
total_mapped = read.csv("totals_mapped_input.reformat.txt", sep = "\t", header = F)
colnames(total_mapped) = c("Sample", "totalReads", "readsMappedCov", "totalBases", "targetsLength", "basesMapped", "avgCov")
total_mapped = total_mapped %>% as_tibble
total_mapped = total_mapped %>% filter(Sample != "P2T31_metaT")

#metadata table
metadata = read.table("metadata.txt", header = T, sep = "\t")
metadata = metadata %>% as.tbl
metadata = metadata %>% filter(Sample != "P2T31_metaT")

metadata$timeOfDay.RNA = (as.numeric(as.POSIXlt(metadata$timeOfDay.RNA, format = "%H:%M") %>% hour) + as.numeric((as.POSIXlt(metadata$timeOfDay.RNA, format = "%H:%M") %>% minute)/60) ) %>% signif(2)

mapped.metadata = full_join(total_mapped, metadata, by = "Sample")

#RNA seq data, ref_len is the total length of reference contigs per category
ref_len.20_gene_taxonomy = read.table("read_counts.20_gene_taxonomy.metaG.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len.20_gene_taxonomy$P2T31_metaT = NULL
ref_len.20_gene_taxonomy$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/20_gene_rarefied_counts.metaG.rds")

#Divide read count by refLen to normalize reads recruited by contig length
reads_per_len = data.frame(read_count.rarefied[,1:4], read_count.rarefied[,5:111]/ref_len.20_gene_taxonomy[,5:111] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

nrow(reads_per_len)
reads_per_len$Edge_num %>% unique %>% length

#split tax_strings
rank_mat = data.frame(reads_per_len$Tax_string %>% str_split(., ";", simplify = T))[,1:5]
colnames(rank_mat) = c("Phylum", "Class", "Order", "Family", "Genus")
reads_per_len.tax = data.frame(
  Edge_num = reads_per_len$Edge_num,
  rank_mat[, c("Phylum", "Genus")],
  reads_per_len[,5:ncol(reads_per_len)]
  )

reads_per_len = reads_per_len.tax

reads_per_len$Edge_num = reads_per_len$Edge_num %>% as.factor
reads_per_len$Phylum = reads_per_len$Phylum %>% as.factor
reads_per_len$Genus = reads_per_len$Genus %>% as.factor

############################
#Checking counts on summary
############################

#how many tax names and how many edge nums
reads_per_len$Edge_num %>% unique %>% length
#1773
reads_per_len$Genus %>% unique %>% length
#245
reads_per_len$Phylum %>% unique %>% length
#38

#look at totals by Genus
reads_per_len.Genus = reads_per_len %>% group_by(Genus) %>% summarize_if(is.numeric,sum,na.rm = TRUE)
nrow(reads_per_len.Genus)
#245
reads_per_len.Genus.rowSum = data.frame(Genus = reads_per_len.Genus$Genus, sum = reads_per_len.Genus[2:length(colnames(reads_per_len.Genus))] %>% rowSums)
reads_per_len.Genus.rowSum[order(reads_per_len.Genus.rowSum$sum),]

top_50_genera = (reads_per_len.Genus.rowSum[order(reads_per_len.Genus.rowSum$sum, decreasing = T),]) %>%
  filter(Genus != "NA") %>% #Not sure why the NA is coming through as text but it is
  .[1:50,]
top_50_genera$global_RA = top_50_genera$sum/sum(reads_per_len.Genus.rowSum$sum)

#rejoin with Phylum to get taxonomic names
#retaining phylum in group_by has odd entries

rank_mat.uni = rank_mat %>% unique
#It looks like the proteobacteria clades are not listed at phylum level but instead at class level
rank_mat.uni %>% filter(Phylum == 1224) %>% select(Phylum, Class) %>% unique

#######################
#Write tables for editing and adding names from NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy)
top_50_genera.genus_strings = top_50_genera %>% pull(Genus)
top_50_genera.phylum_strings = left_join(top_50_genera, rank_mat.uni %>% select(Genus, Phylum, Class), by = "Genus") %>% select(Phylum, Class) %>% unique
top_50_genera.tax_strings = left_join(top_50_genera, rank_mat.uni %>% select(Genus, Phylum, Class), by = "Genus") %>% select(Phylum, Class, Genus)

write.table(top_50_genera.genus_strings, "top_50_genera.genus_strings.txt", quote = F, sep = "\t", row.names = F, col.names = F)
#note that only the 1224 phylum needs to becategorized by class
write.table(top_50_genera.phylum_strings, "top_50_genera.phylum_strings.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(top_50_genera.tax_strings, "top_50_genera.tax_strings.txt", quote = F, sep = "\t", row.names = F, col.names = F)

plot(top_50_genera$global_RA)

#Use top 50 taxa list to filter table for NLS fitting
#var names in mapped meta Plot, hours.cumulative.RNA



#Stick with "Edge_num" summary for now
top_50_OTUs = (reads_per_len.Edge_num.rowSum[order(reads_per_len.Edge_num.rowSum$sum, decreasing = T),])[1:50,]
top_50_OTUs$global_RA = top_50_OTUs$sum/sum(reads_per_len.Tax_name.rowSum$sum)

write.table(top_50_OTUs, "top_50_OTUs.txt", quote = F, sep = "\t", row.names = F)

plot(top_50_OTUs$global_RA)

top_50_OTUs$cat = paste(top_50_OTUs$Tax_name, top_50_OTUs$Edge_num, sep = ".")

#if wanting to summarize by a specific taxonomic level (e.g., genus) can look at markers2taxonomy.txt and use the taxnomicString

############################################
#Filter summarized tax_name table to top 50
############################################

reads_per_len.Tax_name.top_50 = reads_per_len.Tax_name %>% filter(Tax_name %in% top_50_genera$Tax_name)
reads_per_len.Edge_num.top_50 = reads_per_len.Edge_num %>% filter(Edge_num %in% top_50_OTUs$Edge_num)

#For lowess, split into groups of samples per plot. Can just do this instead with select(starts_with())
#P1_samples = mapped.metadata %>% filter(Plot == "P1") %>% pull(Sample)
#P2_samples = mapped.metadata %>% filter(Plot == "P2") %>% pull(Sample)
#P3_samples = mapped.metadata %>% filter(Plot == "P3") %>% pull(Sample)

##########################################
#Overall PERMANOVA. performed on Edge_num
##########################################

P1_times = mapped.metadata %>% filter(Plot == "P1") %>% pull(hours.cumulative.RNA)
P2_times = mapped.metadata %>% filter(Plot == "P2") %>% pull(hours.cumulative.RNA)
P3_times = mapped.metadata %>% filter(Plot == "P3") %>% pull(hours.cumulative.RNA)

reads_per_len.Edge_num.P1 = reads_per_len.Edge_num[,3:ncol(reads_per_len.Edge_num)] %>% select(starts_with("P1"))
reads_per_len.Edge_num.P2 = reads_per_len.Edge_num[,3:ncol(reads_per_len.Edge_num)] %>% select(starts_with("P2"))
reads_per_len.Edge_num.P3 = reads_per_len.Edge_num[,3:ncol(reads_per_len.Edge_num)] %>% select(starts_with("P3"))

lowess.list.full = list()

for(i in 1:nrow(reads_per_len.Edge_num.P1)){
  lowess.p1 = lowess(x = P1_times, y = reads_per_len.Edge_num.P1[i,], f = 6/36)
  lowess.p2 = lowess(x = P2_times, y = reads_per_len.Edge_num.P2[i,], f = 6/36)
  lowess.p3 = lowess(x = P3_times, y = reads_per_len.Edge_num.P3[i,], f = 6/36)
  lowess.list.full[[i]] = data.frame(
    plot = c(rep("P1", length(P1_times)), rep("P2", length(P2_times)), rep("P3", length(P3_times)) ),
    x = c(lowess.p1$x, lowess.p2$x, lowess.p3$x), 
    y = c(lowess.p1$y, lowess.p2$y, lowess.p3$y)
  )
}

lowess.full.df = bind_rows(lowess.list.full, .id = "OTU")
#lowess.list.full.df = lowess.list.full.df %>% mutate(row = row_number())
lowess.OTU_table = pivot_wider(lowess.full.df, id_cols = c(plot,x), names_from = OTU, values_from = y)

#There are negative values. need to add a constant
lowess.OTU_table[,3:ncol(lowess.OTU_table)] %>% range

lowess.adonis = adonis2(lowess.OTU_table[,3:ncol(lowess.OTU_table)]+0.04899737 ~ sin(2*pi*x/24) + cos(2*pi*x/24) + plot + x, data =lowess.OTU_table) #putting the diel signal variates in front because adonis is sensitive to order of input and these are quite small. Want to pick up max diel

lowess.adonis
lowess.adonis %>% str
lowess.adonis$R2

lowess.adonis.coefs = data.frame(cat = "Overall", plot = lowess.adonis$R2[3], trend = lowess.adonis$R2[4], cycle = lowess.adonis$R2[1]+lowess.adonis$R2[2])
#transform for plotting
lowess.adonis.coefs.long = lowess.adonis.coefs %>% pivot_longer(cols = -cat, names_to = "coefs", values_to = "var") 

###################################################
#Set up for lowess and aov on individual summaries
###################################################


#testing loop
#i = 1
#temp_tax = reads_per_len.Tax_name.top_50[i,1] %>% data.frame %>% unname
#temp_dat.P1 = temp_dat %>% select(starts_with("P1")) %>% as.data.frame %>% unname 
#temp_dat.P2 = temp_dat %>% select(starts_with("P2")) %>% as.data.frame %>% unname 
#temp_dat.P3 = temp_dat %>% select(starts_with("P3")) %>% as.data.frame %>% unname 

#plot(P1_times, temp_dat.P1)
#length(temp_dat.P1)

#foo = lowess(x = P1_times, y = temp_dat.P1, f = 6/36)
#lines(foo, col = "blue")
#foo = lowess(x = P1_times, y = temp_dat.P1, f = 3/36)
#lines(foo, col = "red")
#foo = lowess(x = P1_times, y = temp_dat.P1, f = 12/36)
#lines(foo, col = "black")

#lowess.p1 = lowess(x = P1_times, y = temp_dat.P1, f = 6/36)
#lowess.p2 = lowess(x = P2_times, y = temp_dat.P2, f = 6/36)
#lowess.p3 = lowess(x = P3_times, y = temp_dat.P3, f = 6/36)

#lowess_dat = data.frame(x = c(lowess.p1$x, lowess.p2$x, lowess.p3$x), y = c(lowess.p1$y, lowess.p2$y, lowess.p3$y))

#full_dat.lowess = cbind(mapped.metadata %>% select(Plot), lowess_dat)

#lowess.aov = aov(y ~ Plot + x + sin(2*pi*x/24) + cos(2*pi*x/24), data = full_dat.lowess)
#coefs.lowess = summary(lowess.aov)

#data.frame(
#  cat = temp_tax, 
#  plot = coefs.lowess[[1]][1,2], 
#  trend = coefs.lowess[[1]][2,2],
#  cycle = coefs.lowess[[1]][3,2] + coefs.lowess[[1]][4,2],
#  resid = coefs.lowess[[1]][5,2]
#  )


##################
#Tax_name summary
##################

coefs.list = list()

for(i in 1:50){
  temp_tax = reads_per_len.Tax_name.top_50[i,1] %>% data.frame %>% unlist %>% unname %>% as.character
  temp_dat = reads_per_len.Tax_name.top_50[i,2:ncol(reads_per_len.Tax_name.top_50)] 
  
  #filter by plot
  temp_dat.P1 = temp_dat %>% select(starts_with("P1")) %>% as.data.frame %>% unlist %>% unname
  temp_dat.P2 = temp_dat %>% select(starts_with("P2")) %>% as.data.frame %>% unlist%>% unname 
  temp_dat.P3 = temp_dat %>% select(starts_with("P3")) %>% as.data.frame %>% unlist%>% unname 
  #perform lowess smooth
  lowess.p1 = lowess(x = P1_times, y = temp_dat.P1, f = 6/36)
  lowess.p2 = lowess(x = P2_times, y = temp_dat.P2, f = 6/36)
  lowess.p3 = lowess(x = P3_times, y = temp_dat.P3, f = 6/36)
  #combine dat
  lowess_dat = data.frame(x = c(lowess.p1$x, lowess.p2$x, lowess.p3$x), y = c(lowess.p1$y, lowess.p2$y, lowess.p3$y))
  full_dat.lowess = cbind(mapped.metadata %>% select(Plot), lowess_dat)
  #run anova
  lowess.aov = aov(y ~ Plot + x + sin(2*pi*x/24) + cos(2*pi*x/24), data = full_dat.lowess)
  coefs.lowess = summary(lowess.aov)
  total.var = sum(coefs.lowess[[1]][,2])
  #add to list 
  coefs.list[[temp_tax]] = data.frame(
    #cat = temp_tax, 
    plot = coefs.lowess[[1]][1,2]/total.var, 
    trend = coefs.lowess[[1]][2,2]/total.var,
    cycle = (coefs.lowess[[1]][3,2] + coefs.lowess[[1]][4,2])/total.var,
    resid = coefs.lowess[[1]][5,2]/total.var,
    total.ssq = total.var,
    plot.p = ifelse(coefs.lowess[[1]][1,5] < 0.05, "sig", "n.s"),
    trend.p = ifelse(coefs.lowess[[1]][2,5] < 0.05, "sig", "n.s"),
    cycle.p = ifelse(coefs.lowess[[1]][3,5] < 0.05, "sig", "n.s")
  )
}

coefs.lowess.df = bind_rows(coefs.list, .id = "cat")

##################
#Edge_num summary
##################

coefs.list.en = list()

for(i in 1:50){
  temp_tax = paste(
    reads_per_len.Edge_num.top_50[i,2] %>% data.frame %>% unlist %>% unname %>% as.character,
    reads_per_len.Edge_num.top_50[i,1] %>% data.frame %>% unlist %>% unname %>% as.character,
    sep = "."
  )
  temp_dat = reads_per_len.Edge_num.top_50[i,3:ncol(reads_per_len.Edge_num.top_50)] 
  
  #filter by plot
  temp_dat.P1 = temp_dat %>% select(starts_with("P1")) %>% as.data.frame %>% unlist %>% unname
  temp_dat.P2 = temp_dat %>% select(starts_with("P2")) %>% as.data.frame %>% unlist%>% unname 
  temp_dat.P3 = temp_dat %>% select(starts_with("P3")) %>% as.data.frame %>% unlist%>% unname 
  #perform lowess smooth
  lowess.p1 = lowess(x = P1_times, y = temp_dat.P1, f = 6/36)
  lowess.p2 = lowess(x = P2_times, y = temp_dat.P2, f = 6/36)
  lowess.p3 = lowess(x = P3_times, y = temp_dat.P3, f = 6/36)
  #combine dat
  lowess_dat = data.frame(x = c(lowess.p1$x, lowess.p2$x, lowess.p3$x), y = c(lowess.p1$y, lowess.p2$y, lowess.p3$y))
  full_dat.lowess = cbind(mapped.metadata %>% select(Plot), lowess_dat)
  #run anova
  lowess.aov = aov(y ~ Plot + x + sin(2*pi*x/24) + cos(2*pi*x/24), data = full_dat.lowess)
  coefs.lowess = summary(lowess.aov)
  total.var = sum(coefs.lowess[[1]][,2])
  #add to list 
  coefs.list.en[[temp_tax]] = data.frame(
    #cat = temp_tax, 
    plot = coefs.lowess[[1]][1,2]/total.var, 
    trend = coefs.lowess[[1]][2,2]/total.var,
    cycle = (coefs.lowess[[1]][3,2] + coefs.lowess[[1]][4,2])/total.var,
    resid = coefs.lowess[[1]][5,2]/total.var,
    total.ssq = total.var,
    plot.p = ifelse(coefs.lowess[[1]][1,5] < 0.05, "sig", "n.s"),
    trend.p = ifelse(coefs.lowess[[1]][2,5] < 0.05, "sig", "n.s"),
    cycle.p = ifelse(coefs.lowess[[1]][3,5] < 0.05, "sig", "n.s")
  )
}

coefs.en.lowess.df = bind_rows(coefs.list.en, .id = "cat")
coefs.lowess.df = coefs.en.lowess.df

########################################
#Set up for plotting (cluster and table)
########################################

#perform clustering for variable order
coefs.lowess.df.ptc = coefs.lowess.df %>% select(plot, trend, cycle)
coefs.lowess.df.ptc.c = hclust(coefs.lowess.df.ptc %>% scale %>% dist, method = "ward.D")
#ordering vector for cat
clust_order = coefs.lowess.df.ptc.c$order
cat_order = coefs.lowess.df[clust_order,] %>% pull(cat)

#set up df for heatmap
coefs.lowess.df.ptc.p_vals = coefs.lowess.df %>% select(cat, plot.p, trend.p, cycle.p)
coefs.lowess.df.ptc.coefs = coefs.lowess.df %>% select(cat, plot, trend, cycle)

coefs_df = left_join(
    coefs.lowess.df.ptc.p_vals %>% pivot_longer(cols = c(-cat), names_to = "coef_sig", values_to = "sig"),
    coefs.lowess.df.ptc.coefs %>% pivot_longer(cols = c(-cat), names_to = "coef_var", values_to = "var"),
    by = c("cat") 
  )

#ordering vector for coefs
level_order = c("plot", "trend", "cycle")

####################
#Plot heatmap
####################

#Need to check ranges of var vals in order to *set* ranges in fill statement to make sure the top 50 and adonis match
#Use the same scale_fill_gradient2 call for both p1 and p4 (i.e., limits, colors, midpoint)
lowess.adonis.coefs.long$var %>% range
coefs_df$var %>% range

p1 = ggplot(coefs_df, 
       aes(
         x = factor(coef_var, level = level_order), 
         y = factor(cat, level = cat_order), 
         fill = var#, 
         #color = sig
         )
       ) +
  geom_tile() +
  #scale_color_manual(values = c("white", "black"), guide = "none") +
  scale_fill_gradient2(low = "white", mid = "#3690c0", high = "#034e7b", midpoint = 0.35, limits = c(0, 0.76)) +
  labs(x = "Proportion variance") +
  my_gg_theme +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2, 'cm'),
    #legend.justification='left',
    legend.text = element_text(size = 15),
    #axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p1
####################
#Plot RNA rel abd
####################


p2 = ggplot(top_50_OTUs,
       aes(
         x = factor(cat, level = cat_order),
         y = global_RA*100
       )
    ) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(y = "RNA rel. abd. (%)") +
  my_gg_theme +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )
p2
####################
#Plot cat colors
####################

#a table of higher level classifications of the cats, e.g., kingdom for phyla, KEGG B class. for KO, substrate for GH
cat_type = read.table("top_50_genera.phylum.txt", header = T, sep = "\t")

OTUs.cat_type = left_join(top_50_OTUs, cat_type, by = "Tax_name")
OTUs.cat_type$Phylum %>% as.factor %>% levels()

#order phylq
phylum_order = c(
  "Actinobacteria", 
  "Bacteroidetes", 
  "alpha-Proteobacteria",
  "beta-Proteobacteria",
  "gamma-Proteobacteria"#,
#  "Bacteria candidate phyla"
  )

p3 = ggplot(OTUs.cat_type,
       aes(
         x = factor(cat, level = cat_order),
         y = rep(1, length(cat_order)),
         fill = factor(Phylum, level = phylum_order)
       )
    ) +
  geom_raster() +
  coord_flip() +
  scale_fill_manual(values = cbPalette) +
  labs(y = "Phylum") +
  my_gg_theme +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    #axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 20)#,
    #legend.position = "bottom",
    #legend.justification=c(1,-1)
  )
p3 

##################
#Plot adonis
##################

p4 = ggplot(lowess.adonis.coefs.long,
            aes(
              x = factor(coefs, level = level_order),
              y = cat,
              fill = var
              )
            ) + 
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "#3690c0", high = "#034e7b", midpoint = 0.35, limits = c(0, 0.76),
                       guide = "none"
                       ) +
  #coord_flip() +
  my_gg_theme +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 18,face="bold"),
    axis.title.y = element_blank()
  )
p4

##################
#Set up the grid!
##################

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)
gp4 = ggplotGrob(p4)

nrow(gp1) #14
nrow(gp2) #12 # Need to add to rows to the *bottom* of this to match gp1
nrow(gp3) #12 # Need to add to rows to the *bottom* of this to match gp1
nrow(gp4) #12

ncol(gp1)
ncol(gp4)

#also need to set height of the figure grob from gp4 to reduce relative to gp1

gp1$heights
gp4$heights
gtable_show_layout(gp1)

gp4$heights[7] = unit(0.025, "null")
gp4$heights[12] = unit(1, "points")
#First joining 1 and 4 with 4 on top
gp14 = rbind(gp4, gp1)
plot(gp14)
gp14$heights
#check number rows for join to gp2 and gp3
nrow(gp14) #26 #after this join will need to add the difference in rows () (26-14 = 12) to the *top* of gp2 and gp3 to match the joined fig


#need to have the same number of rows per object
#pos = 0 to add to top, pos = -1 to add to bottom
gp2 = gtable_add_rows(gp2, heights = unit(c(1,1), "null"), pos = -1) 
gp3 = gtable_add_rows(gp3, heights = unit(c(1,1), "null"), pos = -1) 
#add 12 to the top
gp2 = gtable_add_rows(gp2, heights = unit(rep(0.01,12), "null"), pos = 0) 
gp3 = gtable_add_rows(gp3, heights = unit(rep(0.01,12), "null"), pos = 0) 

#gtable_add()
plot(gp2)
nrow(gp2)
nrow(gp3)

#Set panel sizes to make phylum coloring smaller
#or can use egg::set_panel_size to set all widths
gtable_show_layout(gp2)
gp2$widths
gp2$widths[5] = unit(0.5, "null")


gtable_show_layout(gp3)
gp3$widths
gp3$widths[5] = unit(0.2, "null")

nrow(gp3)
#gp1 = set_panel_size(g = gp1, width = unit(1, "null"))

gp12 = cbind(gp14, gp2)
plot(gp12)
gp123 = cbind(gp12, gp3)

plot(gp123)

pdf("summary_figs.20_gene.pdf", width = 16, height = 12)
grid::grid.draw(gp123)
dev.off()



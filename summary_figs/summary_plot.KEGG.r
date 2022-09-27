require(tidyverse)
require(gridExtra)
require(scales)
require(vegan)
require(lubridate)
require(gtable)
require(rlang)
require(RColorBrewer)

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
#KO terms
catToSum = "KO"

ref_len = read.table("read_counts.KO.phylodist.readCountRefLen.join.cleanStrings", header = T, sep = "\t") %>% as.tbl
ref_len$P2T31_metaT = NULL
ref_len$X = NULL

read_count.rarefied = readRDS(file = "intermediate_RDS/KO_rarefied_count.rds")

reads_per_len = data.frame(read_count.rarefied[,1:8], read_count.rarefied[,9:115]/ref_len[,9:115] ) %>% as.tbl
reads_per_len[is.na(reads_per_len)] = 0

reads_per_len.sumCat = reads_per_len %>% group_by(sumCat = !!parse_expr(catToSum)) %>% summarize_if(is.numeric,sum,na.rm = TRUE)

print("summarized count tbl:")
print(reads_per_len.sumCat)

reads_per_len = reads_per_len.sumCat
colnames(reads_per_len)[1] = "KO"

############################
#Checking counts on summary
############################

#how many tax names and how many edge nums
reads_per_len$KO %>% unique %>% length
#9706
reads_per_len %>% nrow
#9706

#Top 50 and relative abundance
reads_per_len.rowSum = data.frame(KO = reads_per_len$KO, sum = reads_per_len[2:length(colnames(reads_per_len))] %>% rowSums)
reads_per_len.rowSum[order(reads_per_len.rowSum$sum),]

top_50 = (reads_per_len.rowSum[order(reads_per_len.rowSum$sum, decreasing = T),])[1:50,]
top_50$global_RA = top_50$sum/sum(reads_per_len.rowSum$sum)

#Write top 50 KO terms
write.table(top_50$KO, "top_50_KO_terms.txt", quote = F, row.names = F, sep = "\t", col.names = F)

#Use top 50 list to filter table for NLS fitting
#var names in mapped meta Plot, hours.cumulative.RNA

plot(top_50$global_RA)

############################################
#Filter summarized tax_name table to top 50
############################################

reads_per_len.top_50 = reads_per_len %>% filter(KO %in% top_50$KO)

##########################################
#Overall PERMANOVA.
##########################################

P1_times = mapped.metadata %>% filter(Plot == "P1") %>% pull(hours.cumulative.RNA)
P2_times = mapped.metadata %>% filter(Plot == "P2") %>% pull(hours.cumulative.RNA)
P3_times = mapped.metadata %>% filter(Plot == "P3") %>% pull(hours.cumulative.RNA)

reads_per_len.P1 = reads_per_len[,2:ncol(reads_per_len)] %>% select(starts_with("P1"))
reads_per_len.P2 = reads_per_len[,2:ncol(reads_per_len)] %>% select(starts_with("P2"))
reads_per_len.P3 = reads_per_len[,2:ncol(reads_per_len)] %>% select(starts_with("P3"))

lowess.list.full = list()

for(i in 1:nrow(reads_per_len.P1)){
  lowess.p1 = lowess(x = P1_times, y = reads_per_len.P1[i,], f = 6/36)
  lowess.p2 = lowess(x = P2_times, y = reads_per_len.P2[i,], f = 6/36)
  lowess.p3 = lowess(x = P3_times, y = reads_per_len.P3[i,], f = 6/36)
  lowess.list.full[[i]] = data.frame(
    plot = c(rep("P1", length(P1_times)), rep("P2", length(P2_times)), rep("P3", length(P3_times)) ),
    x = c(lowess.p1$x, lowess.p2$x, lowess.p3$x), 
    y = c(lowess.p1$y, lowess.p2$y, lowess.p3$y)
  )
}

lowess.full.df = bind_rows(lowess.list.full, .id = "KO")
#lowess.list.full.df = lowess.list.full.df %>% mutate(row = row_number())
lowess.KO_table = pivot_wider(lowess.full.df, id_cols = c(plot,x), names_from = KO, values_from = y)

#If there are negative values. need to add a constant
lowess_range = lowess.KO_table[,3:ncol(lowess.KO_table)] %>% range

lowess.adonis = adonis2(lowess.KO_table[,3:ncol(lowess.KO_table)]+2.48622 ~ sin(2*pi*x/24) + cos(2*pi*x/24) + plot + x, data =lowess.KO_table) #putting the diel signal variates in front because adonis is sensitive to order of input and these are quite small. Want to pick up max diel

lowess.adonis
lowess.adonis %>% str
lowess.adonis$R2

lowess.adonis.coefs = data.frame(cat = "Overall", plot = lowess.adonis$R2[3], trend = lowess.adonis$R2[4], cycle = lowess.adonis$R2[1]+lowess.adonis$R2[2])
#transform for plotting
lowess.adonis.coefs.long = lowess.adonis.coefs %>% pivot_longer(cols = -cat, names_to = "coefs", values_to = "var") 

###################################################
#Set up for lowess and aov on individual summaries
###################################################

coefs.list = list()

for(i in 1:50){
  temp_ko = reads_per_len.top_50[i,1] %>% data.frame %>% unlist %>% unname %>% as.character
  temp_dat = reads_per_len.top_50[i,2:ncol(reads_per_len.top_50)] 
  
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
  coefs.list[[temp_ko]] = data.frame(
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

########################################
#Set up for plotting (cluster and table)
########################################

KO_gene_names = read.table("KO_term_names_rarefied_table.txt")
colnames(KO_gene_names) = c("cat", "name")

coefs.lowess.df = left_join(coefs.lowess.df, KO_gene_names, by = "cat")
#perform clustering for variable order
coefs.lowess.df.ptc = coefs.lowess.df %>% select(plot, trend, cycle)
coefs.lowess.df.ptc.c = hclust(coefs.lowess.df.ptc %>% scale %>% dist, method = "ward.D")
#ordering vector for cat
clust_order = coefs.lowess.df.ptc.c$order


cat_order = coefs.lowess.df[clust_order,] %>% pull(cat)
name_order = coefs.lowess.df[clust_order,] %>% pull(name)

#set up df for heatmap
coefs.lowess.df.ptc.p_vals = coefs.lowess.df %>% select(cat, name, plot.p, trend.p, cycle.p)
coefs.lowess.df.ptc.coefs = coefs.lowess.df %>% select(cat, name, plot, trend, cycle)

coefs_df = left_join(
    coefs.lowess.df.ptc.p_vals %>% pivot_longer(cols = c(-cat, -name), names_to = "coef_sig", values_to = "sig"),
    coefs.lowess.df.ptc.coefs %>% pivot_longer(cols = c(-cat, -name), names_to = "coef_var", values_to = "var"),
    by = c("cat", "name") 
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
         y = factor(name, level = name_order), 
         fill = var#, 
         #color = sig
         )
       ) +
  geom_tile() +
  #scale_color_manual(values = c("white", "black"), guide = "none") +
  scale_fill_gradient2(low = "white", mid = "#3690c0", high = "#034e7b", midpoint = 0.3, limits = c(0, 0.64)) +
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


p2 = ggplot(top_50,
       aes(
         x = factor(KO, level = cat_order),
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
cat_type = read.table("top_50_KEGG_map.txt", header = F, sep = "\t")
colnames(cat_type) = c("A", "B", "C", "KO")
cat_type$KO = paste("KO:", cat_type$KO, sep = "")
cat_type$KO %>% unique %>% length
cat_type$A %>% unique()


cat_type %>% filter(A %in% "A09150 Organismal Systems")
cat_type %>% filter(A %in% "A09190 Not Included in Pathway or Brite")
cat_type %>% filter(A %in% "A09180 Brite Hierarchies")
cat_type %>% filter(A %in% "A09120 Genetic Information Processing")

filter_out_A = c(
  "A09160 Human Diseases",
  #"A09180 Brite Hierarchies",
  #"A09190 Not Included in Pathway or Brite",
  "A09150 Organismal Systems",
  "A09130 Environmental Information Processing"
)

cat_type %>% filter(!A %in% filter_out_A) %>% pull(KO) %>% unique %>% length
cat_type %>% filter(!A %in% filter_out_A) %>% nrow
cat_type.filtered = cat_type %>% filter(!A %in% filter_out_A)
cat_type.filtered$A %>% unique()
cat_type.filtered$B %>% unique() %>% length
cat_type.filtered$C %>% unique() %>% length
cat_type.filtered$C %>% unique()

cat_type.filtered$A = gsub("A[0-9]+ ", "", KOs.cat_type$A, perl = T)
cat_type.filtered$B = gsub("B [0-9]+ ", "", KOs.cat_type$B, perl = T)

A_level_order = cat_type.filtered$A %>% unique()
#A_level_order = gsub("A[0-9]+ ", "", A_level_order, perl = T)
A_level_order
A_level_order = c(
  "Cellular Processes",
  "Metabolism",
  "Genetic Information Processing",
  "Brite Hierarchies",
  "Not Included in Pathway or Brite" 
)

#paste(cat_type.filtered$A, cat_type.filtered$B, sep = ".") %>% unique

B_level_order = c(
  cat_type.filtered %>% filter(A %in% A_level_order[1]) %>% pull(B) %>% unique,
  cat_type.filtered %>% filter(A %in% A_level_order[2]) %>% pull(B) %>% unique,
  cat_type.filtered %>% filter(A %in% A_level_order[3]) %>% pull(B) %>% unique,
  cat_type.filtered %>% filter(A %in% A_level_order[4]) %>% pull(B) %>% unique,
  cat_type.filtered %>% filter(A %in% A_level_order[5]) %>% pull(B) %>% unique
)

KOs.cat_type = left_join(top_50, cat_type.filtered, by = "KO")

#KOs.cat_type$A = gsub("A[0-9]+ ", "", KOs.cat_type$A, perl = T)
#KOs.cat_type$B = gsub("B [0-9]+ ", "", KOs.cat_type$B, perl = T)

p3 = ggplot(KOs.cat_type,
       aes(
         x = factor(KO, level = cat_order),
         y = factor(A, level = A_level_order),
         fill = factor(B, level = B_level_order)
       )
    ) +
  geom_raster() +
  coord_flip() +
  scale_fill_manual(values = c(brewer.pal(12, "Paired"), "grey", "black")) +
  scale_y_discrete(labels = scales::label_wrap(width = 15)) +
  labs(y = "Level A KEGG", fill = "Level B KEGG") +
  my_gg_theme +
  theme(
    legend.title = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    #axis.text.x = element_text(angle = 90, hjust = 1),
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
  scale_fill_gradient2(low = "white", mid = "#3690c0", high = "#034e7b", midpoint = 0.3, limits = c(0, 0.64),
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
gp3$widths[5] = unit(1.2, "null") #Setting this to larger than 0.2 because there are 5 panels

gp14$widths
gp14$widths[5] = unit(0.8, "null")
nrow(gp3)

gp12 = cbind(gp14, gp2)
plot(gp12)
gp123 = cbind(gp12, gp3)

plot(gp123)

pdf("summary_figs.KEGG.pdf", width = 31, height = 14)
grid::grid.draw(gp123)
dev.off()



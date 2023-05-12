rm(list=ls())

##plot packages##
library(tidyverse) ## V1.3.2
library(cowplot) ## V1.1.1
library(ggsci) ## V2.9
##analysis
library(microbiome) ## V1.18.0
library(microbiomeSeq) ## V0.1
library(phyloseq) ## V1.40.0
library(ANCOMBC) ## V1.64
library(factoextra) ## V1.07
library(vegan) ## V2.6-4

###ancom###
count_file = "count_file.txt"
tax_file = "tax_file.txt"
sample_file = "sample_file.txt"
pathway_file = "pathway_file.tsv"

formula = "group"
count_mat <- as.matrix(read.table(count_file,header=T,row.names = 1))
tax_mat <- as.matrix(read.table(tax_file,header = T,row.names = 1))
sample_mat <- read.table(sample_file,header=T,row.names = 1)

OTU = otu_table(count_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
SAM = sample_data(sample_mat)
SAM$group <- as.factor(SAM$group)
physeq = phyloseq(OTU, TAX, SAM)

sp_data = aggregate_taxa(physeq, "Species")
sp_out = ancombc(phyloseq = sp_data, formula = formula, 
                 p_adj_method = "BH", prv_cut = 0.1, lib_cut = 0, 
                 group = formula, struc_zero = F, neg_lb = F, tol = 1e-5, 
                 max_iter = 100, conserve = TRUE, alpha = 0.05, global = F)

sp_res = sp_out$res
tab <- bind_cols(
  row.names(sp_res$se),
  sp_res$lfc[,2],
  as.data.frame(sp_res$se[,2],row.names = F),
  as.data.frame(sp_res$p_val[,2],row.names = F),
  as.data.frame(sp_res$q_val[,2],row.names = F),
  sp_res$diff_abn[,2])
c_name = c("taxon_id","lfc","se","p","q","diff")
colnames(tab) <-  c_name
#tab %>% 
#  datatable(caption = "IN-N Result") %>%
#  formatRound(c_name[2:4],digits = 4)
tab$N_count <- 0
tab$IN_count <- 0
class(tab$se)
j = 1;
for(i in 1:nrow(count_mat)){
  if(rownames(count_mat)[i] %in% tab$taxon_id){
    tab[j,7] = sum(count_mat[i,1:13])
    tab[j,8] = sum(count_mat[i,14:26])
    j = j + 1 
  }
}
tab$count <- tab$IN_count+tab$N_count
sum_all <- sum(tab$count)
tab <-  tab %>%
  mutate(IN_rate = IN_count/count,
         N_rate = N_count/count,
         rate = (count*100)/sum_all)
tab$count = log2(tab$count)
tab$log_p <- -log10(tab$p)
tab <- tab %>%
  filter(se != 0)
tab$lfc <- -tab$lfc
tab$class <- "non"
tab <- tab %>%
  mutate(
    class = case_when(
      #lfc > 0 & diff == "TRUE" ~ "up",
      #lfc < 0 & diff == "TRUE" ~ "down",
      lfc > 1 & p<0.05 ~ "up",
      lfc < -1 & p<0.05 ~ "down",
      TRUE ~ "non"
    )
  )
write_csv(tab,"ancom.diff.csv")
tab <- read_csv("ancom.diff.csv")
ggplot(tab)+
  geom_point(aes(x=lfc,y=log_p,size=count,color=class))+
  labs(x="Log2FoldChange",y="-lg(P-value)",color = "Diffexpressed",size="Log2readscount")+
  scale_color_manual(
    values=c("#DC0000B2", "grey", "#00A087B2")
  )+
  scale_x_continuous(breaks = seq(-10,10,2))+
  scale_y_continuous(breaks = seq (0,20,5))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())

###alpha analysis##

ab_mat <- as.matrix(read.table("ab.species.level.txt",header=T,row.names = 1))

OTU = otu_table(ab_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
SAM = sample_data(sample_mat)
SAM$group <- as.factor(SAM$group)
physeq_ab = phyloseq(OTU, TAX, SAM)

sp_data_ab = aggregate_taxa(physeq_ab, "Species")

plot_anova_diversity(sp_data_ab,method = c("shannon"),grouping_column = "group",pValueCutoff=0.05)+
  theme_cowplot()+ 
  theme(axis.title.x=element_blank(), 
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), 
        text = element_text(size = 9))+
  scale_color_manual(
    values = c("#00B5E2","#CC0C00")
  )

###PCA analysis ###

data <- read.table ("ab.species.level.txt",row.names = 1,header = T)
data_t <- as.data.frame(t(data))
data_t <- data_t*100/rowSums(data_t)
data_t <- data_t %>% select_if(~max(.)>=0.01)
Group <- as.factor(as.character(ifelse(grepl("IN", rownames(data_t), perl=TRUE), "Infected", "Naive")))
Group <- factor(Group,levels = c("Naive","Infected"))
res.pca <- prcomp(data_t, scale = T)
fviz_pca_ind(res.pca, label="none",axes = c(1, 2), alpha.ind =1,
             habillage=Group,invisible="quali",pointsize = 2.5, addEllipses = TRUE,
             ellipse.level=0.66,palette = c("#00B5E2","#CC0C00"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.1)
  )
dist = vegdist(data_t, method = 'bray')
site = data.frame(sample = rownames(data_t),group = Group)
adonis_result_dis = adonis2(dist~group, site)

###pathway pca
data <- read_delim (pathway_file)

data2 <- data %>%
  filter(!grepl("UNM", data$`# Pathway`, perl=TRUE)&!grepl("\\|", data$`# Pathway`, perl=TRUE))%>%
  #filter(!grepl("\\|", data$`# Pathway`, perl=TRUE))%>%
  #filter(!grepl("\\|", data$`# Pathway`, perl=TRUE))%>%
  select(!`# Pathway`)
data_t <- as.data.frame(t(data2))
res.pca <- prcomp(data_t, scale = T)
Group <- as.factor(as.character(ifelse(grepl("IN", rownames(data_t), perl=TRUE), "Infected", "Naive")))
Group <- factor(Group,levels = c("Naive","Infected"))
fviz_pca_ind(res.pca, label="none",axes = c(1, 2), alpha.ind =1,
             habillage=Group,invisible="quali",pointsize = 2.5, addEllipses = TRUE,
             ellipse.level=0.66,palette = c("#00B5E2","#CC0C00"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.1)
  )
dist = vegdist(data_t, method = 'bray')
site = data.frame(sample = rownames(data_t),group = Group)
adonis_result_dis = adonis2(dist~group, site)

###bar plot###

data_file = "count.phylum.stat.txt"
data <- read.table(data_file,sep = "\t",header = T,row.names = 1)
data_t <- (t(data))
data_t <- data_t*100/rowSums(data_t)
data_p <- as.data.frame(t(data_t[1:2,]))
data_p$phylum = row.names(data_p)
data_l <- gather(data_p,group,rate,IN,N)
data_l$group <- (as.character(ifelse(grepl("IN", data_l$group, perl=TRUE), "Infected", "Naive")))
data_l$group <- factor(data_l$group,levels = c("Naive","Infected"))
mypal = pal_npg("nrc", alpha = 0.9)(9)
ggplot(data_l)+
  geom_bar(aes(x=group, fill=phylum, y= rate),stat='identity',position='stack')+
  scale_fill_manual(values=mypal)+
  ylab("Relative_abundance")+
  ylim(c(0,101)) + 
  theme_cowplot()+ theme(legend.key=element_blank())

###pathway test

data_path <- read_delim(pathway_file)
res <- data.frame(data_path$`# Pathway`)
colnames(res) <- c("pathway")
res$pvalue <- 0
for(i in 1:nrow(data)){
  datain <- data_path[,2:14]
  datan <- data_path[,15:27]
  t <- wilcox.test(as.numeric(datain[i,]),as.numeric(datan[i,]))
  res[i,2] <- as.numeric(t[3])
}
res$padjust <- p.adjust(res$pvalue,method ="BH")
res <- as.data.frame(res)
write_csv(res,"pathway_test.csv")
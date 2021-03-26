#cytof single clusterwise cell comparison
library(tidyverse)
library(ggpubr)
library(flowCore)
library(FSA)
library(broom)
library(caret)
library(RColorBrewer)

source(file.path("R", "functions.R"))

date <- Sys.Date()

if (!file.exists("data/panela.RData") | !file.exists("data/panelb.RData")) {
  panela <- readRDS("data/panela_2.rds")
  panelb <- readRDS("data/panelb_2.rds")
  
  #add diagnosis
  colnames(panela)[37] <- "Cluster"
  panela$Cluster <- as.factor(panela$Cluster)
  panela$Diagnosis <- ifelse(grepl("^HD", panela$sample), "Ctrl", 
                             ifelse(grepl("^IDH", panela$sample), "IDHmut", "IDHwt"))
  
  panela$Diagnosis <- factor(panela$Diagnosis, levels = c("Ctrl", "IDHmut", "IDHwt"))
  
  # clusters_below_promille <- as.data.frame(table(panela$Cluster, panela$Diagnosis))[,1][which(as.data.frame(table(panela$Cluster, panela$Diagnosis))[,3]< nrow(panela)/2000)]
  
  #  panela <- panela[!panela$Cluster %in% clusters_below_promille,]
  #  panela$Cluster <- droplevels(panela$Cluster)
  
  panela <- panela[, !colnames(panela) %in% "CD3"]
  levels(panela$Cluster) <- as.character(1:length(unique(panela$Cluster)))
  
  #panelb
  #add diagnosis
  colnames(panelb)[37] <- "Cluster"
  panelb$Cluster <- as.factor(panelb$Cluster)
  panelb$Diagnosis <- ifelse(grepl("^HD", panelb$sample), "Ctrl", 
                             ifelse(grepl("^IDH", panelb$sample), "IDHmut", "IDHwt"))
  
  panelb$Diagnosis <- factor(panelb$Diagnosis, levels = c("Ctrl", "IDHmut", "IDHwt"))
  
  #  clusters_below_promille <- as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,1][which(as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,3]< nrow(panelb)/2000)]
  
  #  panelb <- panelb[!panelb$Cluster %in% clusters_below_promille,]
  #  panelb$Cluster <- droplevels(panelb$Cluster)
  
  panelb <- panelb[, !colnames(panelb) %in% c("CD8a", "CD19")]
  
  levels(panelb$Cluster) <- as.character(1:length(unique(panelb$Cluster)))
  
  save(panela, file = "data/panela.RData")
  save(panelb, file = "data/panelb.RData")
} else {
  load("data/panela.RData")
  load("data/panelb.RData")
}

#plot clusters
panela %>% ggplot(aes(tSNE1, tSNE2, color=Cluster)) + geom_point()

#export data for cluster and diagnosis tsne
write_csv(panela[, c("tSNE1", "tSNE2", "Cluster", "Diagnosis")], "data/fig1g_data_panela_cluster_diagnosis_tsne.csv")

#export data for violin plot
write_csv(panela[, -1], "data/figS1d_data_panela_violin_plot.csv")

#cell number per sample
table(panela$Cluster, panela$sample) %>% colSums()


######
#### Plot viSNE with clusters & expression per cluster  
#######

#plot clusters - panel a
ggplot(panela, aes(tSNE1, tSNE2, color = Cluster)) +
  geom_point(pch = 19, size = .5) +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text=element_text(size=17),
        axis.text.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  scale_color_brewer(palette = "Set1") +
  guides(colour = guide_legend(override.aes = list(size=5)))


#plot diagnoses
ggplot(panela, aes(tSNE1, tSNE2, color = Diagnosis)) +
  geom_point(pch = 19, size = .5) +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text=element_text(size=17),
        axis.text.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  scale_color_brewer(palette = "Accent") +
  guides(colour = guide_legend(override.aes = list(size=5)))


#plot cell signatures
#signatures

#homeostatic
plot_expmap(gene=c("TREM2", "CX3CR1","TMEM119", "P2RY12"))
#log transformed
plot_expmap(gene=c("TREM2", "CX3CR1","TMEM119", "P2RY12"), logsc = T)

#apc
plot_expmap(gene=c("CD86","HLADR","CD74"))
#log transformed
plot_expmap(gene=c("CD86","HLADR","CD74"), logsc = T)

#violin plots
panela_red <- panela[,2:36]

panela_red <- panela_red %>%
  pivot_longer(CD45:CD61 ,"Protein", "Expression")

proteins <- hclust(dist(t(scale(panela[,2:35]))))
colnames(panela[,2:35])[proteins$order]

panela_red$Protein <- factor(panela_red$Protein, levels = colnames(panela[,2:35])[proteins$order])

ggplot(panela_red, aes(Protein, value, fill=Cluster)) +
  geom_violin(scale = "width", lwd=0.5) +
  facet_grid(facets = ~Cluster,
             drop = TRUE, 
             #space = "free", 
             #scales = "free_x", 
             #switch = "x",
             space = "free_x") +
  coord_flip() +
  theme_pubr() +
  labs(y="Expression (A.U.)",x=element_blank()) +
  scale_fill_brewer(palette = "Set1", guide=F) + 
  scale_y_continuous(limits=c(0,9), breaks = c(0,5)) +
  theme(
    strip.text.x = element_text(
      size = 12, color = "black", face="bold"
    ))

#marimekko chart
#plot marimekko plot
mosaicGG2(panela, "Cluster", "Diagnosis") +
  scale_fill_brewer(palette = "Accent")



#####  

#
panela <- panela

#plot clusterwise marker expression
terms <- colnames(panela[,-c(1,38)])
stats <- list()

#stat testing and plotting of individual genes            
for (i in terms) {
  panela2 <- panela[[i]]
  panela3 <- panela[,c("Cluster", "Diagnosis")] %>%
    bind_cols(data.frame("gene" = panela2))
  plt <- ggplot(na.omit(panela3), aes(Cluster, gene, fill = Diagnosis)) +
    geom_violin(scale = 'width', lwd=0.25, draw_quantiles = 0.5) +
    scale_fill_brewer(palette = "Accent") +
    theme_minimal() +
    labs(y='log2(normalized counts)', title=i) 
  print(plt)
  ggsave(paste0('plots/others/', date,"-",i ,'-cluster-Diagnosis.pdf'))
  
  list2 <- list()
  
  for (j in unique(panela$Cluster)) {
    mod <-dunnTest(gene ~ Diagnosis, data = panela3[panela3$Cluster == j,], method = "bh")  #tidy(kruskal.test(gene ~ Diagnosis, data = panela3[panela3$Cluster == i,]))
    print(j)
    print(mod)
  }
  
  
}


#plot ApoE and EMR1 signature genes
panela2 <- rowSums(panela[,c("ApoE","EMR1")])
panela3 <- panela[,c("Cluster", "Diagnosis")] %>%
  bind_cols(data.frame("gene" = panela2))

#list2 <- list()

df2 <- data.frame()

for (j in unique(panela$Cluster)) {
  mod <-dunnTest(gene ~ Diagnosis, data = panela3[panela3$Cluster == j,], method = "bh")  #tidy(kruskal.test(gene ~ Diagnosis, data = panela3[panela3$Cluster == i,]))
  mod2 <- as.data.frame(mod$res)
  mod2$Cluster <- j
  #list2[j] <- as.data.frame(mod$res)
  df2 <- bind_rows(df2, mod2)
}

clusters_below_promille <- as.data.frame(table(panela$Cluster, panela$Diagnosis))[,1][which(as.data.frame(table(panela$Cluster, panela$Diagnosis))[,3] > nrow(panela)/1000)]

df2 <- df2 %>% 
  mutate(padj = p.adjust(df2$P.unadj)) #%>%


#plot apc signature genes
panela2 <- rowSums(panela[,c("CD86","HLADR","CD74")])
panela3 <- panela[,c("Cluster", "Diagnosis")] %>%
  bind_cols(data.frame("gene" = panela2))

#list2 <- list()

df2 <- data.frame()

for (j in unique(panela$Cluster)) {
  mod <-dunnTest(gene ~ Diagnosis, data = panela3[panela3$Cluster == j,], method = "bh")  #tidy(kruskal.test(gene ~ Diagnosis, data = panela3[panela3$Cluster == i,]))
  mod2 <- as.data.frame(mod$res)
  mod2$Cluster <- j
  #list2[j] <- as.data.frame(mod$res)
  df2 <- bind_rows(df2, mod2)
}

clusters_below_promille <- as.data.frame(table(panela$Cluster, panela$Diagnosis))[,1][which(as.data.frame(table(panela$Cluster, panela$Diagnosis))[,3] > nrow(panela)/1000)]

df2 <- df2 %>% 
  mutate(padj = p.adjust(df2$P.unadj)) %>%
  dplyr::filter(padj < 0.01) #, !Cluster %in% clusters_below_promille

panela2 <- panela[panela$Cluster %in% clusters_below_promille,]

data_panela_apc <- na.omit(panela3)[panela$Cluster %in% unique(df2$Cluster),]
write_csv(data_panela_apc, "data/fig1h_right_data_panela_violin_apc.csv")

plt <- ggplot(data_panela_apc, aes(Cluster, gene, fill = Diagnosis)) +
  geom_violin(scale = 'width', lwd=0.25, draw_quantiles = 0.5) +
  scale_fill_brewer(palette = "Accent") +
  theme_minimal() +
  labs(y='log2(normalized counts)', title="APC signature") 
print(plt)

#homeostatic sig
panela2 <- rowSums(panela[,c("TREM2", "CX3CR1","TMEM119", "P2RY12" )])
panela3 <- panela[,c("Cluster", "Diagnosis")] %>%
  bind_cols(data.frame("gene" = panela2))


df2 <- data.frame()

for (j in unique(panela$Cluster)) {
  mod <-dunnTest(gene ~ Diagnosis, data = panela3[panela3$Cluster == j,], method = "bh")  #tidy(kruskal.test(gene ~ Diagnosis, data = panela3[panela3$Cluster == i,]))
  mod2 <- as.data.frame(mod$res)
  mod2$Cluster <- j
  #list2[j] <- as.data.frame(mod$res)
  df2 <- bind_rows(df2, mod2)
}

clusters_below_promille <- as.data.frame(table(panela$Cluster, panela$Diagnosis))[,1][which(as.data.frame(table(panela$Cluster, panela$Diagnosis))[,3]< nrow(panela)/1000)]

df2 <- df2 %>% 
  mutate(padj = p.adjust(df2$P.unadj)) %>%
  dplyr::filter(padj < 0.01) #, !Cluster %in% clusters_below_promille

write_csv(df2, "data/stat-testing_panela_homeo.csv")

data_panela_homeo <- na.omit(panela3)[panela$Cluster %in% unique(df2$Cluster),]
write_csv(data_panela_homeo, "data/fig1h_left_data_panela_violin_homeo.csv")

plt <- ggplot(data_panela_homeo, aes(Cluster, gene, fill = Diagnosis)) +
  geom_violin(scale = 'width', lwd=0.25, draw_quantiles = 0.5) +
  scale_fill_brewer(palette = "Accent") +
  theme_minimal() +
  labs(y='log2(normalized counts)', title="Homeostatic signature") 
print(plt)

#panelb

# Define the number of colors you want
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

clusters_below_promille <- as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,1][which(as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,3]< nrow(panelb)/1000)]

#plot clusters

tsne_plot <- ggplot(panelb, aes(tSNE1, tSNE2, color = Cluster)) +
  geom_point(pch = 19, size = .5) +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text=element_text(size=17),
        axis.text.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  scale_color_manual(values = mycolors) +
  guides(colour = guide_legend(override.aes = list(size=5)))

tsne_plot 


#plot diagnoses

tsne_plot <- ggplot(panelb, aes(tSNE1, tSNE2, color = Diagnosis)) +
  geom_point(pch = 19, size = .5) +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        text=element_text(size=17),
        axis.text.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  scale_color_brewer(palette = "Accent") +
  guides(colour = guide_legend(override.aes = list(size=5)))

tsne_plot 


#all genes violin plots
panelb_long <- panelb %>%
  pivot_longer(CD45:MRP14 ,"Protein", "Expression") 

proteins <- hclust(dist(t(scale(panelb[,2:34]))))
colnames(panelb[,colnames(panelb[,2:34])])[proteins$order]

panelb_long$Protein <- factor(panelb_long$Protein, levels =colnames(panelb[,colnames(panelb[,2:34])])[proteins$order])


ggplot(panelb_long, aes(Protein, value+.1, fill=Cluster)) +
  geom_violin(scale = "width", lwd=0.5) +
  facet_grid(facets = ~Cluster,
             drop = TRUE,
             space = "free_x") +
  coord_flip() +
  theme_pubr() +
  labs(y="Expression (A.U.)",x=element_blank()) +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(limits=c(0,15), breaks = c(0,5,10)) +
  theme(
    strip.text.x = element_text(
      size = 12, color = "black", face="bold"
    ))


## Apply cluster size cutoff cutoff

clusters_below_promille <- as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,1][which(as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,3]< nrow(panelb)/2000)]

panelb <- panelb[!panelb$Cluster %in% clusters_below_promille,]
panelb$Cluster <- droplevels(panelb$Cluster)


#export data for violin plot
write_csv(panelb[, -c(1, 29)], "data/figS1k_data_panelb_violin_plot.csv")

#export data for cluster and diagnosis tsne
write_csv(panelb[, c("tSNE1", "tSNE2", "Cluster")], "data/figS1e_data_panelb_cluster_tsne.csv")

write_csv(panelb[, c("tSNE1", "tSNE2", "Diagnosis")], "data/figS1f_data_panelb_diagnosis_tsne.csv")


#plot apc signature genes
panelb2 <- rowSums(as.data.frame(panelb[,c("HLADR")]))
panelb3 <- panelb[,c("Cluster", "Diagnosis")] %>%
  bind_cols(data.frame("gene" = panelb2))

#export data for tsne plot
write_csv(bind_cols(panelb3, panelb[, c( "tSNE1", "tSNE2")]), "data/figS1g_data_panelb_apc_tsne.csv")
#list2 <- list()

df2 <- data.frame()

for (j in unique(panelb$Cluster)) {
  mod <-dunnTest(gene ~ Diagnosis, data = panelb3[panelb3$Cluster == j,], method = "bh")  #tidy(kruskal.test(gene ~ Diagnosis, data = panelb3[panelb3$Cluster == i,]))
  mod2 <- as.data.frame(mod$res)
  mod2$Cluster <- j
  #list2[j] <- as.data.frame(mod$res)
  df2 <- bind_rows(df2, mod2)
}

clusters_below_promille <- as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,1][which(as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,3] < nrow(panelb)/1000)]

df2 <- df2 %>% 
  mutate(padj = p.adjust(df2$P.unadj)) %>%
  dplyr::filter(padj < 0.01, !Cluster %in% clusters_below_promille)

#cluster 3 shows up as one of the comparison clusters and it has extremely high hla-dr expression for a control - I examine the protein expression in this cluster
a <- panelb[panelb$Cluster  == 3,-1] %>%
  group_by(Diagnosis) %>%
  summarise_all(.funs="mean") %>%
  as.data.frame()
rownames(a) <- a$Diagnosis

pheatmap::pheatmap(t(a[,2:34,]), scale="row")
#the heatmap shows that the control cells of Cluster 3 have a high expression of CCR2 compared to IDHwt and IDHmut cells. 
#These cells appear to be monocytes. Therefore I exclude this cluster from the comparison. 

#df2 <- df2[df2$Cluster != "3"]
panelb2 <- panelb[panelb$Cluster %in% clusters_below_promille,] 

data_panelb_apc <- na.omit(panelb3)[panelb$Cluster %in% unique(df2$Cluster),]
write_csv(data_panelb_apc, "data/figS1j_data_panelb_violin_apc.csv")

plt <- ggplot(data_panelb_apc, aes(Cluster, gene, fill = Diagnosis)) +
  geom_violin(scale = 'width', lwd=0.25, draw_quantiles = 0.5) +
  scale_fill_brewer(palette = "Accent") +
  theme_minimal() +
  labs(y='log2(normalized counts)', title="APC signature") 
print(plt)
ggsave(paste0('plots/others/', date,"-" ,'-cluster-Diagnosis-APC-sig-panelb.pdf'))

df2_apc <- df2

write_csv(df2_apc, "data/stat-testing-apc-sig-panelb.csv")

#homeostatic sig
panelb2 <- rowSums(panelb[,c("CD115", "P2RY12" )])
panelb3 <- panelb[,c("Cluster", "Diagnosis")] %>%
  bind_cols(data.frame("gene" = panelb2))

write_csv(bind_cols(panelb3, panelb[, c("tSNE1", "tSNE2")]), "data/figS1i_data_panelb_homeo_tsne.csv")

#list2 <- list()

df2 <- data.frame()

for (j in unique(panelb$Cluster)) {
  mod <-dunnTest(gene ~ Diagnosis, data = panelb3[panelb3$Cluster == j,], method = "bh")  #tidy(kruskal.test(gene ~ Diagnosis, data = panelb3[panelb3$Cluster == i,]))
  mod2 <- as.data.frame(mod$res)
  mod2$Cluster <- j
  #list2[j] <- as.data.frame(mod$res)
  df2 <- bind_rows(df2, mod2)
}

clusters_below_promille <- as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,1][which(as.data.frame(table(panelb$Cluster, panelb$Diagnosis))[,3]< nrow(panelb)/1000)]

df2 <- df2 %>% 
  mutate(padj = p.adjust(df2$P.unadj)) %>%
  dplyr::filter(padj < 0.01)

data_panelb_homeo <- na.omit(panelb3)[panelb$Cluster %in% unique(df2$Cluster),]
write_csv(data_panelb_homeo, "data/figS1h_data_panelb_violin_homeo.csv")

plt <- ggplot(data_panelb_homeo, aes(Cluster, gene, fill = Diagnosis)) +
  geom_violin(scale = 'width', lwd=0.25, draw_quantiles = 0.5) +
  scale_fill_brewer(palette = "Accent") +
  theme_minimal() +
  labs(y='log2(normalized counts)', title="Homeostatic signature") 
print(plt)

df2_homeo <- df2

write_csv(df2_homeo, "data/stat-testing-homeostatic-sig-panelb.csv")

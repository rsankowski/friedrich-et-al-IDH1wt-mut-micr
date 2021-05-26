library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)
library(Seurat)
library(ggrepel)
library(MASS)
library(emmeans)

date = Sys.Date()
#load('data/sc.Robj')
load('data/seurat-integration-10x-wt-rh-gbm.RData')
source("bin/functions.R")

#export sessionInfo
writeLines(capture.output(sessionInfo()), "data/sessionInfo.txt")

if (!file.exists("data/metadata.RData")) {
  #build data frame
  metadata <- data.frame('Cluster' = all@active.ident, all@reductions$umap@cell.embeddings) 
  metadata$ID <- rownames(metadata)
  metadata$Timepoint <- factor(all$timepoint, levels = c('early', 'late'))
  metadata$Genotype <- factor(all$genotype, levels = c("WT", "RH"))
  metadata$Sample <- all$id
  
  if (!file.exists("data/retain_cl.RData")){
    
    retain_cl <- levels(all)
    
    save(retain_cl, file = "data/retain_cl.RData")
  } else {
    load("data/retain_cl.RData")
  }
  
  
  #assign cell type etc.
  load("data/df-umap.Robj")
  metadata <- metadata %>%
    left_join(df[, c("ID", "Cell_type", "Cell_type_simpl","Cell_family")])
  
  levels(metadata$Cell_type_simpl) <- c(levels(metadata$Cell_type_simpl), "Granulocytes", "Mast_cells", "Glial_cells")
  metadata$Cell_type_simpl[metadata$Cluster =="40"] <- "Granulocytes"
  metadata$Cell_type_simpl[metadata$Cluster =="41"] <- "Mo_der_Cells"
  metadata$Cell_type_simpl[metadata$Cluster =="42"] <- "Glial_cells"
  metadata$Cell_type_simpl[metadata$Cluster =="43"] <- "Monocytes"
  metadata$Cell_type_simpl[metadata$Cluster =="44"] <- "Mast_cells"
  metadata$Cell_type_simpl[metadata$Cluster =="45"] <- "Microglia"
  
  all$Cell_type <- metadata$Cell_type_simpl
  #save(all, file = 'data/seurat-integration-10x-wt-rh-gbm.RData')
  
  if (!file.exists('data/ord_clust.RData')) {
    ord_clust <- retain_cl
    save(ord_clust, file = 'data/ord_clust.RData')
  } else {
    load('data/ord_clust.RData')
  }
  
  #metadata <- metadata[metadata$Cluster %in% retain_cl,]
  
  save(metadata, file="data/metadata.RData")
} else {
  load("data/metadata.RData")
  load("data/retain_cl.RData")
  load('data/ord_clust.RData')
}

#export raw counts
raw_counts <- all@assays$RNA@counts
save(raw, file = "data/raw_counts.RData")

#export normalized counts
scaled_counts <- all[["SCT"]]@data
save(scaled_counts, file = "data/scaled_counts.RData")

#Plot timepoint umap
umap <- DimPlot(all, group.by = "timepoint") +
  scale_color_brewer(palette = "Accent") +
  theme_void()

umap

ggsave(paste0('plots/umap/', date, '-timepoint-umap-plot.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-timepoint-umap-plot.ps'), width = 8.57, height = 5.79)
umap
dev.off()

#timepoint umap no legends
umap <- DimPlot(all, group.by = "timepoint") +
  scale_color_brewer(palette = "Accent") +
  NoLegend() +
  theme_void()

umap

ggsave(paste0('plots/umap/', date, '-timepoint-umap-plot-no-legend.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-timepoint-umap-plot-no-legend.ps'), width = 8.57, height = 5.79)
umap
dev.off()

#marimekko plot - patients
mosaicGG2(metadata, "Cluster", "Timepoint", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = "Accent")
ggsave(paste0('plots/others/', date,'-compartment-marimekko-cluster-plot.pdf'))

#marimekko stat plot
mosaicGG(metadata, "Cluster", "Timepoint", rect_col = 'black', line_width = 0.1)
ggsave(paste0('plots/others/', date,'-timepoints-marimekko-cluster-stat-plot.pdf'))

#timepoints separated by facets
umap <- DimPlot(all, group.by = "timepoint") +
  scale_color_brewer(palette = "Accent") +
  theme_void() +
  facet_wrap(~ timepoint) +
  NoLegend()

umap

ggsave(paste0('plots/umap/', date, '-timepoint-umap-plot-separated-by-facets.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-timepoint-umap-plot-separated-by-facets.ps'), width = 8.57, height = 5.79)
umap
dev.off()

#plot late timepoint only
#late timepoint umap no legends
umap <- DimPlot(all[,colnames(all)[all$timepoint == "late"]], label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)) +
  theme_void() +
  NoLegend() 

umap

ggsave(paste0('plots/umap/', date, '-late-timepoint-umap-plot-no-legend.pdf'), width = 8.57, height = 5.79, useDingbats = F)  

#Plot genotype map
umap <- DimPlot(all, group.by = "genotype") +
  scale_color_brewer(palette = "Set1") +
  theme_void()

umap

ggsave(paste0('plots/umap/', date, '-Genotype-umap-plot.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-Genotype-umap-plot.ps'), width = 8.57, height = 5.79)
umap
dev.off()

#Genotype umap no legends
umap <- DimPlot(all, group.by = "genotype") +
  scale_color_brewer(palette = "Set1")  +
  theme_void() +
  NoLegend()

umap

ggsave(paste0('plots/umap/', date, '-Genotype-umap-plot-no-legend.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-Genotype-umap-plot-no-legend.ps'), width = 8.57, height = 5.79)
umap
dev.off()

#marimekko plot - patients
mosaicGG2(metadata, "Cluster", "Genotype", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = "Accent")
ggsave(paste0('plots/others/', date,'-compartment-marimekko-cluster-plot.pdf'))

#marimekko stat plot
mosaicGG(metadata, "Cluster", "Genotype", rect_col = 'black', line_width = 0.1)
ggsave(paste0('plots/others/', date,'-Genotypes-marimekko-cluster-stat-plot.pdf'))

#Genotypes separated by facets
umap <- DimPlot(all, group.by = "genotype") +
  scale_color_brewer(palette = "Set1") +
  theme_void() +
  facet_wrap(~ all$timepoint) +
  NoLegend()

umap

ggsave(paste0('plots/umap/', date, '-Genotype-umap-plot-separated-by-facets.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-Genotype-umap-plot-separated-by-facets.ps'), width = 8.57, height = 5.79)
umap
dev.off()


#Plot cluster umap map
umap <- DimPlot(all, label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)) +
  theme_void()
umap

ggsave(paste0('plots/umap/', date, '-clusters-umap-plot.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-clusters-umap-plot.ps'), width = 8.57, height = 5.79)
umap
dev.off()

#Plot cluster umap map
umap <- DimPlot(all, label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)) +
  theme_void() +
  NoLegend()

umap

ggsave(paste0('plots/umap/', date, '-clusters-umap-plot-no-legend.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-clusters-umap-plot-no-legend.ps'), width = 8.57, height = 5.79)
umap
dev.off()

#cluster late
#Plot cluster umap map
umap <- DimPlot(all[, colnames(all)[all$timepoint == "late"]], label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)) +
  theme_void() +
  NoLegend()

umap

ggsave(paste0('plots/umap/', date, '-late-clusters-umap-plot-no-legend.pdf'), width = 8.57, height = 5.79, useDingbats = F)  


#plot cell types
all$Cell_type <- metadata$Cell_type_simpl 
umap <- DimPlot(all, group.by = "Cell_type", label = T) +
  #scale_color_brewer(type = "qual", palette = "Set3") +
  scale_color_manual(values = colors_many) +
  theme_void()
umap

ggsave(paste0('plots/umap/', date, '-cell-types-umap-plot.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-cell-types-umap-plot.ps'), width = 8.57, height = 5.79)
umap
dev.off()

#plot cell types
umap <- DimPlot(all, group.by = "Cell_type", label = T) +
  #scale_color_brewer(type = "qual", palette = "Set3") +
  scale_color_manual(values = colors_many) +
  theme_void() +
  NoLegend()

umap

ggsave(paste0('plots/umap/', date, '-cell-types-umap-plot-no-legend.pdf'), width = 8.57, height = 5.79)  

postscript(paste0('plots/umap/', date, '-cell-types-umap-plot-no-legend.ps'), width = 8.57, height = 5.79)
umap
dev.off()



#marimekko myeloid celltypes per genotype
#marimekko plot
a <-metadata[metadata$Cell_type_simpl %in% c("Microglia", "Monocytes", "Mo_der_Cells", "DCs", "Granulocytes") & metadata$Timepoint == "early",] %>%
  mosaicGG2("Cell_type_simpl", "Genotype", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = "Set1") 
a

ggsave(paste0('plots/others/', date,'-genotype-celltype-marimekko-early.pdf'))

a <-metadata[metadata$Cell_type_simpl %in% c("Microglia", "Monocytes", "Mo_der_Cells", "DCs", "Granulocytes") & metadata$Timepoint == "late",] %>%
  mosaicGG2("Cell_type_simpl", "Genotype", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = "Set1") 
a

ggsave(paste0('plots/others/', date,'-genotype-celltype-marimekko-late.pdf'))



#plot cell signatures
signature_genes <-  read_excel('/home/roman/Documents/Single cell analysis/EAE_Final/Cluster-information-dimRed.xlsx', 'Core signature',skip = 2)

#take out calm1
signature_genes$T_cells[4] <- NA
#signature_genes <- data.frame("ahr"=c("Ahr", "Ahrr", "Cyp1a1", "Cyp1b1", "Entpd1", "Hif1a", "Il1b", "Il27", "Klf4", "Pparg", "Stat1", "Stat3", "Tiparp", "Vegfa" ), stringsAsFactors = F)


for (i in colnames(signature_genes)) {
  tryCatch({
    pdf(paste0('plots/umap/', date, as.character(i), '-gene_signature.pdf'), width = 8.57, height = 5.79, useDingbats = F)
    pl <- plot_expmap_seurat(features=c(na.omit(signature_genes[[i]])), point_size = .05)
    print(pl)
    dev.off()
    
    pdf(paste0('plots/umap/', date, as.character(i), '-gene_signature-logsc.pdf'), width = 8.57, height = 5.79, useDingbats = F)
    pl <- plot_expmap_seurat(features=c(na.omit(signature_genes[[i]])), point_size = .05, logsc = T)
    print(pl)
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
}     

for (i in colnames(signature_genes)) {
  tryCatch({
    pdf(paste0('plots/umap/', date, as.character(i), '-late-gene_signature.pdf'), width = 8.57, height = 5.79, useDingbats = F)
    pl <- plot_expmap_seurat(object=all[, colnames(all)[all$timepoint == "late"]] , features=c(na.omit(signature_genes[[i]])), point_size = .05)
    print(pl)
    dev.off()
    
    pdf(paste0('plots/umap/', date, as.character(i), '-late-gene_signature-logsc.pdf'), width = 8.57, height = 5.79, useDingbats = F)
    pl <- plot_expmap_seurat(object=all[, colnames(all)[all$timepoint == "late"]] , features=c(na.omit(signature_genes[[i]])), point_size = .05, logsc = T)
    print(pl)
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  on.exit(dev.off())
}     

#find cluster markers
all.markers <- FindAllMarkers(all, only.pos = T)

save(all.markers, file= "data/all_markers.RData")
write.csv(all.markers, "data/all_markers.csv")

#plot single cell gene expression
#up_genes <- load_data(file.path('data/Cluster specific genes/Up'))


load("data/cluster-markers.RData") 
plot_genes <- lapply(markers_map, function(x) subset(x, padj<0.05, fc>2)) %>%
  lapply(function(x) x$GENEID = rownames(x)) %>%
  unlist %>%
  unique

#plot_genes <-  as.character(unique(up_genes$GENEID))

#plot_genes <- stringr::str_to_title(read_csv('data/genes-mirco.csv')[['Gene']])

for (i in plot_genes) {
  tryCatch({
    #postscript(paste0('plots/umap/', i, '.ps'), width = 8.57, height = 5.79)
    pdf(paste0('plots/umap/', i, '.pdf'), width = 8.57, height = 5.79)
    pl <- plot_expmap_seurat(features=c(i), point_size = 1.5, line_width = 0)
    print(pl)
    dev.off()
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #on.exit(dev.off())
}     

#coparison between early and late samples
df2 <- metadata %>%
  group_by(Genotype, Sample, Timepoint, Cell_type_simpl) %>%
  summarise(Cell_count = n()) %>%
  group_by(Genotype, Timepoint, .drop = F) %>%
  ungroup() %>%
  group_by(Sample, .drop=F) %>%
  mutate(rel_count = Cell_count/sum(Cell_count))

df3 <- df2 %>%
  group_by(Genotype, Timepoint, Cell_type_simpl) %>%
  summarise(mean_rel_count = mean(rel_count),
            sd_rel_count = sd(rel_count),
            se_rel_count = sd(rel_count)/sqrt(n()-1))

#plot myeloid cells
point_line <- df2[df2$Cell_type_simpl %in% c("Microglia", "Monocytes", "Mo_der_Cells", "DCs", "Granulocytes"),] %>%
  ggplot(aes(x=Timepoint, y=rel_count * 100, color = Cell_type_simpl, group = Genotype)) +
  geom_point(size=5, alpha=.5) +
  facet_grid(Genotype~Cell_type_simpl) +
  #geom_line(data = df3, aes(x=Timepoint, y=mean_rel_count), group=df3$Genotype) +
  #geom_errorbar(data = df3, aes(ymin=mean_rel_count-se_rel_count, ymax=mean_rel_count+se_rel_count), width=0.1) +
  stat_summary(fun=mean, geom='line') +
  stat_summary(fun.data=mean_se, geom='errorbar', width=.5) +
  theme_bw() +
  scale_y_log10() +
  #scale_color_brewer(type = "qual", palette = "Set1") 
  scale_color_manual(values = colors_many[c(1,2,3,4,7)]) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x="Timepoint", y="% of all cells")

point_line

ggsave("plots/others/point_line_plots_celltypes_per_timepoint_myeloid_cells.pdf")

postscript("plots/others/point_line_plots_celltypes_per_timepoint_myeloid_cells.ps")

df2[df2$Cell_type_simpl %in% c("Microglia", "Monocytes", "Mo_der_Cells", "DCs", "Granulocytes"),] %>%
  ggplot(aes(x=Timepoint, y=rel_count, color = Cell_type_simpl, group = Genotype)) +
  geom_point(size=5) +
  facet_grid(Genotype~Cell_type_simpl) +
  #geom_line(data = df3, aes(x=Timepoint, y=mean_rel_count), group=df3$Genotype) +
  #geom_errorbar(data = df3, aes(ymin=mean_rel_count-se_rel_count, ymax=mean_rel_count+se_rel_count), width=0.1) +
  stat_summary(fun=mean, geom='line') +
  stat_summary(fun.data=mean_se, geom='errorbar', width=.5) +
  theme_bw() +
  scale_y_log10() +
  #scale_color_brewer(type = "qual", palette = "Set1") 
  scale_color_manual(values = colors_many) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x="Timepoint", y="% of all cells")

dev.off()


#plot all cells
point_line <- df2 %>%
  ggplot(aes(x=Timepoint, y=rel_count * 100, color = Cell_type_simpl, group = Genotype)) +
  geom_point(size=5, alpha=.5) +
  facet_grid(Genotype~Cell_type_simpl) +
  #geom_line(data = df3, aes(x=Timepoint, y=mean_rel_count), group=df3$Genotype) +
  #geom_errorbar(data = df3, aes(ymin=mean_rel_count-se_rel_count, ymax=mean_rel_count+se_rel_count), width=0.1) +
  stat_summary(fun=mean, geom='line') +
  stat_summary(fun.data=mean_se, geom='errorbar', width=.5) +
  theme_bw() +
  scale_y_log10() +
  #scale_color_brewer(type = "qual", palette = "Set1") 
  scale_color_manual(values = colors_many[c(1,2,3,4,7)]) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x="Timepoint", y="% of all cells")

point_line

ggsave("plots/others/point_line_plots_celltypes_per_timepoint_all_cells.pdf")

postscript("plots/others/point_line_plots_celltypes_per_timepoint_all_cells.ps")

df2 %>%
  ggplot(aes(x=Timepoint, y=rel_count, color = Cell_type_simpl, group = Genotype)) +
  geom_point(size=5) +
  facet_grid(Genotype~Cell_type_simpl) +
  #geom_line(data = df3, aes(x=Timepoint, y=mean_rel_count), group=df3$Genotype) +
  #geom_errorbar(data = df3, aes(ymin=mean_rel_count-se_rel_count, ymax=mean_rel_count+se_rel_count), width=0.1) +
  stat_summary(fun=mean, geom='line') +
  stat_summary(fun.data=mean_se, geom='errorbar', width=.5) +
  theme_bw() +
  scale_y_log10() +
  #scale_color_brewer(type = "qual", palette = "Set1") 
  scale_color_manual(values = colors_many) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x="Timepoint", y="% of all cells")

dev.off()

mod <- glm.nb(Cell_count ~ Timepoint * Cell_type_simpl , data=df2[df2$Genotype=="WT",])
emms <- emmeans(mod,
                ~ Timepoint : Cell_type_simpl)
comps <- pairs(emms, adjust="tukey") %>%
  as.data.frame() %>%
  filter(p.value < .05)

write.csv(comps, "data/wt-significant_stat_comparison_cell_types.csv")

mod <- glm.nb(Cell_count ~ Timepoint * Cell_type_simpl , data=df2[df2$Genotype=="RH",])
emms <- emmeans(mod,
                ~ Timepoint : Cell_type_simpl)
comps <- pairs(emms, adjust="tukey") %>%
  as.data.frame() %>%
  filter(p.value < .05)

write.csv(comps, "data/rh-significant_stat_comparison_cell_types.csv")


#donut plots
df3[df3$Cell_type_simpl %in% c("Microglia", "Monocytes", "Mo_der_Cells", "DCs") & df3$Timepoint == "early",] %>% #, "Granulocytes"
  ggplot(aes(x=2, y=mean_rel_count, fill=Cell_type_simpl)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_manual(values=colors_many[c(1,2,3,4,7)]) +
  facet_wrap(~Genotype) +
  xlim(0.5, 2.5)

ggsave("plots/others/early_donutplot_cell_types.pdf")

df3[df3$Cell_type_simpl %in% c("Microglia", "Monocytes", "Mo_der_Cells", "DCs") & df3$Timepoint == "late",] %>% #, "Granulocytes"
  ggplot(aes(x=2, y=mean_rel_count, fill=Cell_type_simpl)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_manual(values=colors_many[c(1,2,3,4,7)]) +
  facet_wrap(~Genotype) +
  xlim(0.5, 2.5)

ggsave("plots/others/late_donutplot_cell_types.pdf")

#stat testing
metadata[metadata$Cell_type_simpl %in% c("Microglia", "Monocytes", "Mo_der_Cells", "DCs", "Granulocytes") & metadata$Timepoint == "late",] %>%
  hyper_test_n(data=., var1 = "Genotype",var2="Cell_type_simpl")

#make volcano plots
all$Cell_geno <- factor(paste(all$timepoint , all$Cell_type, all$genotype, sep = "_"), levels = c(paste(c("early"), c("Microglia_WT", "Microglia_RH", "Mo_der_Cells_WT", "Mo_der_Cells_RH", "Monocytes_WT", "Monocytes_RH", "DCs_WT", "DCs_RH", "Granulocytes_WT", "Granulocytes_RH","T_cells_WT", "T_cells_RH", "B_cells_WT", "B_cells_RH", "Glial_cells_WT", "Glial_cells_RH", "Mast_cells_WT", "Mast_cells_RH"), sep="_"),
                                                                                                  paste(c("late"), c("Microglia_WT", "Microglia_RH", "Mo_der_Cells_WT", "Mo_der_Cells_RH", "Monocytes_WT", "Monocytes_RH", "DCs_WT", "DCs_RH", "Granulocytes_WT", "Granulocytes_RH","T_cells_WT", "T_cells_RH", "B_cells_WT", "B_cells_RH", "Glial_cells_WT", "Glial_cells_RH", "Mast_cells_WT", "Mast_cells_RH"), sep="_")))
Idents(all) <- all$Cell_geno

micr <- FindMarkers(all, ident.1 = "early_Microglia_WT",
                    ident.2 = "early_Microglia_RH",
                    logfc.threshold = 0.01,
                    min.pct = 0.01) %>%
  rownames_to_column(var="gene")  %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))
#save(micr, file = "data/micr_early_diffgenes_wt_vs_rh.RData")
load("data/micr_early_diffgenes_wt_vs_rh.RData")

top5_wt <- micr %>%
  top_n(5, avg_log2FC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes)

top_5_both <- micr %>%
  top_n(-5, avg_log2FC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes) %>%
  bind_rows(top5_wt) 

micr <- micr %>%
  left_join(top_5_both) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_log2FC) > .25, "sig.", "not sig."))


micr_volcano <- ggplot(micr, aes(x=avg_logFC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  geom_point(size=5) + #aes(size=avg_logFC)
  geom_text_repel(size=7, box.padding=1.15) +
  expand_limits(y=c(0, 350), x=c(-2.25, 2.25)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. logFC", y="-log10 transf. adj. p-value")

micr_volcano 

postscript("plots/others/volcano_plot_early_micr.ps", width = 7, height = 16)
micr_volcano
dev.off()

#late micr
late_micr <- FindMarkers(all, ident.1 = "late_Microglia_WT",
                         ident.2 = "late_Microglia_RH",
                         logfc.threshold = 0.01,
                         min.pct = 0.01) %>%
  rownames_to_column(var="gene")  %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))
save(late_micr, file = "data/micr_late_diffgenes_wt_vs_rh.RData")

top5_wt <- late_micr %>%
  top_n(10, avg_logFC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes)

top_5_both <- late_micr %>%
  top_n(-10, avg_logFC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes) %>%
  bind_rows(top5_wt) 

late_micr <- late_micr %>%
  left_join(top_5_both) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_logFC) > .25, "sig.", "not sig."))


micr_volcano <- ggplot(late_micr, aes(x=avg_logFC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  geom_point(size=5) + #aes(size=avg_logFC)
  geom_text_repel(size=7, box.padding=1.15) +
  expand_limits( x=c(-.55, .55)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. logFC", y="-log10 transf. adj. p-value")

micr_volcano 

pdf("plots/others/volcano_plot_late_micr.pdf", width = 10, height = 10, useDingbats = F)
micr_volcano
dev.off()

#
mfs <- FindMarkers(all, 
                   ident.1 = "late_Mo_der_Cells_WT",
                   ident.2 = "late_Mo_der_Cells_RH",
                   logfc.threshold = 0.01,
                   min.pct = 0.01)  %>%
  rownames_to_column(var="gene")  %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))
#save(mfs, file = "data/monmacs_late_diffgenes_wt_vs_rh.RData")

top5_wt <- mfs %>%
  top_n(10, avg_logFC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes)

top_5_both <- mfs %>%
  top_n(-10, avg_logFC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes) %>%
  bind_rows(top5_wt) 

mfs <- mfs %>%
  left_join(top_5_both) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_logFC) > .25, "sig.", "not sig."))



ggplot(mfs, aes(x=avg_logFC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  geom_point(size=5) + #aes(size=avg_logFC)
  geom_text_repel(size=7, box.padding=1.15) +
  expand_limits( x=c(-1.1, 1.1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. logFC", y="-log10 transf. adj. p-value")

ggsave("plots/others/volcano_plot_late_macs.pdf", width = 10, height = 10, useDingbats = F)

#dcs
dcs <- FindMarkers(all, 
                   ident.1 = "late_DCs_WT",
                   ident.2 = "late_DCs_RH",
                   logfc.threshold = 0.01,
                   min.pct = 0.01)  %>%
  rownames_to_column(var="gene")  %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))
save(dcs, file = "data/modcs_late_diffgenes_wt_vs_rh.RData")

top5_wt <- dcs %>%
  top_n(10, avg_log2FC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes)

top_5_both <- dcs %>%
  top_n(-10, avg_log2FC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes) %>%
  bind_rows(top5_wt) 

dcs <- dcs %>%
  left_join(top_5_both) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_log2FC) > .25, "sig.", "not sig."))

ggplot(dcs, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  geom_point(size=5) + #aes(size=avg_logFC)
  geom_text_repel(size=7, box.padding=1.15) +
  expand_limits( x=c(-1.1, 1.1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. logFC", y="-log10 transf. adj. p-value")

ggsave("plots/others/volcano_plot_late_dcs.pdf", width = 10, height = 10, useDingbats = F)

#modcs 
#reassign modc cluster
levels(all$Cell_geno) <- c(levels(all$Cell_geno), "early_DCs_RH_modc", "early_DCs_WT_modc", "late_DCs_RH_modc", "late_DCs_WT_modc")
all$Cell_geno[all$seurat_clusters %in% c(17, 26)] <- paste(all$Cell_geno[all$seurat_clusters %in% c(17, 26)], "modc", sep="_")
Idents(all) <- all$Cell_geno

modc <- FindMarkers(all, ident.1 = "late_DCs_WT_modc",
                    ident.2 = "late_DCs_RH_modc",
                    logfc.threshold = 0.01,
                    min.pct = 0.01) %>%
  rownames_to_column(var="gene")  %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))
save(modc, file = "data/modc_late_diffgenes_wt_vs_rh.RData")
load("data/modc_late_diffgenes_wt_vs_rh.RData")

top5_wt <- modc %>%
  top_n(0, avg_log2FC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes)

top_5_both <- modc %>%
  top_n(-12, avg_log2FC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes) %>%
  bind_rows(top5_wt) 

modc <- modc %>%
  left_join(top_5_both) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_log2FC) > .25, "sig.", "not sig."))


micr_volcano <- ggplot(modc, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  geom_point(size=5) + #aes(size=avg_logFC)
  geom_text_repel(size=7, box.padding=1.15) +
  expand_limits(x=c(-.75, .75)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. logFC", y="-log10 transf. adj. p-value")

micr_volcano 

postscript("plots/others/volcano_plot_late_modc.ps", width = 7, height = 16)
micr_volcano
dev.off()

#heatmap of the genes in the trajectory analysis
all2 <- all[, colnames(all)[all$timepoint == "early" & all$seurat_clusters %in% c("1","14", "13", "6", "2", "7", "9")]]
levels(all2) <- c("1","14", "13", "6", "2", "7", "9")

markers_early_micr <- FindAllMarkers(all2)          

save(markers_early_micr, file = "data/markers_early_micr.RData")
write_csv(markers_early_micr, "data/markers_early_micr.csv")

#all.markers <- all.markers[!grepl("^(Htra|Lin|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", all.markers$gene),]
#all.markers$cluster <- factor(all.markers$cluster, levels = ord_clust)

top20 <- markers_early_micr %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC) #filter(!cluster %in% c("5","8")) %>% 

heat <- DoHeatmap(all2,features = top20$gene, group.colors = rev(colors_pat)) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(paste0("plots/heatmaps/",date,"-early-micr-top15-gene-heatmap.pdf"), width = 30, height = 20)

DoHeatmap(all2,features = top20$gene, label = F, group.bar = F, size = 7)  + 
  scale_fill_viridis(option = "B") + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text = element_blank()
  ) +
  NoLegend()

ggsave("plots/heatmaps/heatmap_early.png")          

#plot expression of slc3a2 and slc7a5
df <- data.frame(Cell_type=all$Cell_type, Timepoint = all$timepoint, Expr = colSums(as.matrix(all[["SCT"]]@counts[c("Slc3a2", "Slc7a5"),])), Genotype=all$genotype)

#violin plot of the gene expression
df %>%
  ggplot(aes(x=Cell_type, Expr, fill=Cell_type)) +
  geom_boxplot() +
  facet_grid(Genotype~Timepoint) +
  coord_flip() +
  theme_bw()

#plot myeloid cells
point_line <- df %>%
  ggplot(aes(x=Timepoint, y=Expr, color = Cell_type, group = Genotype)) +
  facet_grid(Genotype~Cell_type) +
  stat_summary(fun=mean, geom='point') +
  stat_summary(fun=mean, geom='line') +
  stat_summary(fun.data=mean_se, geom='errorbar', width=.5) +
  theme_bw() +
  #scale_y_log10() +
  #scale_color_brewer(type = "qual", palette = "Set1") 
  scale_color_manual(values = colors_many) + #[c(1,2,3,4,7)]
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x="Timepoint", y="Expression")

point_line

#ggsave("plots/others/point_line_plots_celltypes_per_timepoint_myeloid_cells.pdf")



library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)
library(Seurat)
library(MASS)
library(rebus)
library(emmeans)

date = Sys.Date()
source("R/functions.R")
load("data/seurat-integration-wt-rh-gbm-final.RData")

#get metadata 
    metadata <- all@meta.data

# Figure 1 b left - summary statistics
    fig1b <- metadata[!duplicated(metadata[,c("Age", "ID", "Diagnosis")]),]
    fig1b %>%
      group_by(Diagnosis) %>%
      summarise(Mean_age = mean(Age),
                sd = sd(Age))
    table(fig1b$Gender, fig1b$Diagnosis)

# Figure 1c left - plot cell types (as indicated by circles)
    DimPlot(all, reduction = "umap", group.by = "Celltype", pt.size = 1, label = T) +
      scale_color_manual(values = colors_many[c(1, 2,6, 8,7)]) +
      theme_void() +
      NoLegend()

# Figure 1c left -  plot clusters
    DimPlot(all, reduction = "umap", pt.size = 1, label = T) +
      scale_color_manual(values = rev(colors_pat)) +
      theme_void() +
      NoLegend()

# Figure 1c right top - microglia signature
    signature_genes <- data.frame("monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                                  "macrophages"=c("MRC1", "MS4A7", "CD163", "LYVE1", "STAB1"),
                                  "microglia"= c('P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5'),
                                  "tcell"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                                  "myeloid"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                                  "oligodendrocyte"=c('MBP',  'MOG', 'MAG', 'PLP1', NA),
                                  "bcells"=c('CD79A', 'IGHG4', 'IGLL5', NA, NA),
                                  "astrocyte"=c("GFAP", "HEPACAM","SOX9","AQP4",NA),
                                  "apc"=c("CD74", "CD80", "CD86", "HLA-DRA", "CD40"), stringsAsFactors = F)
    
    plot_expmap_seurat(na.omit(signature_genes[["microglia"]]), .retain_cl = levels(all)) + labs(subtitle = "Microglia signature")

# Figure 1c right bottom - monocyte signature
    plot_expmap_seurat(na.omit(signature_genes[["monocytes"]]), .retain_cl = levels(all)) + labs(subtitle = "Monocyte signature")

# Figure 1d left - diagnosis umap
    DimPlot(all, reduction = "umap", group.by = "Diagnosis", pt.size = 1, label = T) +
      scale_color_brewer(palette = "Accent") +
      theme_void() 

# Figure 1d right - marimekko chart of diagnoses
    mosaicGG2(metadata, "seurat_clusters", "Diagnosis", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
      scale_fill_brewer(palette = "Accent")
    
    #statistical testing of diagnoses per cluster
    hyper_test_n(data= metadata, var1 = "integrated_snn_res.0.6", var2 = "Diagnosis") 

# Figure 1e - top15 gene heatmap
    if (!file.exists(file.path("data","cluster_marker_genes_human.RData"))) {
      all.markers <- FindAllMarkers(all, only.pos = T)
      save(all.markers, file = file.path("data","cluster_marker_genes_human.RData"))
    } else {
      load(file.path("data","cluster_marker_genes_human.RData"))
    }
    
    #remove uninformative/non-microglia/stress-associated genes
    all.markers <- all.markers[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", all.markers$gene),]
    
    #extract top 15 genes per cluster
    top15 <- all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
    
    #plot
    heat <- DoHeatmap(all,features = top15$gene, group.colors = rev(colors_pat)) 
    heat + scale_fill_viridis(option = "B")

# Figure 1f - volcano plot
    all2 <- all[, colnames(all)[all$Celltype == "Micr" & all$Diagnosis != "Ctrl"]]
    Idents(all2) <- all2$Diagnosis
    
    micr <- FindMarkers(all2, ident.1 = "IDHwt_GBM",
                        ident.2 = "IDHmut_GBM",
                        logfc.threshold = 0.01,
                        min.pct = 0.01) %>%
      rownames_to_column(var="gene")  %>%
      mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))
    save(micr, file = "data/micr_diffgenes_wt_vs_mut.RData")
    load("data/micr_diffgenes_wt_vs_mut.RData")
    
    #remove uninformative/non-microglia/stress-associated genes
    micr <- micr[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP|EGR|LYVE1|PARP)", micr$gene),]
    
    top5_wt <- micr %>%
      top_n(15, avg_logFC) %>%
      mutate(show_genes = gene) %>%
      dplyr::select(gene, show_genes)
    
    top_5_both <- micr %>%
      top_n(-15, avg_logFC) %>%
      mutate(show_genes = gene) %>%
      dplyr::select(gene, show_genes) %>%
      bind_rows(top5_wt) 
    
    micr <- micr %>%
      left_join(top_5_both) %>%
      mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_logFC) > .25, "sig.", "not sig."))
    
    
    micr_volcano <- ggplot(micr, aes(x=avg_logFC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
      geom_point(size=5) + #aes(size=avg_logFC)
      ggrepel::geom_text_repel(size=7, box.padding=1.15) +
      #expand_limits(y=c(0, 350), x=c(-2.25, 2.25)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            text = element_text(size=25)) +
      scale_color_manual(values = c("light grey", "black")) +
      NoLegend() +
      labs(x="avg. logFC", y="-log10 transf. adj. p-value")
    
    micr_volcano 

# Extended Data Figure 1a - donut plot
    metadata %>%
      group_by(Diagnosis,  seurat_clusters ) %>%
      summarise(freq = length( seurat_clusters )) %>%
      ggplot(aes(x=2, y=freq,fill= seurat_clusters )) +
      geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
      coord_polar(theta='y', start=0) +
      theme_void() +
      scale_fill_manual(values = rev(colors_pat)) +
      facet_wrap(~Diagnosis) +
      xlim(0.2, 2.5)
    
# Extended Data Figure 1b - violin plots
    data_t <- as.data.frame(t(as.matrix(all[["SCT"]]@counts)[,colnames(all)[all$Celltype =="Micr"]]))
    
    data_t <- data_t %>% merge(metadata, by=0)
    
    signature_genes <- data.frame("ahr"=c("Ahr", "Ahrr", "Cyp1a1", "Cyp1b1", "Entpd1", "Hif1a", "Il1b","Il27", "Klf4", "Pparg", "Stat1", "Stat3", "Tiparp", "Vegfa"),
                                  "APC"=c("Cd74", "Cd80", "Cd86", "H2-Aa", "Cd40", rep(NA, 9)),
                                  "Homeo"=c('P2ry12', 'Cx3cr1', 'Csf1r', 'Tmem119', 'Slc2a5', rep(NA, 9)), stringsAsFactors = F)
    
    signature_genes <- apply(signature_genes, 2, toupper) %>% as.data.frame()
    
    terms <- colnames(signature_genes)
    stats <- list()
    
    #stat testing and plotting            
    for (i in terms) {
      data_t2 <- data_t[, colnames(data_t) %in% na.omit(signature_genes[[i]])] %>% rowSums()
      data_t3 <- data_t[,c("seurat_clusters", "Diagnosis")] %>%
        bind_cols(data.frame("gene" = data_t2))
      plt <- ggplot(na.omit(data_t3), aes(seurat_clusters, log2(gene), fill = Diagnosis)) +
        geom_violin(scale = 'width', lwd=0.25, draw_quantiles = 0.5) +
        #geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(1))+
        scale_fill_brewer(palette = "Accent") +
        theme_minimal() +
        labs(y='log2(normalized counts)', title=i) 
      print(plt)
      
      #from urls: https://www.researchgate.net/post/Should_I_use_post-hoc_tukey_HSD_for_pairwise_comparisons_of_a_factor_on_a_zero-inflated_negative_binomial_mixed_models_ZINB
      #url:https://stats.stackexchange.com/questions/60352/comparing-levels-of-factors-after-a-glm-in-r
      mod <- glm.nb(gene ~ seurat_clusters * Diagnosis, data = na.omit(data_t3))
      
      emms <- emmeans(mod, ~ seurat_clusters : Diagnosis)
      comps <- pairs(emms, adjust="tukey") 
      stats[[i]] <- as.data.frame(comps)
      
    }
    
    stats_all <- bind_rows(stats, .id = "Signature")
    
    stats2 <- stats_all %>% 
      dplyr::filter(Signature %in% c("ahr", "APC", 'Homeo'), p.value < 0.05) 
    
    stats2$comparison <- gsub(",| - ", " ", stats2$contrast)
    
    pattern = START %R% capture(one_or_more(DGT)) %R% SPC %R% one_or_more(WRD) %R% SPC %R% capture(one_or_more(DGT))
    
    ind <- str_match(stats2$comparison, pattern=pattern)
    ind <- which(ind[,2] == ind[,3])
    
    #statistical testing for APC signatures significant within cluster comparisons diagnosis
    stats3 <- stats2[ind,]
    a <- stats3[stats3$Signature == "APC",] %>%
      mutate(Significance = ifelse(.$p.value<0.05 & .$p.value>0.01, '*',
                                   ifelse(.$p.value<0.01 & .$p.value>0.001, '**',
                                          ifelse(.$p.value<0.001, '***','n.s.'))))
     
    stats1 <- stats_all %>% 
      filter(Signature %in% c("ahr", "APC", "Homeo"), p.value < 0.05)
    
    write_csv(stats1, paste0("data/",date,"-statistical testing for APC and ahr signatures across cluster comparisons.csv"))
    
    stats3 <- stats2[ind,]
    
    #plot singnificant clusters
    clusters <- unique(as.numeric(gsub(" .*", "", stats3[stats3$Signature == "APC",]$comparison)))
    
    data_t2 <- rowSums(data_t[, colnames(data_t) %in% na.omit(signature_genes[["APC"]])])
    data_t3 <- data_t[,c("seurat_clusters", "Diagnosis")] %>%
      bind_cols(data.frame("gene" = data_t2))
    plt <- ggplot(na.omit(data_t3[data_t3$seurat_clusters %in% as.character(clusters),]), aes(seurat_clusters, log2(gene), fill = Diagnosis)) + #1:3,7:9,12
      geom_violin(scale = 'width',draw_quantiles =0.5) +
      #geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(0.9))+
      scale_fill_brewer(palette = "Accent") +
      theme_minimal() +
      labs(y='log2(normalized counts)', title="APC signature") 
    print(plt)
    ggsave(paste0('plots/others/',date,'-human-APC-signature-violin-plot-significant-cluster-Diagnosis.pdf'))
    

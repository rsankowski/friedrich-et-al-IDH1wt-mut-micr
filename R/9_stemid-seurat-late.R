library(RaceID)
library(Seurat)
library(tidyverse)
library(assertthat)
library(future)
#StemID
library(FateID)
library(readxl)
library(limma)
library(MASS)
library(RColorBrewer)

load("data/metadata.RData")
rownames(metadata) <- metadata$ID
date <- Sys.Date()

source('bin/functions.R')

begin <- Sys.time()

set.seed(79106)

if (!file.exists("data/ltr_late.RData")){ #data/ltr-larger-clusters.RData
  ltr <- Ltree(sc)
  
  #convert clusters in integers
  ltr@sc@cpart <- as.numeric(as.character(ltr@sc@cpart))
  names(ltr@sc@cpart) <- colnames(ltr@sc@ndata)
  
  ltr <- compentropy(ltr)
  ltr <- projcells(ltr,nmode=TRUE,fr=FALSE) #400
  ltr <- projback(ltr,pdishuf=100)
  ltr <- lineagegraph(ltr)
  ltr <- comppvalue(ltr,pthr=0.01)
  
  save(ltr, file = 'data/ltr.RData')
} else {
  load('data/ltr_late.RData')
}

ltr@sc@fcol <- sample(rainbow(max(ltr@sc@cpart)))

end <- Sys.time()

timediff <- end-begin
plotgraph(ltr,scthr=0.9,showCells=FALSE)

postscript(paste0('plots/umap/lineage-graph-tsne-plot-late.ps'), width = 8.57, height = 5.79)
plotgraph(ltr,scthr=0.9,showCells=FALSE)
dev.off()

pdf(paste0('plots/umap/lineage-graph-tsne-plot-late.pdf'), width = 8.57, height = 5.79)
plotgraph(ltr,scthr=0.9,showCells=FALSE,)
dev.off()

x <- compscore(ltr,scthr=0.8)
#plotdistanceratio(ltr)
#plotspantree(ltr)
#plotprojections(ltr)

#lineage tree for moDCs
micr <- c(2,4,3,8,10)
mon_macs <- c(34,15,11) #c(34,15,19,20)
mon_dcs <- c(34, 15,17,25)

#pseudotemporal micr
                    n <- cellsfromtree(ltr,micr)
                    x <- getfdata(ltr@sc)
                    
                    fs  <- filterset(x,n=n$f, minexpr = 2, minnumber = 10)
                    
                    if (!file.exists("data/s1d-micr_late.Robj")) {
                      s1d <- getsom(fs,nb=1000,alpha=.5)
                      save(s1d, file = "data/s1d-micr_late.Robj")
                    } else {
                      load("data/s1d-micr_late.Robj")
                    }
                    ps  <- procsom(s1d,corthr=.85,minsom=3)
                    y    <- ltr@sc@cpart[n$f]
                    fcol <- ltr@sc@fcol
                    plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                    
                    pdf(paste0('plots/heatmaps/', date, '-micr_late-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
                    plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                    dev.off()
                    
                    postscript(paste0('plots/heatmaps/', date, '-micr_late-trajectory-heatmap.ps'), width = 8.57, height = 5.79)
                    plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                    dev.off()
                    
                    png(paste0('plots/heatmaps/', date, '-micr_late-trajectory-heatmap.png'))
                    image(t(as.matrix(ps$all.z)), col = rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100)), axes = FALSE, ylim = c(-0.02,1))
                    dev.off()
                    
                    postscript(paste0('plots/heatmaps/', date, 'micr_late-trajectory-heatmap-labels-only.ps'), width = 8.57, height = 5.79)
                    plotheatmap2(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                    dev.off()
                    
                    ggplot(metadata[n$f, ], aes(x=1:length(n$f), y=0, color=Cluster)) +
                      geom_point(pch=15) +
                      theme_void()
                    
                    ggsave(paste0('plots/heatmaps/', date, 'micr_late-trajectory-cluster_cols.ps'))
                    
                    #export node genes
                    modules <- data.frame('Node' = NA, 'Genes' = NA)
                    for (i in 1: max(ps$nodes)) {
                      gene_names <- names(ps$nodes)[ps$nodes == i]
                      gene_names <- gsub('_.*', '', gene_names)
                      modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
                      modules2$Node <- rep(as.character(i), nrow(modules2))
                      modules <- rbind(na.omit(modules), modules2)
                    }
                    
                    write_csv(modules, paste0('data/',date,'-nodes-stemid-vector-micr_late.csv'))
                    
                    #plot trajectory composition
                    ps  <- procsom(s1d,corthr=.85,minsom=3)
                    traj <- data.frame(ID = colnames(ps$all.z), Trajectory = 1: length(colnames(ps$all.z))) %>%
                      left_join(metadata) %>%
                      mutate(step = ceiling(Trajectory/(max(Trajectory)/30)))
                    
                    traj %>%
                      ggplot(aes(x=step, fill=Genotype)) +
                      geom_bar(position = "fill", width = .95) +
                      scale_fill_brewer(palette = "Set1") +
                      theme_void()
                    
                    ggsave("plots/others/trajectory-micr_late.pdf", width = 8.57, height = 5.79)
                    
#pseudotemporal ordering - mon_macs
                  n <- cellsfromtree(ltr,mon_macs)
                  x <- getfdata(ltr@sc)
                  
                  
                  fs  <- filterset(x,n=n$f, minexpr = 2, minnumber = 10)
                  if (!file.exists("data/s1d-mon_macs_late_11.Robj")) {
                    s1d <- getsom(fs,nb=1000,alpha=.5)
                    save(s1d, file = "data/s1d-mon_macs_late_11.Robj")
                  } else {
                    load("data/s1d-mon_macs_late_11.Robj")
                  }
                  
                  ps  <- procsom(s1d,corthr=.85,minsom=3)
                  y    <- ltr@sc@cpart[n$f]
                  fcol <- ltr@sc@fcol
                  
                  save(ps, file = "data/late_moMFs_trajectory_heatmap_data.RData")
                  
                  plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  
                  pdf(paste0('plots/heatmaps/', date, '-mon_macs_late-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
                  plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  postscript(paste0('plots/heatmaps/', date, '-mon_macs_late-trajectory-heatmap.ps'), width = 8.57, height = 5.79)
                  plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  png(paste0('plots/heatmaps/', date, '-mon_macs_late-trajectory-heatmap.png'))
                  image(t(as.matrix(ps$all.z)), col = rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100)), axes = FALSE, ylim = c(-0.02,1))
                  dev.off()
                  
                  postscript(paste0('plots/heatmaps/', date, '-mon_macs_late-trajectory-heatmap-labels-only.ps'), width = 8.57, height = 5.79)
                  plotheatmap2(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  
                  ggplot(metadata[n$f, ], aes(x=1:length(n$f), y=0, color=Cluster)) +
                    geom_point(pch=15) +
                    theme_void()
                  
                  ggsave(paste0('plots/heatmaps/', date, '-mon_macs_late-trajectory-cluster_cols.ps'))
                  
                  #export data for heatmap fig 3e
                  write.csv(as.data.frame(ps$all.z), "data/fig3e_data_mon_late_heatmap.csv")
                  
                  #eport node genes
                  modules <- data.frame('Node' = NA, 'Genes' = NA)
                  for (i in 1: max(ps$nodes)) {
                    gene_names <- names(ps$nodes)[ps$nodes == i]
                    gene_names <- gsub('_.*', '', gene_names)
                    modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
                    modules2$Node <- rep(as.character(i), nrow(modules2))
                    modules <- rbind(na.omit(modules), modules2)
                  }
                  
                  write.csv(modules, paste0('data/',date,'-nodes-stemid-vector-mon_macs_late.csv'))
                  write.csv(modules, paste0('data/fig3e_gene_module_data.csv'))
                  
              
                  #plot trajectory composition
                  ps  <- procsom(s1d,corthr=.85,minsom=3)
                  
                  traj <- data.frame(ID = colnames(ps$all.z), Trajectory = 1: length(colnames(ps$all.z))) %>%
                    left_join(metadata) %>%
                    mutate(step = ceiling(Trajectory/(max(Trajectory)/30)))
                  
                  traj %>%
                    ggplot(aes(x=step, fill=Genotype)) +
                    geom_bar(position = "fill", width = .95) +
                    scale_fill_brewer(palette = "Set1") +
                    theme_void()
                  
                  ggsave("plots/others/trajectory-mon_macs_late_34_15_11.pdf", width = 8.57, height = 5.79)
                  
                  #export traj data
                  traj2 <- traj %>%
                    group_by(step) %>%
                    count(Genotype) %>%
                    mutate(n=n/17)
                  
                  write.csv(traj2, "data/fig3e_correlation_plot.csv")
#pseudotemporal ordering - mon_dcs
                  n <- cellsfromtree(ltr,mon_dcs)
                  x <- getfdata(ltr@sc)
                  
                  
                  fs  <- filterset(x,n=n$f, minexpr = 2, minnumber = 10)
                  if (!file.exists("data/s1d-mo_dcs_late.Robj")) {
                    s1d <- getsom(fs,nb=1000,alpha=.5)
                    save(s1d, file = "data/s1d-mo_dcs_late.Robj")
                  } else {
                    load("data/s1d-mo_dcs_late.Robj")
                  }
                  
                  ps  <- procsom(s1d,corthr=.85,minsom=3)
                  y    <- ltr@sc@cpart[n$f]
                  fcol <- ltr@sc@fcol
                  plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  
                  pdf(paste0('plots/heatmaps/', date, '-mo_dcs_late-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
                  plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  postscript(paste0('plots/heatmaps/', date, '-mo_dcs_late-trajectory-heatmap.ps'), width = 8.57, height = 5.79)
                  plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  #export data for fig 3f
                  write.csv(as.data.frame(ps$all.z), "data/fig3f_data_mon_modc_late_heatmap.csv")
                  
                  #eport node genes
                  modules <- data.frame('Node' = NA, 'Genes' = NA)
                  
                  for (i in 1: max(ps$nodes)) {
                    gene_names <- names(ps$nodes)[ps$nodes == i]
                    gene_names <- gsub('_.*', '', gene_names)
                    modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
                    modules2$Node <- rep(as.character(i), nrow(modules2))
                    modules <- rbind(na.omit(modules), modules2)
                  }
                  
                  write.csv(modules, paste0('data/',date,'-nodes-stemid-vector-mo_dcs_late.csv'))
                  write.csv(modules, "data/fig3f_data_mon_modc_late_gene_module_information.csv")
                  
                  png(paste0('plots/heatmaps/', date, '-mo_dcs_late-trajectory-heatmap.png'))
                  image(t(as.matrix(ps$all.z)), col = rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100)), axes = FALSE, ylim = c(-0.02,1))
                  dev.off()
                  
                  postscript(paste0('plots/heatmaps/', date, '-mo_dcs_late-trajectory-heatmap-labels-only.ps'), width = 8.57, height = 5.79)
                  plotheatmap2(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  ggplot(metadata[n$f, ], aes(x=1:length(n$f), y=0, color=Cluster)) +
                    geom_point(pch=15) +
                    theme_void()
                  
                  ggsave(paste0('plots/heatmaps/', date, '-mon_dcs_late-trajectory-cluster_cols.ps'))
                  
                   
                  #plot trajectory composition
                  ps  <- procsom(s1d,corthr=.85,minsom=3)
                  
                  traj <- data.frame(ID = colnames(ps$all.z), Trajectory = 1: length(colnames(ps$all.z))) %>%
                    left_join(metadata) %>%
                    mutate(step = ceiling(Trajectory/(max(Trajectory)/30)))
                  
                  traj %>%
                    ggplot(aes(x=step, fill=Genotype)) +
                    geom_bar(position = "fill", width = .95) +
                    scale_fill_brewer(palette = "Set1") +
                    theme_void()
                  
                  ggsave("plots/others/trajectory-mo_dcs_late.pdf", width = 8.57, height = 5.79)
                  
                  #export traj data
                  traj2 <- traj %>%
                    group_by(step) %>%
                    count(Genotype) %>%
                    mutate(n=n/27)
                  
                  write.csv(traj2, "data/fig3f_correlation_plot.csv")
                  
#plot genes along trajectory; genes downloaded from url: http://www.informatics.jax.org/go/term/GO:0006865 on 24 apr 2020
                  genes <- read_delim("data/GO_term_summary_20200424_110924.txt", delim = "\t")[["Symbol"]] %>%
                    na.omit() %>%
                  str_split( "\\|") %>%
                    unlist() %>%
                    unique()
                  
                  genes <- genes[genes %in% rownames(all[["SCT"]]@counts)]

#plot these genes along the microglia trajectory           
                  load("data/s1d-micr_late.Robj")
                  cells <- colnames(s1d$z)
                  clusters <- droplevels(all$seurat_clusters[cells])
                  
                  cts_micr <- all[["SCT"]]@counts[genes, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters,
                                     trajectory = 1:nrow(.)))
                  
                  cts_micr_l <- cts_micr %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Abat:Xk, names_to = "Gene", values_to = "Expr") 
                  
                  #%>% nest(data=c(ID,Expr)) %>%
                      #loess(Expr ~ 1:length(), data=., span=0.50)
                  
                  for (i in unique(cts_micr_l$Gene)) {
                    traj <- cts_micr_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-mon_late_no_dots.pdf"), useDingbats=F)
                    
                    traj <- cts_micr_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                  print(traj)
                  ggsave(paste0("plots/others/trajectory_plots/",i,"-micr_late.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_micr %>%
                    mutate(all_genes = rowSums(.)) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/cumulative_go_terms_micr_late.pdf", useDingbats=F)
                  
                  #expression of Slc13a3, Slc7a5 and Slc3a2 
                  genes2 <- c("Slc13a3", "Slc7a5", "Slc3a2")
                  
                  cts_micr <- all[["SCT"]]@counts[genes2, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters,
                                     trajectory = 1:nrow(.)))
                  
                  cts_micr_l <- cts_micr %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Slc13a3:Slc3a2, names_to = "Gene", values_to = "Expr") 
                  
                  #plot
                  cts_micr_l %>%
                    ggplot(aes(x=1:nrow(.), y=Expr)) +
                    geom_smooth(aes(color=Gene),se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_micr_late.pdf", useDingbats=F)
                  
                  
                  cts_micr_l %>%
                    ggplot(aes(x=1:nrow(.), y=Expr+0.1)) +
                    geom_smooth(aes(color=Gene),se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() +
                    scale_y_log10() 
                  ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_micr_late.pdf", useDingbats=F)
                  
                  
                  #fit loess
                  mod <- map(genes2, function(x) {
                    loess(Expr ~ trajectory, 
                          data = cts_micr_l[cts_micr_l$Gene == x,],
                          span = 0.75) %>%
                      broom::augment()
                  }) 
                  
                  names(mod) <- genes2
                  mod <- mod %>%
                    bind_rows(.id="gene") 
                  
                  mod %>%
                    write_csv("data/loess_fit_micr_late.csv")
                  
#monmacs
                  load("data/s1d-mon_macs_late_11.Robj")
                  cells <- colnames(s1d$z)
                  
                  clusters <- droplevels(all$seurat_clusters[cells])
                  
                  cts_mon <- all[["SCT"]]@counts[genes, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters,
                                     trajectory = 1:nrow(.)))
                  
                  cts_mon_l <- cts_mon %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Abat:Xk, names_to = "Gene", values_to = "Expr") 
                  
                  #%>% nest(data=c(ID,Expr)) %>%
                  #loess(Expr ~ 1:length(), data=., span=0.50)
                  
                  for (i in unique(cts_mon_l$Gene)) {
                    traj <- cts_mon_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-mon_late_no_dots.pdf"), useDingbats=F)
                    
                    traj <- cts_mon_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-mon_late_with_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_mon %>%
                    mutate(all_genes = rowSums(.[-120])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/cumulative_go_terms_mon_late.pdf", useDingbats=F)
                  
                  #expression of Slc13a3, Slc7a5 and Slc3a2 
                  genes2 <- c("Slc13a3", "Slc7a5", "Slc3a2")
                  
                  cts_mon <- all[["SCT"]]@counts[genes2, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters,
                                     trajectory = 1:nrow(.)))
                  
                  cts_mon_l <- cts_mon %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Slc13a3:Slc3a2, names_to = "Gene", values_to = "Expr") 
                  
                  #plot
                  cts_mon_l %>%
                    ggplot(aes(x=1:nrow(.), y=Expr)) +
                    geom_smooth(aes(color=Gene),se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_macs_late.pdf", useDingbats=F)
                  
                  #fit loess
                  mod <- map(genes2, function(x) {
                    loess(Expr ~ trajectory, 
                          data = cts_mon_l[cts_mon_l$Gene == x,],
                          span = 0.75) %>%
                      broom::augment()
                  }) 
                  
                  names(mod) <- genes2
                  mod <- mod %>%
                    bind_rows(.id="gene") 
                  
                  mod %>%
                    write_csv("data/loess_fit_mon_late.csv")
                  
#modcs
                  load("data/s1d-mo_dcs_late.Robj")
                  cells <- colnames(s1d$z)
                  
                  clusters <- droplevels(all$seurat_clusters[cells])
                  
                  cts_modc <- all[["SCT"]]@counts[genes, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters,
                                     trajectory = 1:nrow(.)))
                  
                  cts_modc_l <- cts_modc %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Abat:Xk, names_to = "Gene", values_to = "Expr") 
                  
                  #%>% nest(data=c(ID,Expr)) %>%
                  #loess(Expr ~ 1:length(), data=., span=0.50)
                  
                  for (i in unique(cts_modc_l$Gene)) {
                    traj <- cts_modc_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-modc_late_no_dots.pdf"), useDingbats=F)
                    
                    traj <- cts_modc_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-modc_late_with_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_modc %>%
                    mutate(all_genes = rowSums(.[,-120])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/cumulative_go_terms_modc_late.pdf", useDingbats=F)
                  
                  #expression of Slc13a3, Slc7a5 and Slc3a2 
                  genes2 <- c("Slc13a3", "Slc7a5", "Slc3a2")
                  
                  cts_modc <- all[["SCT"]]@counts[genes2, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters,
                                     trajectory = 1:nrow(.)))
                  
                  cts_modc_l <- cts_modc %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Slc13a3:Slc3a2, names_to = "Gene", values_to = "Expr") 
                  
                  #plot
                  cts_modc_l %>%
                    ggplot(aes(x=1:nrow(.), y=Expr)) +
                    geom_smooth(aes(color=Gene),se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_modc_late.pdf", useDingbats=F)
                  
                  #heatmap
                  cells <- colnames(s1d$z)
                  cts_modc <- all[["SCT"]]@counts[genes, cells] %>%
                    as.matrix() #%>%
                    #t %>%
                    #as.data.frame()
                  
                  pheatmap::pheatmap(cts_modc[apply(cts_modc,1, max)>0,], cluster_cols = F)
                  
                  fits <- cts_modc[apply(cts_modc,1, max)>0,] %>%
                    split(factor(rownames(cts_modc[apply(cts_modc,1, max)>0,]))) %>%
                    map(function(x) { 
                    mod <- loess(x~1:815)
                    predict(mod, x)
                    })
                  
                  cts_modc_l <- cts_modc %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Abat:Xk, names_to = "Gene", values_to = "Expr") 
                  
                  #%>% nest(data=c(ID,Expr)) %>%
                  #loess(Expr ~ 1:length(), data=., span=0.50)
                  
                  for (i in unique(cts_modc_l$Gene)) {
                    traj <- cts_modc_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-modc_late.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_modc %>%
                    mutate(all_genes = rowSums(.)) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/cumulative_go_terms_modc_late.pdf", useDingbats=F)
                  
                  
                  #fit loess
                  mod <- map(genes2, function(x) {
                    loess(Expr ~ trajectory, 
                          data = cts_modc_l[cts_modc_l$Gene == x,],
                          span = 10) %>%
                      broom::augment()
                  }) 
                  
                  names(mod) <- genes2
                  mod <- mod %>%
                    bind_rows(.id="gene") 
                  
                  mod %>%
                    write_csv("data/loess_fit_modc_late.csv")
                  
                  
#plot ahr signature
                  genes <- c("Ahr", "Ahrr", "Cyp1a1", "Cyp1b1", "Entpd1", "Hif1a", "Il1b", "Il27", "Klf4", "Pparg", "Stat1", "Stat3", "Tiparp", "Vegfa" )
                  
                  genes <- genes[genes %in% rownames(all[["SCT"]]@counts)]
                  
                  #plot these genes along the microglia trajectory           
                  load("data/s1d-micr_late.Robj")
                  cells <- colnames(s1d$z)
                  
                  clusters <- droplevels(all$seurat_clusters[cells])
                  
                  cts_micr <- all[["SCT"]]@counts[genes, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters))
                  
                  cts_micr_l <- cts_micr %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Ahr:Vegfa, names_to = "Gene", values_to = "Expr") 
                  
                  #%>% nest(data=c(ID,Expr)) %>%
                  #loess(Expr ~ 1:length(), data=., span=0.50)
                  
                  for (i in unique(cts_micr_l$Gene)) {
                    traj <- cts_micr_l[cts_micr_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-micr_late_with_dots.pdf"), useDingbats=F)
                    
                    #no dots
                    traj <- cts_micr_l[cts_micr_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-micr_late_no_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_micr %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    #geom_jitter(aes(color=cluster), size=5, alpha=0.5, height = .2, width = .2) +
                    geom_point(aes(color=cluster), size=5, alpha=0.5) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression - ahr", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_micr_late.pdf", useDingbats=F)
                  
                  #no dots
                  cts_micr %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression - ahr", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_micr_late_no_dots.pdf", useDingbats=F)
                  
                  #monmacs
                  load("data/s1d-mon_macs_late.Robj")
                  cells <- colnames(s1d$z)
                  
                  clusters <- droplevels(all$seurat_clusters[cells])
                  
                  cts_mon <- all[["SCT"]]@counts[genes, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters))
                  
                  cts_mon_l <- cts_mon %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Ahr:Vegfa, names_to = "Gene", values_to = "Expr") 
                  
                  #%>% nest(data=c(ID,Expr)) %>%
                  #loess(Expr ~ 1:length(), data=., span=0.50)
                  
                  for (i in unique(cts_mon_l$Gene)) {
                    traj <- cts_mon_l[cts_mon_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-mon_late_with_dots.pdf"), useDingbats=F)
                    
                    #no dots
                    traj <- cts_mon_l[cts_mon_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-mon_late_no_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_mon %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_mon_late_no_dots.pdf", useDingbats=F)
                  
                  #no dots
                  cts_mon %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_point(aes(color=cluster), size=5, alpha=0.5) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression - ahr", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_mon_late_with_dots.pdf", useDingbats=F)
                  
                  #modcs
                  load("data/s1d-mo_dcs_late.Robj")
                  cells <- colnames(s1d$z)
                  
                  clusters <- droplevels(all$seurat_clusters[cells])
                  
                  cts_modc <- all[["SCT"]]@counts[genes, cells] %>%
                    as.matrix() %>%
                    t %>%
                    as.data.frame() %>%
                    cbind(data.frame(cluster=clusters))
                  
                  cts_modc_l <- cts_modc %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Ahr:Vegfa, names_to = "Gene", values_to = "Expr") 
                  
                  #%>% nest(data=c(ID,Expr)) %>%
                  #loess(Expr ~ 1:length(), data=., span=0.50)
                  
                  for (i in unique(cts_modc_l$Gene)) {
                    traj <- cts_modc_l[cts_modc_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-modc_late_with_dots.pdf"), useDingbats=F)
                    
                    #no dots
                    traj <- cts_modc_l[cts_modc_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-modc_late_no_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_modc %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_modc_late_no_dots.pdf", useDingbats=F)
                  
                  #no dots
                  cts_modc %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_point(aes(color=cluster), size=5, alpha=0.5) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression - ahr", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_modc_late_with_dots.pdf", useDingbats=F)
                  
                  #heatmap
                  cells <- colnames(s1d$z)
                  cts_modc <- all[["SCT"]]@counts[genes, cells] %>%
                    as.matrix() #%>%
                  #t %>%
                  #as.data.frame()
                  
                  pheatmap::pheatmap(cts_modc[apply(cts_modc,1, max)>0,], cluster_cols = F)
                  
                  fits <- cts_modc[apply(cts_modc,1, max)>0,] %>%
                    split(factor(rownames(cts_modc[apply(cts_modc,1, max)>0,]))) %>%
                    map(function(x) { 
                      mod <- loess(x~1:815)
                      predict(mod, x)
                    })
                  
                  cts_modc_l <- cts_modc %>% 
                    rownames_to_column(var = "ID") %>%
                    pivot_longer(Abat:Xk, names_to = "Gene", values_to = "Expr") 
                  
                  #%>% nest(data=c(ID,Expr)) %>%
                  #loess(Expr ~ 1:length(), data=., span=0.50)
                  
                  for (i in unique(cts_modc_l$Gene)) {
                    traj <- cts_modc_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-modc_late.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_modc %>%
                    mutate(all_genes = rowSums(.)) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_modc_late.pdf", useDingbats=F)
                  
                  
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


date <- Sys.Date()

source('bin/functions.R')

begin <- Sys.time()

if (!file.exists("data/ltr_early.RData")){ #data/ltr-larger-clusters.RData
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
  load('data/ltr_early.RData')
}

end <- Sys.time()

timediff <- end-begin
plotgraph(ltr,scthr=0.9,showCells=FALSE)

pdf(paste0('plots/umap/lineage-graph-early.pdf'), width = 8.57, height = 5.79, useDingbats=F)
plotgraph(ltr,scthr=0.9,showCells=FALSE)
dev.off()


x <- compscore(ltr,scthr=0.9)
#plotdistanceratio(ltr)
#plotspantree(ltr)
#plotprojections(ltr)

#lineage tree for moDCs
micr <- c(2,15,14,7,3,8,10)
mo_macs <- c(36,16,20,21)
mo_dc <- c(36,16,18,27)

#pseudotemporal micr
                    n <- cellsfromtree(ltr,micr)
                    
                    fs  <- filterset(getfdata(ltr@sc),n=n$f, minexpr = 2, minnumber = 10)
                    
                    if (!file.exists("data/s1d-micr_early.Robj")) {
                      s1d <- getsom(fs,nb=1000,alpha=.5)
                      save(s1d, file = "data/s1d-micr_early.Robj")
                    } else {
                      load("data/s1d-micr_early.Robj")
                    }
                    ps  <- procsom(s1d,corthr=.85,minsom=3)
                    y    <- ltr@sc@cpart[n$f]
                    fcol <- sc@fcol
                    plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                    
                    pdf(paste0('plots/heatmaps/', date, '-micr-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
                    plotheatmap(ps$all_z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                    dev.off()
                    
                    png(paste0('plots/heatmaps/', date, '-micr-trajectory-heatmap.png'), width = 8.57, height = 5.79, res = 600, units = "in")
                    plotheatmap(ps$all_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                    dev.off()
                    
                    #export data for heatmap
                    heat <- as.data.frame(ps$all.z)
                    write.csv(heat, "data/fig_2j_heatmap_data.csv")
                    #export node genes
                    modules <- data.frame('Node' = NA, 'Genes' = NA)
                    for (i in 1: max(ps$nodes)) {
                      gene_names <- names(ps$nodes)[ps$nodes == i]
                      gene_names <- gsub('_.*', '', gene_names)
                      modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
                      modules2$Node <- rep(as.character(i), nrow(modules2))
                      modules <- rbind(na.omit(modules), modules2)
                    }
                    
                    write_csv(modules, paste0('data/',date,'-nodes-stemid-vector-micr-new.csv'))
                    write_csv(modules, paste0('data/fig_2j_gene_modules.csv'))
                    
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
                    
                    ggsave("plots/others/trajectory-micr_early.pdf", width = 8.57, height = 5.79)
                    
                    #export traj data
                    traj2 <- traj %>%
                      group_by(step) %>%
                      count(Genotype) %>%
                      mutate(n=n/143)
                    
                    write.csv(traj2, "data/fig2j_correlation_plot.csv")
                    
                    #plot genes along trajectory
                    #plot these genes along the microglia trajectory           
                    load("data/s1d-micr_early.Robj")
                    cells <- colnames(s1d$z)
                    clusters <- droplevels(all$seurat_clusters[cells])
                    
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
                    ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_micr_early.pdf", useDingbats=F)
                    
                    #export loess fit
                    a <- cts_micr_l %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      stat_smooth(aes( outfit=fit<<-..y.., color=Gene),se=F)
                    
                    cts_micr_l %>%
                      ggplot(aes(x=1:nrow(.), y=Expr+0.1)) +
                      geom_smooth(aes(color=Gene),se=F) +
                      labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                      theme_bw() +
                      scale_y_log10() 
                    ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_micr_early.pdf", useDingbats=F)
                    
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
                      write_csv("data/loess_fit_micr_early.csv")
                    
                    

#pseudotemporal ordering - mo_macs
                  n <- cellsfromtree(ltr,mo_macs)
                  x <- getfdata(ltr@sc)
                  
                  
                  fs  <- filterset(x,n=n$f, minexpr = 2, minnumber = 5)
                  if (!file.exists("data/s1d-mo_macs-new.Robj")) {
                    s1d <- getsom(fs,nb=1000,alpha=.5)
                    save(s1d, file = "data/s1d-mo_macs-new.Robj")
                  } else {
                    load("data/s1d-mo_macs.Robj")
                  }
                  
                  ps  <- procsom(s1d,corthr=.85,minsom=3)
                  y    <- ltr@sc@cpart[n$f]
                  fcol <- ltr@sc@fcol
                  plotheatmap(ps$all_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  
                  pdf(paste0('plots/heatmaps/', date, '-mo_macs-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
                  plotheatmap(ps$all_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  png(paste0('plots/heatmaps/', date, '-mo_macs-trajectory-heatmap.png'), width = 8.57, height = 5.79, res = 600, units = "in")
                  plotheatmap(ps$all_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  #eport node genes
                  modules <- data.frame('Node' = NA, 'Genes' = NA)
                  for (i in 1: max(ps$nodes)) {
                    gene_names <- names(ps$nodes)[ps$nodes == i]
                    gene_names <- gsub('_.*', '', gene_names)
                    modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
                    modules2$Node <- rep(as.character(i), nrow(modules2))
                    modules <- rbind(na.omit(modules), modules2)
                  }
                  
                  write.csv(modules, paste0('data/',date,'-nodes-stemid-vector-mo_macs.csv'))
                  
#pseudotemporal ordering - mo_dc
                  n <- cellsfromtree(ltr,mo_dc)
                  x <- getfdata(ltr@sc)
                  
                  
                  fs  <- filterset(x,n=n$f, minexpr = 2, minnumber = 5)
                  if (!file.exists("data/s1d-mo_dc-new.Robj")) {
                    s1d <- getsom(fs,nb=1000,alpha=.5)
                    save(s1d, file = "data/s1d-mo_dc-new.Robj")
                  } else {
                    load("data/s1d-mo_dc.Robj")
                  }
                  
                  ps  <- procsom(s1d,corthr=.85,minsom=3)
                  y    <- ltr@sc@cpart[n$f]
                  fcol <- ltr@sc@fcol
                  plotheatmap(ps$all_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  
                  pdf(paste0('plots/heatmaps/', date, '-mo_dc-trajectory-heatmap.pdf'), width = 8.57, height = 5.79)
                  plotheatmap(ps$all_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  png(paste0('plots/heatmaps/', date, '-mo_dc-trajectory-heatmap.png'), width = 8.57, height = 5.79, res = 600, units = "in")
                  plotheatmap(ps$all_ea.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
                  dev.off()
                  
                  #export node genes
                  modules <- data.frame('Node' = NA, 'Genes' = NA)
                  for (i in 1: max(ps$nodes)) {
                    gene_names <- names(ps$nodes)[ps$nodes == i]
                    gene_names <- gsub('_.*', '', gene_names)
                    modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
                    modules2$Node <- rep(as.character(i), nrow(modules2))
                    modules <- rbind(na.omit(modules), modules2)
                  }
                  
                  write.csv(modules, paste0('data/',date,'-nodes-stemid-vector-mo_dc.csv'))
                  
                 
                  #plot genes along trajectory; genes downloaded from url: http://www.informatics.jax.org/go/term/GO:0006865 on 24 apr 2020
                  genes <- read_delim("data/GO_term_summary_20200424_110924.txt", delim = "\t")[["Symbol"]] %>%
                    na.omit() %>%
                    str_split( "\\|") %>%
                    unlist() %>%
                    unique()
                  
                  genes <- genes[genes %in% rownames(all[["SCT"]]@counts)]
                  
                  #plot these genes along the microglia trajectory           
                  load("data/s1d-micr_early.Robj")
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
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-micr_early_no_dots.pdf"), useDingbats=F)
                    
                    traj <- cts_micr_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-micr_early.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_micr %>%
                    mutate(all_genes = rowSums(.)) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/cumulative_go_terms_micr_early.pdf", useDingbats=F)
                  
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
                  ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_micr_early.pdf", useDingbats=F)
                  
                  
                  cts_micr_l %>%
                    ggplot(aes(x=1:nrow(.), y=Expr+0.1)) +
                    geom_smooth(aes(color=Gene),se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() +
                    scale_y_log10() 
                  ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_micr_early.pdf", useDingbats=F)
                  
                  #fit loess
                  mod <- map(genes2, function(x) {
                    loess(Expr ~ trajectory, 
                          data = cts_micr_l[cts_micr_l$Gene == x,],
                          span = .75) %>%
                      broom::augment()
                  }) 
                  
                  names(mod) <- genes2
                  mod <- mod %>%
                    bind_rows(.id="gene") 
                  
                  mod %>%
                    write_csv("data/loess_fit_micr_early.csv")
                  
                  
                  #monmacs
                  load("data/s1d-mon_macs_early.Robj")
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
                  
                 
                  for (i in unique(cts_mon_l$Gene)) {
                    traj <- cts_mon_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-mon_early_no_dots.pdf"), useDingbats=F)
                    
                    traj <- cts_mon_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-mon_early_with_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_mon %>%
                    mutate(all_genes = rowSums(.[-120])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/cumulative_go_terms_mon_early.pdf", useDingbats=F)
                  
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
                  ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_macs_early.pdf", useDingbats=F)
                  
                  #fit loess
                  mod <- map(genes2, function(x) {
                    loess(Expr ~ trajectory, 
                          data = cts_mon_l[cts_mon_l$Gene == x,],
                          span = .75) %>%
                      broom::augment()
                  }) 
                  
                  names(mod) <- genes2
                  mod <- mod %>%
                    bind_rows(.id="gene") 
                  
                  mod %>%
                    write_csv("data/loess_fit_mon_early.csv")
                  
                  #modcs
                  load("data/s1d-mo_dcs_early.Robj")
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
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-modc_early_no_dots.pdf"), useDingbats=F)
                    
                    traj <- cts_modc_l %>%
                      filter(Gene == i) %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_point(aes(color=cluster), size=5, alpha=0.5) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-modc_early_with_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_modc %>%
                    mutate(all_genes = rowSums(.[,-120])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/cumulative_go_terms_modc_early.pdf", useDingbats=F)
                  
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
                  ggsave("plots/others/trajectory_plots/expr_slc13a3_etc_modc_early.pdf", useDingbats=F)
                  
                  #fit loess
                  mod <- map(genes2, function(x) {
                    loess(Expr ~ trajectory, 
                          data = cts_modc_l[cts_modc_l$Gene == x,],
                          span = 0.75) %>%
                      broom::augment()
                  }) 
                  
                  names(mod) <- genes2
                  mod <- mod %>%
                    bind_rows(.id="gene") 
                  
                  mod %>%
                    write_csv("data/loess_fit_modc_early.csv")
                  
                  
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
                    ggsave(paste0("plots/others/trajectory_plots/",i,"-modc_early.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_modc %>%
                    mutate(all_genes = rowSums(.)) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/trajectory_plots/cumulative_go_terms_modc_early.pdf", useDingbats=F)
                  
                  #plot ahr signature
                  genes <- c("Ahr", "Ahrr", "Cyp1a1", "Cyp1b1", "Entpd1", "Hif1a", "Il1b", "Il27", "Klf4", "Pparg", "Stat1", "Stat3", "Tiparp", "Vegfa" )
                  
                  genes <- genes[genes %in% rownames(all[["SCT"]]@counts)]
                  
                  #plot these genes along the microglia trajectory           
                  load("data/s1d-micr_early.Robj")
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
                    ggsave(paste0("plots/others/ahr/",i,"-micr_early_with_dots.pdf"), useDingbats=F)
                    
                    #no dots
                    traj <- cts_micr_l[cts_micr_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-micr_early_no_dots.pdf"), useDingbats=F)
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
                  ggsave("plots/others/ahr/cumulative_ahr_micr_early.pdf", useDingbats=F)
                  
                  #no dots
                  cts_micr %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression - ahr", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_micr_early_no_dots.pdf", useDingbats=F)
                  
                  #monmacs
                  load("data/s1d-mon_macs_early.Robj")
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
                    ggsave(paste0("plots/others/ahr/",i,"-mon_early_with_dots.pdf"), useDingbats=F)
                    
                    #no dots
                    traj <- cts_mon_l[cts_mon_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-mon_early_no_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_mon %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_mon_early_no_dots.pdf", useDingbats=F)
                  
                  #no dots
                  cts_mon %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_point(aes(color=cluster), size=5, alpha=0.5) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression - ahr", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_mon_early_with_dots.pdf", useDingbats=F)
                  
                  #modcs
                  load("data/s1d-mo_dcs_early.Robj")
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
                    ggsave(paste0("plots/others/ahr/",i,"-modc_early_with_dots.pdf"), useDingbats=F)
                    
                    #no dots
                    traj <- cts_modc_l[cts_modc_l$Gene ==i,] %>%
                      ggplot(aes(x=1:nrow(.), y=Expr)) +
                      geom_smooth(se=F) +
                      labs(title = i, x="trajectory", y="Expression") +
                      theme_bw() 
                    
                    print(traj)
                    ggsave(paste0("plots/others/ahr/",i,"-modc_early_no_dots.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_modc %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_modc_early_no_dots.pdf", useDingbats=F)
                  
                  #no dots
                  cts_modc %>%
                    mutate(all_genes = rowSums(.[,-13])) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_point(aes(color=cluster), size=5, alpha=0.5) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression - ahr", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_modc_early_with_dots.pdf", useDingbats=F)
                  
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
                    ggsave(paste0("plots/others/ahr/",i,"-modc_early.pdf"), useDingbats=F)
                  }
                  
                  #cumulative
                  cts_modc %>%
                    mutate(all_genes = rowSums(.)) %>% 
                    ggplot(aes(x=1:nrow(.), y=all_genes)) +
                    geom_smooth(se=F) +
                    labs(title = "Cumulative gene expression", x="trajectory", y="Expression") +
                    theme_bw() 
                  ggsave("plots/others/ahr/cumulative_ahr_modc_early.pdf", useDingbats=F)
                  
                  
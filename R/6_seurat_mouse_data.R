#RaceID4 
library(tidyverse)
library(viridis)
library(RaceID)
library(Seurat)
library(Matrix)
library(clustree)
library(cowplot)
library(assertthat)

date = Sys.Date()

doublets <- list()
doublets <- map(list.files(path = "data/", pattern = "_doublet_scores.csv"),function(x)  read_csv(file.path("data",x)))
names(doublets) <- gsub("_.*", "",list.files(path = "data/", pattern = "_doublet_scores.csv"))

first <- Read10X_h5('data/counts/more-late/Mm_CD45_IDH_SMAR_6mice_filtered_feature_bc_matrix.h5') 

#exclude doublets
plot(density(doublets[["first"]][[2]]))
assert_that(ncol(first) == nrow(doublets[["first"]]))
summary(doublets[["first"]])
first <- first[,which(doublets[["first"]][[2]] < quantile(doublets[["first"]][[2]], 0.94))]

first <- first %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)

#load new data
rh1 <- Read10X_h5("data/counts/Early CD45 GL261 10x run/ea_RH_1_v3_filtered_feature_bc_matrix.h5") 

#remove doublets
#exclude doublets
plot(density(doublets[["rh1"]][[2]]))
assert_that(ncol(rh1) == nrow(doublets[["rh1"]]))
summary(doublets[["rh1"]])
rh1 <- rh1[,which(doublets[["rh1"]][[2]] < quantile(doublets[["rh1"]][[2]], 0.94))]

rh1 <- rh1 %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)
rh1$id <- "nRH1"
rh1$genotype <- "RH"

rh2 <- Read10X_h5("data/counts/Early CD45 GL261 10x run/ea_RH_2_v3_filtered_feature_bc_matrix.h5") 

#exclude doublets
plot(density(doublets[["rh2"]][[2]]))
assert_that(ncol(rh2) == nrow(doublets[["rh2"]]))
summary(doublets[["rh2"]])
rh2 <- rh2[,which(doublets[["rh2"]][[2]] < quantile(doublets[["rh2"]][[2]], 0.94))]

rh2 <- rh2 %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)
rh2$id <- "nRH2"
rh2$genotype <- "RH"

rh3 <- Read10X_h5("data/counts/Early CD45 GL261 10x run/ea_RH_3_v3_filtered_feature_bc_matrix.h5") 

#exclude doublets
plot(density(doublets[["rh3"]][[2]]))
assert_that(ncol(rh3) == nrow(doublets[["rh3"]]))
summary(doublets[["rh3"]])
rh3 <- rh3[,which(doublets[["rh3"]][[2]] < quantile(doublets[["rh3"]][[2]], 0.94))]

rh3 <- rh3 %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)
rh3$id <- "nRH3"
rh3$genotype <- "RH"

rh4 <- Read10X_h5("data/counts/Early CD45 GL261 10x run/ea_RH_4_v3_filtered_feature_bc_matrix.h5") 

#exclude doublets
plot(density(doublets[["rh4"]][[2]]))
assert_that(ncol(rh4) == nrow(doublets[["rh4"]]))
summary(doublets[["rh4"]])
rh4 <- rh4[,which(doublets[["rh4"]][[2]] < quantile(doublets[["rh4"]][[2]], 0.94))]

rh4 <- rh4 %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)
rh4$id <- "nRH4"
rh4$genotype <- "RH"

wt1 <- Read10X_h5("data/counts/Early CD45 GL261 10x run/ea_WT_1_v3_filtered_feature_bc_matrix.h5") 

#exclude doublets
plot(density(doublets[["wt1"]][[2]]))
assert_that(ncol(wt1) == nrow(doublets[["wt1"]]))
summary(doublets[["wt1"]])
wt1 <- wt1[,which(doublets[["wt1"]][[2]] < quantile(doublets[["wt1"]][[2]], 0.94))]

wt1 <- wt1 %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)
wt1$id <- "nWT1"
wt1$genotype <- "WT"

wt2 <- Read10X_h5("data/counts/Early CD45 GL261 10x run/ea_WT_2_v3_filtered_feature_bc_matrix.h5") 

#exclude doublets
plot(density(doublets[["wt2"]][[2]]))
assert_that(ncol(wt2) == nrow(doublets[["wt2"]]))
summary(doublets[["wt2"]])
wt2 <- wt2[,which(doublets[["wt2"]][[2]] < quantile(doublets[["wt2"]][[2]], 0.94))]

wt2 <- wt2 %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)
wt2$id <- "nWT2"
wt2$genotype <- "WT"

wt3 <- Read10X_h5("data/counts/Early CD45 GL261 10x run/ea_WT_3_v3_filtered_feature_bc_matrix.h5") 

#exclude doublets
plot(density(doublets[["wt3"]][[2]]))
assert_that(ncol(wt3) == nrow(doublets[["wt3"]]))
summary(doublets[["wt3"]])
wt3 <- wt3[,which(doublets[["wt3"]][[2]] < quantile(doublets[["wt3"]][[2]], 0.94))]

wt3 <- wt3 %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)
wt3$id <- "nWT3"
wt3$genotype <- "WT"

wt4 <- Read10X_h5("data/counts/Early CD45 GL261 10x run/ea_WT_4_v3_filtered_feature_bc_matrix.h5") 

#exclude doublets
plot(density(doublets[["wt4"]][[2]]))
assert_that(ncol(wt4) == nrow(doublets[["wt4"]]))
summary(doublets[["wt4"]])
wt4 <- wt4[,which(doublets[["wt4"]][[2]] < quantile(doublets[["wt4"]][[2]], 0.94))]

wt4 <- wt4 %>%
  CreateSeuratObject(project="wt", min.cells = 5, min.features = 500) %>%
  subset(subset=nFeature_RNA>500 & nFeature_RNA < 3000)
wt4$id <- "nWT4"
wt4$genotype <- "WT"

new1 <- merge(x=rh1, y=c(rh2,wt1,wt2), add.cell.ids = c("ea_rh1", "ea_rh2", "ea_wt1","ea_wt2"),project="early2")
new1$timepoint <- "early"
new1$batch <- as.factor("round2")
new1[["percent.mt"]] <-PercentageFeatureSet(new1,pattern="^mt-")
new1[["dissociation"]] <-PercentageFeatureSet(new1,pattern="^(Fos|Jun|Hsp|Atf3|Zfp36|Socs3)")

new2 <- merge(x=rh3, y=c(rh4,wt3,wt4), add.cell.ids = c("ea_rh3", "ea_rh4", "ea_wt3", "ea_wt4"),project="early2")
new2$timepoint <- "early"
new2$batch <- as.factor("round2")
new2[["percent.mt"]] <-PercentageFeatureSet(new2,pattern="^mt-")
new2[["dissociation"]] <-PercentageFeatureSet(new2,pattern="^(Fos|Jun|Hsp|Atf3|Zfp36|Socs3)")

#exclude doublets between lymphocytes and microglia, granulocytes and microglia
cells_to_exclude <- c(read_csv("data/microglia-lymphocyte-doublets.csv")[[1]],
                      read_csv("data/microglia-granulocyte-doublets.csv")[[1]])

first <- first[, !colnames(first) %in% cells_to_exclude]
new1 <- new1[, !colnames(new1) %in% cells_to_exclude]
new2 <- new2[, !colnames(new2) %in% cells_to_exclude]


if (F){existing <- merge(early,y=late,add.cell.ids=c("ea","la"),project="wt")
existing$batch <- as.factor("round1")
existing[["percent.mt"]] <-PercentageFeatureSet(existing,pattern="^mt-")
existing$genotype <- "WT"
}

first$batch <- as.factor("round0")
first[["percent.mt"]] <-PercentageFeatureSet(first,pattern="^mt-")
first[["dissociation"]] <-PercentageFeatureSet(first,pattern="^(Fos|Jun|Hsp|Atf3|Zfp36|Socs3)")
first$genotype <- as.factor(gsub('[A-Z-]', '', colnames(first)))
levels(first$genotype) <- c('WT1', 'WT2', 'WT4', 'RH1', 'RH2', 'RH3')
first$genotype <- gsub("[0-9]", "", as.character(first$genotype))
first$id <- as.factor(gsub('[A-Z-]', '', colnames(first)))
levels(first$id) <- c('fWT1', 'fWT2', 'fWT4', 'fRH1', 'fRH2', 'fRH3')
first$timepoint <- "late"

  #remove single objects
  rm(rh1,rh2,rh3,rh4,wt1,wt2,wt3,wt4)
  
  #sctransform
  first <-  SCTransform(first, vars.to.regress = c("percent.mt"), #,"dissociation"
                        variable.features.n = 10000)
  #save.image("data/checkpoint1.RData")
  new1 <-  SCTransform(new1, vars.to.regress = c("percent.mt"), #,"dissociation"
                       variable.features.n = 10000)
  #save.image("data/checkpoint2.RData")
  new2 <-  SCTransform(new2, vars.to.regress = c("percent.mt"), #,"dissociation"
                       variable.features.n = 10000)
  save.image("data/checkpoint3.RData")
  
  
  load("data/checkpoint3.RData")
  
  #data integration
  #set parallelization
 
  anchor_features <- grep("^(Fos|Jun|Gm|Rps|Rpl|Atf3|Zfp36|AY|Egr|Malat1|Xist|Hsp|mt-|Hist|Socs3)", SelectIntegrationFeatures(list(first, new1, new2)), invert = T, value = T)
  
  immune.anchors <- FindIntegrationAnchors(object.list = list(first, new1, new2),
                                           dims = 1:20,
                                           anchor.features = anchor_features) 
  immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
  
  #integrated analysis
  DefaultAssay(immune.combined) <- "integrated"
  
  rm(new1, new2, first)
  
  save.image("data/checkpoint4.RData")
  
  load("data/checkpoint4.RData")
  
  all <- immune.combined
  rm(immune.anchors, immune.combined)
 
  all <- ScaleData(all, verbose = T,
                   features = grep("^(Fos|Jun|Gm|Rps|Rpl|Atf3|Zfp36|AY|Egr|Malat1|Xist|Hsp|mt-|Hist|Socs3|Lars2)", VariableFeatures(all), invert = T, value = T))
  
  all <- RunPCA(all, 
                npcs = 20, 
                pc.genes =  grep("^(Fos|Jun|Gm|Rps|Rpl|Atf3|Zfp36|AY|Egr|Malat1|Xist|Hsp|mt-|Hist|Socs3|Lars2)", VariableFeatures(all), invert = T, value = T)) #, pc.genes =  grep("^(Fos|Jun|Gm|Rps|Rpl|Atf3|Zfp36|AY|Egr|Malat1|Xist|Hsp|mt-|Hist|Socs3)", VariableFeatures(all), invert = T, value = T)
  ElbowPlot(all)
  
  # UMAP and Clustering
  all <- RunUMAP(all, reduction = "pca", dims = 1:15)
  all <- FindNeighbors(all, reduction = "pca", dims = 1:15)
  all <- FindClusters(all, resolution = c(.1,.5,1,1.5,2,2.5,3))
  
  #url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
  clustree(all)
  
  all <- FindClusters(all, resolution = 2.5)
  
  #plot
  p1 <- DimPlot(all, reduction = "umap", group.by = "genotype")
  p2 <- DimPlot(all, reduction = "umap", label = TRUE)
  plot_grid(p1, p2)
  
  save(all, file = "data/seurat-integration-10x-wt-rh-gbm.RData")
  
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

load('data/prdata.RData')

#load metadata
metadata <- read.csv("data/human_metadata.csv")

#build batch effect table
batch1 <- prdata[,metadata$Batch =="Sequencing2"] %>%
  CreateSeuratObject(min.cells = 3, min.features = 500)
batch1$percent.mt <- PercentageFeatureSet(batch1,pattern="^MT-")
batch1$percent.rp <- PercentageFeatureSet(batch1,pattern="^(RPS|RPL)")
batch1$Age <- metadata[colnames(batch1), "Age"] 

batch2 <- prdata[,metadata$Batch =="Sequencing3"] %>%
  CreateSeuratObject(min.cells = 3, min.features = 500)
batch2$percent.mt <- PercentageFeatureSet(batch2,pattern="^MT-")
batch2$percent.rp <- PercentageFeatureSet(batch2,pattern="^(RPS|RPL)") 
batch2$Age <- metadata[colnames(batch2), "Age"] 

batch3 <- prdata[,metadata$Batch =="Sequencing5"] %>%
  CreateSeuratObject(min.cells = 3, min.features = 500)
batch3$percent.mt <- PercentageFeatureSet(batch3,pattern="^MT-")
batch3$percent.rp <- PercentageFeatureSet(batch3,pattern="^(RPS|RPL)")
batch3$Age <- metadata[colnames(batch3), "Age"] 

batch4 <- prdata[,metadata$Batch =="Sequencing6"] %>%
  CreateSeuratObject(min.cells = 3, min.features = 500)
batch4$percent.mt <- PercentageFeatureSet(batch4,pattern="^MT-")
batch4$percent.rp <- PercentageFeatureSet(batch4,pattern="^(RPS|RPL)")
batch4$Age <- metadata[colnames(batch4), "Age"] 

batch5 <- prdata[,metadata$Batch =="Sequencing9"] %>%
  CreateSeuratObject(min.cells = 3, min.features = 500)
batch5$percent.mt <- PercentageFeatureSet(batch5,pattern="^MT-")
batch5$percent.rp <- PercentageFeatureSet(batch5,pattern="^(RPS|RPL)")
batch5$Age <- metadata[colnames(batch5), "Age"] 

batch1$Batch <- "A"
batch2$Batch <- "B"
batch3$Batch <- "C"
batch4$Batch <- "D"
batch5$Batch <- "E"

micr_list <- list(batch1, batch2, batch3, batch4, batch5)

rm(batch1, batch2, batch3, batch4, batch5, prdata)

  
  #sctransform
  for (i in 1:length(micr_list)) {
    micr_list[[i]] <- SCTransform(micr_list[[i]], 
                                  variable.features.n = 10000,
                                  verbose = T,
                                  return.only.var.genes = F,
                                  vars.to.regress = c("percent.mt")
    )
  }
  
  #data integration
  #following url: https://satijalab.org/seurat/v3.0/integration.html sctransform workflow from April 16th, 2020
  anchor_features <- grep("^(FOS|JUN|RP|ZFP36|EGR|MALAT1|XIST|HSP|MT-|HIST)", SelectIntegrationFeatures(micr_list, nfeatures = 10000), invert = T, value = T)
  
  micr_list <- PrepSCTIntegration(object.list = micr_list, 
                                  anchor.features = anchor_features, 
                                  verbose = T)
  
  immune.anchors <- FindIntegrationAnchors(object.list = micr_list, 
                                           normalization.method = "SCT",
                                           dims = 1:20,
                                           anchor.features = anchor_features) 
  immune.combined <- IntegrateData(anchorset = immune.anchors, 
                                   normalization.method = "SCT", 
                                   dims = 1:20)
  
  #integrated analysis
  DefaultAssay(immune.combined) <- "integrated"
  
  
  all <- immune.combined
  rm(immune.anchors, immune.combined, micr_list)
  
  all <- RunPCA(all, 
                npcs = 20, 
                pc.genes =  grep("^(FOS|JUN|RP|ZFP36|EGR|MALAT1|XIST|HSP|MT-|HIST)", VariableFeatures(all), invert = T, value = T)
  ) 
  ElbowPlot(all)
  
  # t-SNE and Clustering
  all <- RunUMAP(all, reduction = "pca", dims = 1:15)
  all <- FindNeighbors(all, reduction = "pca", dims = 1:15)
  all <- FindClusters(all, resolution = c(.1,.2,.6,.8,1,1.5,2))
  
  #url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
  clustree(all)
  
  all <- FindClusters(all, resolution = .6)
  
  DefaultAssay(all) <- "SCT"
  
  #add metadata
  all$Diagnosis <- factor(ifelse(grepl("IDHwt",colnames(all)), "IDHwt_GBM", 
                                 ifelse(grepl("IDHmut",colnames(all)), "IDHmut_GBM", "Ctrl")),
                          levels = c("Ctrl", "IDHmut_GBM", "IDHwt_GBM"))
  
  all$Celltype <- all@active.ident
  levels(all$Celltype) <- c(rep("Micr", 7), "T cells", "NKT cells", "Glial cells", "Monocytes")
  all$Celltype <- factor(all$Celltype, levels = c("Micr", "Monocytes", "T cells",  "NKT cells", "Glial cells"))
  all$Age <- metadata[colnames(all),]$Age
  all$Localization <- as.factor(metadata[colnames(all),]$Lobe)
  all$Gender <- as.factor(metadata[colnames(all),]$Gender)
  all$ID <- as.factor(gsub("_.*", "", colnames(all)))
  
  #plot
  DimPlot(all, reduction = "umap", group.by = "Diagnosis")
  DimPlot(all, reduction = "umap", label = TRUE)
  DimPlot(all, reduction = "umap", group.by = "Celltype")
  
  
  save(all, file = "data/seurat-integration-wt-rh-gbm-final.RData")
 


#Sankowski et al -functions and plots

#colors
colors_many <- toupper(read_csv('data/20180313_trubetskoy_colors.csv')[[3]])
colors <- toupper(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'))
colors_pat <- toupper(c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5'))
colors_fig <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3', "light grey", "grey", "dark grey", "#696969")


#run RaceID
run_raceid <- function(.prdata = prdata,
                       .mintotal = 1500, 
                       .CGenes = c("HSP90AA1__chr14", "HSPA1A__chr6", "MTRNR2L1__chrX", "MTRNR2L12__chr3", "MTRNR2L8__chr11", "FOS__chr14", "MALAT1__chr11", "JUN__chr1", "DUSP1__chr5"),
                       .FGenes = c("HSPB1__chr7", "HSPA6__chr1", "HSPH1__chr13", "HSPA1B__chr6", "HSP90AB1__chr6", "FOSB__chr19", "HSP90B1__chr12", "RPS16__chr19","DNAJB1__chr19", "H3F3B__chr17", "HERPUD1__chr16", "NEAT1__chr11", "RP11-212I21.4__chr16", "IVNS1ABP__chr1", "HSPA8__chr11", "HIST1H2BG__chr6", "HSPA5__chr9", "JUNB__chr19","ZFP36L1__chr14"),
                       location_batch_info='/home/roman/Documents/Single cell analysis/20170921 Healthy microglia/20171210-batch-information.csv',
                       PCs_for_clustering = NULL,
                       .cln=cln
                       ) {
  require(tidyverse)
   
  data <- SCseq(.prdata)
  # filtering of expression data
  data <- filterdata(data, mintotal=.mintotal, 
                   minexpr=5, 
                   minnumber=1, 
                   maxexpr=Inf, 
                   downsample=FALSE, 
                   sfn=FALSE, 
                   hkn=FALSE, 
                   dsn=1, 
                   rseed=17000, 
                   CGenes=.CGenes,  
                   FGenes=.FGenes,
                   ccor=.4)
  
  # regress out the batch effect
  # optional:
  #build a table with cell Ids, cluster and conditions
  data_t <- as.data.frame(t(data@fdata))
  
  data_t$ID <- rownames(data_t)
  data_t$cell_ID <- rownames(data_t)
  data_t$Region_gm_wm <- ifelse(grepl('WM', data_t$ID), 'WM', 
                                ifelse(grepl('GM', data_t$ID), 'GM', 'Both')) 
  
  data_t$ID <- gsub('_.*', '', data_t$ID)
  data_t$ID <- gsub('-.*', '', data_t$ID)
  
  data_t$ID <- gsub('GM|WM|all|micr|pos|17Pl1|17Pl2', '', data_t$ID)
  
  batch_info <- read_csv(location_batch_info)[,-1]
  data_t <- merge(data_t, batch_info)
  table(data_t$ID, data_t$Batch)
  vars <- as.data.frame(data_t$Batch[data_t$cell_ID %in% colnames(data@fdata)])
  data@fdata <- varRegression(data@fdata,vars)
  
  # correct for cell cycle, proliferation, and expression of degradation markers by PCA
  # optional:
  x <- CCcorrect(data@fdata,
                 vset=NULL,
                 CGenes=.CGenes,
                 ccor=.4,
                 nComp=,
                 pvalue=.05,
                 quant=.01,
                 mode="pca")
  # number of principal components that have been removed
  x$n
  # loadings of the first principal component that has been removed
  y <- x$pca$rotation[,x$n[1]]
  # genes from vset are either enriched in the head or the tail of this list
  tail(y[order(y,decreasing=TRUE)],10)
  # reassign the corrected expression matrix to data@fdata
  data@fdata <- x$xcor
  
  # k-medoids clustering
  data <- clustexp(data,
                 clustnr=30,
                 bootnr=50,
                 metric="pearson",
                 do.gap=FALSE,
                 sat=TRUE,
                 SE.method="Tibs2001SEmax",
                 SE.factor=.25,
                 B.gap=50,
                 cln=.cln,
                 rseed=17000,
                 FUNcluster="kmedoids",
                 FSelect=TRUE)
  # compute t-SNE map
  data <- comptsne(data,
                 rseed=15555,
                 sammonmap=FALSE,
                 initial_cmd=TRUE,
                 fast=TRUE,
                 perplexity=30)
  # detect outliers and redefine clusters
  data <- findoutliers(data, 
                     outminc=5,
                     outlg=2,
                     probthr=1e-8,
                     thr=2**-(1:40),
                     outdistquant=.95)
  # reassign clusters based on random forest
  data <- rfcorrect(data,rfseed=12345,
                  final=TRUE,
                  nbfactor=5)
  
  data
}

#differential gene expression
differential_test <- function(data, pval.cutoff = 0.05, homeDirectory = home, subDirectory = "Cluster specific genes") {
  mainDir <- homeDirectory
  subDir <- subDirectory
  subsubDir_up <- 'Up'
  subsubDir_down <- 'Down'
  
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    dir.create(file.path(mainDir, subDir, subsubDir_up))
    dir.create(file.path(mainDir, subDir, subsubDir_down))
    setwd(file.path(mainDir, subDir))
    
  }
  for (i in c(1: max(data@cpart))) {
    cl <- names(data@cpart[data@cpart %in% c(i)])
    rest <- names(data@cpart[data@cpart %in% data@cpart[data@cpart != i]])
    diffexp <- diffexpnb(data@ndata,rest,cl,norm=F,vfit=data@background$vfit) 
    diffexpgenes <- diffexp[["res"]]
    diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < pval.cutoff)
    diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
    diffexpgenes$GENEID <- rownames(diffexpgenes)
    diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
    #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
    diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
    diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
    diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
    write.csv(diffgene_up, file = paste0('Up/', "diffgenes_cl", as.character(i), "_up_rest.csv"))
    write.csv(diffgene_down, file = paste0('Down/', "diffgenes_cl", as.character(i), "_down_rest.csv"))
    
    
    tryCatch({
      png(paste0('Up/', date,'-MA-plot-Cl', as.character(i) ,'_rest.png'), res = 300, width =7, height = 7, units = 'in')
      plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
      dev.off()
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    for (j in c(1: max(data@cpart))) {
      cl <- names(data@cpart[data@cpart %in% i])
      rest <- names(data@cpart[data@cpart %in% j])
      diffexp <- diffexpnb(data@ndata,rest,cl,norm=F,vfit=data@background$vfit) 
      diffexpgenes <- diffexp[["res"]]
      diffexpgenes <- subset(diffexpgenes, diffexpgenes$pval < pval.cutoff)
      diffexpgenes <- subset(diffexpgenes, abs(diffexpgenes$log2FoldChange) > 0)
      diffexpgenes$GENEID <- rownames(diffexpgenes)
      diffexpgenes <- diffexpgenes %>% dplyr::arrange(padj)
      #diffexpgenes$GENEID <- gsub('__.*', '', diffexpgenes$GENEID)
      diffexpgenes$Cluster <- rep(i, nrow(diffexpgenes))
      diffgene_up <- diffexpgenes[diffexpgenes$log2FoldChange > 0, ]
      diffgene_down <- diffexpgenes[diffexpgenes$log2FoldChange < 0, ]
      write.csv(diffgene_up, file = paste0('Up/', "diffgenes_cl", as.character(i), "_up_vs_cl", as.character(j), ".csv"))
      write.csv(diffgene_down, file = paste0('Down/', "diffgenes_cl", as.character(i), "_down_vs_cl", as.character(j), ".csv"))
      
      tryCatch({
        png(paste0('Up/', date,'-MA-plot-Cl', as.character(i) , "_vs_cl", as.character(j), ".png"), res = 300, width =7, height = 7, units = 'in')
        plotdiffgenesnb(diffexp,show_names = T, pthr=.01, lthr = .0001, mthr = .0001, padj=T)
        dev.off()
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      
    }
    
  }
  
}

#build a dataset
build_tidy_dataset <- function(data = sc, percent_cluster = 1, location_batch_info='/home/roman/Documents/Single cell analysis/20170921 Healthy microglia/20171210-batch-information.csv') {
  
  cell_numbers <-as.numeric()
  for (i in 1:length(unique(data@cpart)))
  {
    cell_numbers[i] <- length(data@cpart[data@cpart==i])
  }
  names(cell_numbers) <- c(1:length(unique(data@cpart)))
  retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(data@ndata)[2]/(100/percent_cluster)]))
  order_clusters <- clustheatmap(data,final=TRUE)
  order_clusters <- order_clusters[order_clusters %in% retain_cl]
  df <- data.frame('Cluster' = data@cpart, data@tsne) 
  df$ID <- rownames(df)
  df$cell_ID <- rownames(df)
  df$Condition <- gsub('_.*', '', rownames(df))
  df <- df[df$Cluster %in% retain_cl,] 
  df$Cluster <- factor(df$Cluster, levels = order_clusters)
  df$Cluster <- droplevels(df$Cluster)
  
  #add batch info
  batch_info <- read_csv(location_batch_info)[,-1]
  colnames(batch_info)[c(1,4)] <- c('Condition', 'Diagnosis')
  df <- left_join(df, batch_info)
  df$Condition <- as.factor(df$Condition)
  df$Condition <- reorder(df$Condition, df$Age)
  df$anon_ID <- df$Condition
  levels(df$anon_ID) <- paste0("Pat", 1:length(unique(df$Condition)))
  df$Region <- ifelse(grepl('WM',df$ID), 'WM', 
                      ifelse(grepl('GM',df$ID), 'GM', 'Mixed'))
  df$Region <- as.factor(df$Region)
  
  df
}

#enrichment test
hyper_test <- function(data1 = sc, data2 = df, var1 = "Cluster", var2 = "Region") {
  require(tidyverse)            
  clusters <- as_tibble(table(data1@cpart))
              colnames(clusters) <- c(var1, 'cluster_size')
              vars <- as_tibble(table(data2[,var1], data2[,var2]))
              colnames(vars) <- c(var1, var2, "freq_var2")
              vars_wide <- spread(vars, var2, freq_var2)
              
              vars_df <- vars_wide %>%
                left_join(clusters)
              
              
              #hypergeometric test
              #Male
              test_df<- data.frame(q=vars_df[,3], 
                                                   m=sum(vars_df[,3]), 
                                                   n=sum(vars_df[,2]),
                                                   k=vars_df[,4])
              
              p_hyper <- apply(test_df, MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
              padj <- p.adjust(p_hyper, method="BH")
              
              
              test_df$p_hyper <- p_hyper
              test_df$padj <- padj
              test_df$Cluster <- vars_df$Cluster
              test_df$Significance <- ifelse(test_df$padj<0.05, '*',
                                                            ifelse(test_df$padj<0.01, '**',
                                                                   ifelse(test_df$padj<0.001, '***','n.s.')))
              test_df$enrichment_var <- colnames(test_df)[1]
              
              #Female
              test_df2 <- data.frame(q=vars_df[,2], 
                                                     m=sum(vars_df[,2]), 
                                                     n=sum(vars_df[,3]),
                                                     k=vars_df[,4])
              p_hyper <- apply(test_df2, MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
              padj <- p.adjust(p_hyper, method="BH")
              
              test_df2$p_hyper <- p_hyper
              test_df2$padj <- padj
              test_df2$Cluster <- vars_df$Cluster
              test_df2$Significance <- ifelse(test_df2$padj<0.05, '*',
                                                              ifelse(test_df2$padj<0.01, '**',
                                                                     ifelse(test_df2$padj<0.001, '***','n.s.')))
              
              test_df2$enrichment_var <- colnames(test_df2)[1]
              
              return(bind_rows(test_df, test_df2))
}

#plots
plotexptsne2 <-  function(object,g,n="",logsc=FALSE){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
            if ( !is.numeric(logsc) & !is.logical(logsc) ) stop("argument logsc has to be logical (TRUE/FALSE)")
            if ( n == "" ) n <- g[1]
            l <- apply(object@ndata[g,names(object@ndata) %in% df$ID] - .1,2,sum) + .1
            if (logsc) {
              f <- l == 0
              l <- log2(l)
              l[f] <- NA
            }
            mi <- min(l,na.rm=TRUE)
            ma <- max(l,na.rm=TRUE)
            ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
            #ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            v <- round((l - mi)/(ma - mi)*99 + 1,0)
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            par(mar = c(3,5,2.5,2))
            data <- object@tsne[which(names(object@cpart) %in% df$ID),]
            plot(data,xlab="",ylab="",main=n,pch=19,cex=0,col="grey", xaxt='n', yaxt='n', bty="n") #, xaxt='n', yaxt='n', bty="n",ann=FALSE is to remove labels and axes, i also removed 'Dim 1' and 'Dim 2'
            kk <- order(v,decreasing=F)
            points(data[kk,1],data[kk,2],bg=ColorRamp[v[kk]],pch=21,cex=1.5, col=ColorRamp[v[kk]], lwd=0) #bg=ColorRamp[v[kk]], cex=1.75
            par(mar = c(3,2.5,2.5,2))
            image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="",
                  xaxt="n")
            layout(1)
          }

#

cluster_stacked_barplot <- function(data=df,FILL = Condition, fill_colors = colors_pat, fill_variable_name = "Condition") {
  #Stacked cluster plot
  cluster_stack_plot <- ggplot(data, aes(Cluster, fill = FILL, group = FILL)) + 
    geom_bar(position = 'fill', color = 'black', lwd = 0.25) +
    theme_minimal() +
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          text=element_text(size=17)) + 
    scale_fill_manual(fill_variable_name, values = fill_colors) +
    labs(title = 'Conditions in Cluster', y='#cells in Cluster / total #cells', x='Cluster')
  
  cluster_stack_plot
}


normalized_stack_plot <- function(data=df, Var1 = 'Cluster', Var2 = "Diagnosis", fill_colors = toupper(c('#f1a340','#998ec3')), line_width=0.25) {
  #normalized stacked map
  pats <- table(data[,Var1], data[,Var2])
  pats_prop <- as.data.frame(prop.table(pats, 2))
  
  ggplot(pats_prop, aes(Var1, Freq, fill = Var2)) +
    geom_bar(stat = 'identity' , position = "fill", color = 'black', lwd=line_width) +
    theme_minimal() +
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          text=element_text(size=17)) +
    scale_fill_manual(Var2, values = fill_colors) +
    labs(title = paste0(Var2,' in ', Var1), y='normalized % of cells', x=Var1)
}


#tsne plots

                tsne_plot <- function(data = df, FILL = Condition, fill_colors = colors_pat, point_outline = "black", point_size = 3, line_width = 0.25, point_shape=21, alpha_value = 1) {
                  tsne_plot <- ggplot(data, aes(V1, V2, fill = FILL, color = FILL)) +
                    geom_point(pch = point_shape, size = point_size, stroke = line_width, color = point_outline, alpha = alpha_value) +
                    theme(panel.background=element_blank(),
                          panel.border=element_blank(),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          plot.background=element_blank(),
                          text=element_text(size=17),
                          axis.text.x=element_blank(), 
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          legend.title = element_blank(),
                          legend.key = element_blank()) +
                    scale_fill_manual(values = fill_colors) +
                    scale_color_manual(values = fill_colors, guide = FALSE) 
                  
                  tsne_plot 
                }

#tsne plot without outline
              tsne_plot_no_outline <- function(data = df, FILL = Condition, fill_colors = colors_pat, point_outline = "black", point_size = 3, line_width = 0.25, point_shape=21, alpha_value = 1) {
                tsne_plot <- ggplot(data, aes(V1, V2, fill = FILL, color = FILL)) +
                  geom_point(pch = point_shape, size = point_size, stroke = line_width, alpha = alpha_value) +
                  theme(panel.background=element_blank(),
                        panel.border=element_blank(),
                        panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),
                        plot.background=element_blank(),
                        text=element_text(size=17),
                        axis.text.x=element_blank(), 
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        legend.title = element_blank(),
                        legend.key = element_blank()) +
                  scale_fill_manual(values = fill_colors) +
                  scale_color_manual(values = fill_colors, guide = FALSE) 
                
                tsne_plot 
              }

#Marimekko plot with stats
            mosaicGG <- function(data, X, FILL, rect_col = 'white', line_width = 0.25) {
              require(dplyr)
              require(reshape2)
              #require(ggthemes)
              # Proportions in raw data
              DF <- as.data.frame.matrix(table(data[[X]], data[[FILL]]))
              DF$groupSum <- rowSums(DF)
              DF$xmax <- cumsum(DF$groupSum)
              DF$xmin <- DF$xmax - DF$groupSum
              DF$X <- row.names(DF)
              DF$groupSum <- NULL
              DF_melted <- melt(DF, id = c("X", "xmin", "xmax"), variable.name = "FILL")
              DF_melted <- DF_melted %>%
                group_by(X) %>%
                mutate(ymax = cumsum(value/sum(value)),
                       ymin = ymax - value/sum(value))
              
              # Chi-sq test
              results <- chisq.test(table(data[[FILL]], data[[X]])) # fill and then x
              resid <- melt(results$residuals)
              names(resid) <- c("FILL", "X", "residual")
              
              # Merge data
              DF_all <- merge(DF_melted, resid)
              
              # Positions for labels
              DF_all$xposn <- DF_all$xmin + (DF_all$xmax - DF_all$xmin)/2
              index <- DF_all$xmax == max(DF_all$xmax)
              #DF_all$yposn <- DF_all$ymin[index] + (DF_all$ymax[index] - DF_all$ymin[index])/2
              yposn = 0
              # Plot
              g <- ggplot(DF_all, aes(ymin = ymin,  ymax = ymax, xmin = xmin,
                                      xmax = xmax, fill = residual)) +
                geom_rect(col = rect_col, lwd = line_width) +
                geom_text(aes(x = xposn, label = X),
                          y = 1, size = 3, angle = 90, hjust = 1, show.legend = FALSE,check_overlap = T) +
                geom_text(aes(x = max(xmax),  y = yposn, label = FILL),
                          size = 3, hjust = 1, show.legend = FALSE,check_overlap = T) +
                scale_fill_gradient2("Residuals") +
                scale_x_continuous(X, expand = c(0,0)) +
                scale_y_continuous("Proportion", expand = c(0,0)) +
                theme_minimal() +
                theme(legend.position = "bottom")
              print(g)
            }

#Marimekko plot without stats
            mosaicGG2 <- function(data, X, FILL, colors = colors_many, rect_col = 'white', line_width = 0.25) {
              require(dplyr)
              require(reshape2)
              #require(ggthemes)
              # Proportions in raw data
              DF <- as.data.frame.matrix(table(data[[X]], data[[FILL]]))
              DF$groupSum <- rowSums(DF)
              DF$xmax <- cumsum(DF$groupSum)
              DF$xmin <- DF$xmax - DF$groupSum
              DF$X <- row.names(DF)
              DF$groupSum <- NULL
              DF_melted <- melt(DF, id = c("X", "xmin", "xmax"), variable.name = "FILL")
              DF_melted <- DF_melted %>%
                group_by(X) %>%
                mutate(ymax = cumsum(value/sum(value)),
                       ymin = ymax - value/sum(value))
              
              # Chi-sq test
              results <- chisq.test(table(data[[FILL]], data[[X]])) # fill and then x
              resid <- reshape2::melt(results$residuals)
              names(resid) <- c("FILL", "X", "residual")
              
              # Merge data
              DF_all <- merge(DF_melted, resid)
              
              # Positions for labels
              DF_all$xposn <- DF_all$xmin + (DF_all$xmax - DF_all$xmin)/2
              index <- DF_all$xmax == max(DF_all$xmax)
              #DF_all$yposn <- DF_all$ymin[index] + (DF_all$ymax[index] - DF_all$ymin[index])/2
              yposn = 0
              # Plot
              g <- ggplot(DF_all, aes(ymin = ymin,  ymax = ymax, xmin = xmin,
                                      xmax = xmax, fill = FILL)) +
                geom_rect(col = rect_col, lwd = line_width) +
                geom_text(aes(x = xposn, label = X),
                          y = 1, size = 3, angle = 90, hjust = 1, show.legend = FALSE,check_overlap = T) +
                geom_text(aes(x = max(xmax),  y = yposn, label = FILL),
                          size = 3, hjust = 1, show.legend = FALSE,check_overlap = T) +
                scale_fill_manual(FILL, values = colors) +
                scale_x_continuous(X, expand = c(0,0)) +
                scale_y_continuous("Proportion", expand = c(0,0)) +
                theme_minimal() +
                theme(legend.position = "bottom")
              print(g)
            }
            

#GO terms enrichment analysis

enrichment_analysis <- function(data1 = sc, 
                                data2 = df_upgenes, 
                                ontology = "BP",
                                stress_genes = c("HSPA1A","MTRNR2L8","MTRNR2L12","HSP90AA1","MALAT1","ZFP36L1","ZFP36","MTRNR2L1","FOS","MALAT1","HSPB*","DUSP1","HSPH1","HSPA*","JUN","HSP90B1","RPS16","DNAJB1","H3F3B","HERPUD1","NEAT1","IVNS1ABP","HIST1H2BG","RP*","XIST")) 
  {
  #source("https://bioconductor.org/biocLite.R")
  #biocLite('org.Hs.eg.db')
  require(clusterProfiler)
  require(org.Hs.eg.db)
  require(tidyverse)
  keytypes(org.Hs.eg.db)
  require(pheatmap)
  
  data2 <- data2[!grepl(paste0("^(", paste(stress_genes, collapse='|'), ")"), data2$GENEID),  ]
  
  #load background genes
  back_genes <- rownames(data1@ndata)[which(apply(data1@ndata > .1, 1, sum)>0)]
  background <- bitr(sub('_.*', '',back_genes), fromType = "SYMBOL",
                     toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                     OrgDb = org.Hs.eg.db)
  
  #define empty data frame to collect data
  enrich_up <- data.frame(matrix(ncol = 10))
  colnames(enrich_up) <- c('ID','Description', 'GeneRatio', 'BgRatio' ,'pvalue', 'p.adjust', 'qvalue', 'geneID','Count' , 'Cluster')
  
  for (i in seq_along(unique(data2$Cluster)))  {
    
    tryCatch({
      gene <- gsub('_.*', '',data2$GENEID[data2$Cluster == i])
      gene.data2 <- bitr(gene, fromType = "SYMBOL",
                         toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                         OrgDb = org.Hs.eg.db)
     
      ggo <- groupGO(gene     = gene.data2[,3],
                     OrgDb    = org.Hs.eg.db,
                     ont      = ontology,
                     level    = 3,
                     readable = TRUE)
      
      
      ego <- enrichGO(gene          = gene.data2[,3],
                      universe      = background[,3],
                      OrgDb         = org.Hs.eg.db,
                      minGSSize     = 1,
                      ont           = ontology,
                      pool          = TRUE,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
      
      ego_simpl <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min, measure = "Wang")
      ego_simpl2 <- ego_simpl[!duplicated(ego_simpl@result$geneID)]
      ego_simpl2$Cluster <- rep(as.character(i, nrow(ego_simpl2)))
      enrich_up <- rbind(enrich_up, ego_simpl2)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  na.omit(enrich_up)
  
  
}

run_gsea <- function(data = df_upgenes, 
                     seed = 79106,
                     term2gene = "GO",
                     stress_genes = c("HSPA1A","MTRNR2L8","MTRNR2L12","HSP90AA1","MALAT1","ZFP36L1","ZFP36","MTRNR2L1","FOS","MALAT1","HSPB*","DUSP1","HSPH1","HSPA*","JUN","HSP90B1","RPS16","DNAJB1","H3F3B","HERPUD1","NEAT1","IVNS1ABP","HIST1H2BG","RP*","XIST")) {
  TERM2GENE <- read.gmt("/home/roman/Documents/Single cell analysis/msigdb.v6.2.symbols.gmt")
  
  if (term2gene == "GO") {
  TERM2GENE <- read.gmt("/home/roman/Documents/Single cell analysis/msigdb.v6.2.symbols.gmt")}
  if (term2gene == "immuno") {
  TERM2GENE <- read.gmt("/home/roman/Documents/Single cell analysis/c7.all.v6.2.symbols.gmt")
  }
  data <- data[!grepl(paste0("^(", paste(stress_genes, collapse='|'), ")"), data$GENEID),  ]
  df_gsea <- list()
  for (i in seq_along(unique(data$Cluster))) {
    tryCatch({
            data1 <- data %>%
            filter(Cluster == i) %>%
            group_by(GENEID) %>%
            filter(foldChange == max(foldChange))
            
            colnames(gene_symbols)[1] <- "GENEID"
            data1 <- left_join(data1, gene_symbols)
            
            gene_names <- data1$foldChange
            names(gene_names) <- data1$GENEID
            geneList <- sort(gene_names, decreasing = T)
            
            df_gsea1 <- GSEA(geneList=geneList, TERM2GENE = TERM2GENE)
            df_gsea[[i]] <- df_gsea1 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  df_gsea
}

#GO dot plot
go_dot_plot <- function(data = enrich_up, var1 = "Cluster", var2 = "Description", point_size = "GeneCount", FILL = "qvalue", line_width = 0.25, fill_name = "q-value",...) {
  require(tidyverse)
  require(viridis)
  dot_plot <- ggplot(enrich_up, aes_string(x=var1, y=var2, size = point_size, fill= -log10(data[[FILL]]))) + 
  geom_point(pch=21, stroke=line_width) +
  scale_fill_viridis(fill_name) +
  theme_light() +
  theme(text=element_text(size=10),
        axis.title.y=element_blank())
print(dot_plot)
}

#Histo cell counting
#count images in the target directory
          count_images <- function(source_path) {
            
            #counting existing images in dataset
            #count all the available images
            
            distance = data.frame()
            distance_summary_images = data.frame()
            nndistance = data.frame()
            delaunay = data.frame()
            derichelet = data.frame()
            all_coordinates = data.frame()
            delaunay_dist = data.frame()
            
            
            for (i in list.files(source_path)) {
              path=list.files(paste0(source_path, '/', as.character(i)))
              for (file in path) {
                for (tables in list.files(paste0(source_path, '/', as.character(i), '/', as.character(file)))) {
                  table_name = paste0(source_path, '/', as.character(i), '/', as.character(file), '/', as.character(tables))
                  
                  #generate table with summary stats
                  distance_summary2 = data.frame(ID = i, 
                                                 Condition = file, 
                                                 Image = tables
                  )
                  distance_summary_images = rbind(distance_summary_images, distance_summary2)
                  
                  
                }
                
                
              }
              
              
              
              
              
              
            }
            
            distance_summary_images$Image <- gsub('image', '', distance_summary_images$Image)
            distance_summary_images$Image <- gsub('.tif', '.csv', distance_summary_images$Image)
            
            return(distance_summary_images)
          }

          #pattern analysis
          pattern_analysis <- function(source_path, spatstat = FALSE) {
            #First version :20181015
            #last update :20181031
            #only usable for 20x magnification images
            
            #list.files()
            distance = data.frame()
            distance_summary = data.frame()
            delaunay = data.frame()
            derichelet = data.frame()
            all_coordinates = data.frame()
            delaunay_dist = data.frame()
            
            
            for (i in list.files(source_path)) {
              path=list.files(paste0(source_path, '/', as.character(i)))
              for (file in path) {
                tryCatch({
                  for (tables in list.files(paste0(source_path, '/', as.character(i), '/', as.character(file)))) {
                    if (grepl(".csv$", tables)) {table_name = paste0(source_path, '/', as.character(i), '/', as.character(file), '/', as.character(tables))
                    b = read.csv(table_name, 1, sep = ',', colClasses='numeric')[,6:7]
                    #calibrate the imagej output to microns
                    b$X = ((2560/8.51)*0.17) * b$X
                    b$Y = ((1920/6.39)*0.17) * b$Y
                    
                    #Generate data frame with all coordinates
                    coord = data.frame(ID = rep(i, nrow(b)), 
                                       Region = rep(file, nrow(b)), 
                                       Image = rep(tables, nrow(b)), 
                                       X = b$X,
                                       Y = b$Y)
                    all_coordinates = rbind(all_coordinates, coord)
                    
                    #Compute distances
                    distances = as.vector(dist(as.matrix(b)))
                    distances2 = data.frame(ID = rep(i, length(distances)), 
                                            Region = rep(file, length(distances)), 
                                            Image = rep(tables, length(distances)), 
                                            Distance = distances)
                    distance = rbind(distance, distances2) 
                    
                    #generate table with summary stats
                    distance_summary2 = data.frame(ID = i, 
                                                   Condition = file, 
                                                   Cell_count = nrow(b), 
                                                   Cell_density = nrow(b)/((2560*1920*(0.17)^2/1000000)), 
                                                   Mean_cell_dist = mean(distances), 
                                                   Median_cell_dist = median(distances), 
                                                   Range_cell_dist = max(distances)-min(distances),
                                                   Image = tables,
                                                   Antigen = strsplit(source_path, "/")[[1]][6])
                    distance_summary = rbind(distance_summary, distance_summary2)
                    
                    if (spatstat) { 
                      require(spatstat)
                      #Spatstat package anaysis
                      if (nrow(b) > 2) {
                        #Voronoi tessallation
                        derich = ppp(b$X,b$Y, xrange = range(0,(2560 * 0.17)), yrange = range(0,(1920 * 0.17)))
                        derich2  = dirichletAreas(derich)
                        #build data frame
                        derich3= data.frame(ID = rep(i, length(derich2)), 
                                            Region = rep(file, length(derich2)), 
                                            Image = rep(tables, length(derich2)), 
                                            Derichelet_tiles = derich2)
                        derichelet = rbind(derichelet, derich3) 
                        
                        if (nrow(b) > 3) {
                          #delaunay
                          delaun = delaunay(derich)
                          delaun2 = tiles(delaun)
                          delaun3 = unlist(lapply(delaun2, area.owin))
                          
                          #Build data frame
                          delaun4 = data.frame(ID = rep(i, length(delaun3)), 
                                               Region = rep(file, length(delaun3)), 
                                               Image = rep(tables, length(delaun3)), 
                                               Delaunay_tiles = delaun3)
                          delaunay = rbind(delaunay, delaun4)
                          
                          #delaunay distance
                          delaun10 = data.frame(ID = i,
                                                Region = file, 
                                                Image = tables,
                                                Modus = 'Patient', 
                                                Del_distance = mean(delaunayDistance(derich)))
                          delaunay_dist = rbind(delaunay_dist, delaun10)
                        }
                        
                        
                      }
                      
                    }
                    
                    }
                    
                  }
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                
              }
              
            }
            
            return(distance_summary)
            
          }
          
          
#make long dataset 
          make_data_long <- function(data1 = df, data2 = sc, order_clusters = levels(df$Cluster)) {
            data_t <- data.frame(t(data2@ndata), 'ID' = colnames(data2@ndata))
            
            data_long <- data_t %>% gather('Gene', 'Expression', -ID) %>%
              left_join(df[,c(1,4)])
            
            #order clusters based on hierarchical clustering
            data_long$Cluster <- factor(data_long$Cluster, levels = order_clusters)
            
            return(data_long)
          }
          
          
#line plot
          gene_line_plot <- function(data = data_long, gene) {
            line_plot <- ggplot(data[data$Gene %in% gene,], aes(x=ID, y=Expression, color = Cluster, fill = Cluster)) +
              geom_bar(stat = 'identity') + #width = 0.1, 
              facet_grid(facets = ~Cluster, 
                         drop = TRUE, 
                         #space = "free", 
                         scales = "free", 
                         switch = "x",
                         space = "free_x") +
              labs(title = gene, y = 'Gene Expression', x = 'Cluster') +
              theme_minimal() +
              theme(axis.line = element_blank(), 
                    #axis.title.y = element_blank(),
                    #axis.title.x = element_blank(),
                    axis.ticks.y = element_blank(), 
                    strip.text.x = element_text(),
                    #axis.text.y = element_text(size = 10), 
                    axis.text.x = element_blank(),#element_text(size = cex.col), 
                    strip.background = element_blank(), 
                    panel.grid = element_blank(),
                    panel.spacing.x = unit(0.2, units = 'line'),
                    legend.position = 'None') +
              scale_fill_manual(values =c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3')) +
              scale_color_manual(values = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3'))
            print(line_plot)
            
          }
          
          gene_line_plot_go <- function(data = go_term_exp, gene) {
            line_plot <- ggplot(data, aes(x=ID, y=Expression, color = Cluster, fill = Cluster)) +
              geom_bar(stat = 'identity') + #width = 0.1, 
              facet_grid(facets = ~Cluster, 
                         drop = TRUE, 
                         #space = "free", 
                         scales = "free", 
                         switch = "x",
                         space = "free_x") +
              labs(title = gene, y = 'Gene Expression', x = 'Cluster') +
              theme_minimal() +
              theme(axis.line = element_blank(), 
                    #axis.title.y = element_blank(),
                    #axis.title.x = element_blank(),
                    axis.ticks.y = element_blank(), 
                    strip.text.x = element_text(),
                    #axis.text.y = element_text(size = 10), 
                    axis.text.x = element_blank(),#element_text(size = cex.col), 
                    strip.background = element_blank(), 
                    panel.grid = element_blank(),
                    panel.spacing.x = unit(0.2, units = 'line'),
                    legend.position = 'None') +
              scale_fill_manual(values =c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3')) +
              scale_color_manual(values = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3'))
            print(line_plot)
            
          }
#generate cumulative gene expression values from enrich_up file
          
          go_term_gene_exp <- function(data1 = enrich_up, data2 = sc, data3 = df) {
            data_t <- data.frame(t(data2@ndata), 'ID' = colnames(data2@ndata))
            data <- data.frame(row.names = colnames(data2@ndata))
            data1 <- na.omit(data1[!duplicated(data1$Description),])
            
            for (i in 1:nrow(data1)) {
              tryCatch({
                n=data1[i,'Description']
                data_ <- dplyr::select(data_t, matches(paste0("^(", c(gsub('/', '|', data1[i,'geneID'])), ")")))
                data_ <- data.frame(rowSums(data_))
                colnames(data_) <- as.character(n)
                data <- cbind(data, data_)
              }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            }
            
            #add cluster and other data
            data$Cluster <- sc@cpart
            data$Cluster <- factor(data$Cluster, levels = levels(df$Cluster))
            data$ID <- rownames(data)
            data <- na.omit(data)
            data <- data[order(data$Cluster),]
            
            colnames(data) <- gsub("/", "_", colnames(data))
            return(data)
          }
          
          
          
#other functions
id2name <- function(x) sub("\\_\\_chr\\w+","",x)

name2id <- function(x,id) {
  ##  id[sub("\\_\\_chr\\w+","",id) %in% x]
  n <- c()
  for ( j in x ){ n <- append(n,id[grep(paste(j,"(\\_\\_chr\\w+|$|\\|)",sep=""),id)])
  }
  n
}

#loading genes
load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  do.call(rbind, tables)
}

#mean heatmap
mean_heatmap <- function(.sc = sc, 
                         .genes = genes, 
                         .retain_cl = retain_cl, 
                         .color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)) {
            mat <- as.matrix(.sc@ndata)
            nd <- mat[.genes,]
            nd[is.na(nd)] <- 0
            nd  <- t(nd[complete.cases(nd),])
            clust_n <- .sc@cpart[.sc@cpart %in% .retain_cl]
            mnd <- as.data.frame(cbind(clust_n,nd[rownames(nd) %in% names(clust_n),]))
            mean_mnd <- aggregate(mnd[, 2:dim(mnd)[2]], list(mnd$clust_n), mean)
            row.names(mean_mnd) <- paste("C",mean_mnd$Group.1,sep = "")
            mean_mnd$Group.1 <- NULL
            mean_mnd <- as.matrix(mean_mnd) + 0.1
            gene <- as.data.frame(t(log10(mean_mnd)))
            gene <- gene[complete.cases(gene),]
            row.names(gene) <- gsub('__.*', '', row.names(gene))
            pheatmap(gene, color = .color, cluster_cols = T, cluster_rows=T,fontsize_row = 8, border_color = F, show_rownames = T,show_colnames = T, scale = "row")
            
}


sc_heatmap <- function(.sc = sc, 
                       .genes = genes, 
                       .retain_cl = retain_cl, 
                       .color = colorRampPalette(c('blue', 'white', 'red'))(1000)) {
              exp <- as.data.frame(t(as.matrix(.sc@ndata[.genes,])+0.1))
              nam <- as.data.frame(.sc@cpart[.sc@cpart %in% .retain_cl])
              nam$newid <- paste(rownames(nam),nam$`.sc@cpart[.sc@cpart %in% .retain_cl]`,sep = "_cl")
              exp.new <- merge(exp, nam, by = "row.names")
              rownames(exp.new) <- exp.new$newid 
              exp.new <- exp.new[,c(3:ncol(exp.new)-1)] # here you have to increase 11 to one more number per gene you will add
              exp.new$`.sc@cpart[.sc@cpart %in% .retain_cl]` <- factor(exp.new$`.sc@cpart[.sc@cpart %in% .retain_cl]`, levels = ord_clust)
              exp.new <- exp.new[order(exp.new$`.sc@cpart[.sc@cpart %in% .retain_cl]`),]
              log10exp <- log10(exp.new[,-ncol(exp.new)])
              colnames(log10exp) <- gsub('_.*', '', colnames(log10exp))
              
              #find position of column breaks for heatmap (basically where the new cluster starts)
              breaks <- c(table(exp.new$`.sc@cpart[.sc@cpart %in% .retain_cl]`))
              breaks <- as.numeric(breaks)
              a <- c()
              for (i in 1:length(breaks)) {
                if (i==1) {
                  a <- breaks[1]
                } else {a <- c(a, a[i-1]+breaks[i])}
              }
              
              #define annotation colors from url: https://stackoverflow.com/questions/33292067/pheatmap-annotation-colors-and-border
              annotation_col = data.frame(ID = factor(exp.new$`.sc@cpart[.sc@cpart %in% .retain_cl]`))
              rownames(annotation_col)<-rownames(exp.new)
              cols <- c(colors_pat, colors_many)[1:length(.retain_cl)] #colorRampPalette(c(brewer.pal(length(breaks) ,name = 'Set3'), brewer.pal(length(breaks) - 10,name = 'Set2')))
              names(cols) <- unique(annotation_col$ID)
              annotation_colors <- list(ID = cols)
              
              
              pheatmap(t(log10exp),
                       show_colnames = F, 
                       color = .color,
                       cluster_cols = F, 
                       annotation_legend = T, 
                       cluster_rows = T, 
                       fontsize_row = 10, 
                       scale = 'column',
                       gaps_col=a[-length(a)],
                       annotation_col = annotation_col, 
                       annotation_colors = annotation_colors[1])
              
}

#load csv data
        load_data <- function(path) { 
          files <- dir(path, pattern = '\\.csv', full.names = TRUE)
          tables <- lapply(files, read.csv)
          do.call(rbind, tables)
        }

#plot_expmap
        plot_expmap <- function(gene, .sc=sc, point_size=2.5, logsc=FALSE, line_width=0.25, .retain_cl = retain_cl) {
          l <- colSums(as.data.frame(as.matrix(.sc@ndata))[rownames(.sc@ndata) %in% gene,] * min(.sc@counts)) + 0.1
          l <- l[names(.sc@cpart)[sc@cpart %in% .retain_cl]]
          mi <- min(l)
          ma <- max(l)
          ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
          ColorLevels <- seq(mi, ma, length = length(ColorRamp))
          v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
          
          kk <- bind_cols(data.frame('l'=l), .sc@tsne[which(.sc@cpart %in% .retain_cl),]) %>% arrange(l)
          
          if(logsc) {
            plot <- ggplot(kk, aes(V1, V2, fill = log(l))) +
              geom_point(size = point_size, pch = 21, stroke=line_width) +
              scale_fill_gradientn('', colors = ColorRamp) +
              theme_void() +
              labs(title = paste(gene, collapse = ',')) 
            return(plot)
          }
          else {
            plot <- ggplot(kk, aes(V1, V2, fill = l)) +
              geom_point(size = point_size, pch = 21, stroke=line_width) +
              scale_fill_gradientn('', colors = ColorRamp) +
              theme_void() +
              labs(title = paste(gene, collapse = ','))
            return(plot)
          }
          
        }
       
# hyper_test
        hyper_test_new <- function(data1 = sc, data2 = df, var1 = "Cluster", var2 = "Region") {
          require(tidyverse) 
          require(broom)
          clusters <- as_tibble(table(data1@cpart), .name_repair = 'unique')
          colnames(clusters) <- c(var1, 'cluster_size')
          vars <- as_tibble(table(data2[,var1], data2[,var2]), .name_repair = 'unique')
          colnames(vars) <- c(var1, var2, "freq_var2")
          vars_wide <- spread(vars, var2, freq_var2)
          
          vars_df <- vars_wide %>%
            left_join(clusters)
          
          
          #hypergeometric test
          #option a
          test_df<- data.frame(q=vars_df[,3], 
                               m=sum(vars_df[,3]), 
                               n=sum(vars_df[,2]),
                               k=vars_df[,4])
          
          p_hyper <- apply(test_df, MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
          padj <- p.adjust(p_hyper, method="BH")
          
          
          test_df$p_hyper <- p_hyper
          test_df$padj <- padj
          test_df$Cluster <- vars_df$Cluster
          test_df$Significance <- ifelse(test_df$padj<0.05 & test_df$padj>0.01, '*',
                                         ifelse(test_df$padj<0.01 & test_df$padj>0.001, '**',
                                                ifelse(test_df$padj<0.001, '***','n.s.')))
          test_df$enrichment_var <- colnames(test_df)[1]
          
          #option b
          test_df2 <- data.frame(q=vars_df[,2], 
                                 m=sum(vars_df[,2]), 
                                 n=sum(vars_df[,3]),
                                 k=vars_df[,4])
          p_hyper <- apply(test_df2, MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
          padj <- p.adjust(p_hyper, method="BH")
          
          test_df2$p_hyper <- p_hyper
          test_df2$padj <- padj
          test_df2$Cluster <- vars_df$Cluster
          test_df2$Significance <- ifelse(test_df2$padj<0.05 & test_df2$padj>0.01, '*',
                                          ifelse(test_df2$padj<0.01 & test_df2$padj>0.001, '**',
                                                 ifelse(test_df2$padj<0.001, '***','n.s.')))
          
          test_df2$enrichment_var <- colnames(test_df2)[1]
          
          return(bind_rows(test_df, test_df2))
        }
        
# hyper_test
        hyper_test_seurat <- function(data = df, var1 = "Cluster", var2 = "Region") {
          require(tidyverse) 
          require(broom)
          clusters <- as_tibble(table(data$Cluster), .name_repair = 'unique')
          colnames(clusters) <- c(var1, 'cluster_size')
          vars <- as_tibble(table(data[,var1], data[,var2]), .name_repair = 'unique')
          colnames(vars) <- c(var1, var2, "freq_var2")
          vars_wide <- spread(vars, var2, freq_var2)
          
          vars_df <- vars_wide %>%
            left_join(clusters)
          
          
          #hypergeometric test
          #option a
          test_df<- data.frame(q=vars_df[,3], 
                               m=sum(vars_df[,3]), 
                               n=sum(vars_df[,2]),
                               k=vars_df[,4])
          
          p_hyper <- apply(test_df, MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
          padj <- p.adjust(p_hyper, method="BH")
          
          
          test_df$p_hyper <- p_hyper
          test_df$padj <- padj
          test_df$Cluster <- vars_df$Cluster
          test_df$Significance <- ifelse(test_df$padj<0.05 & test_df$padj>0.01, '*',
                                         ifelse(test_df$padj<0.01 & test_df$padj>0.001, '**',
                                                ifelse(test_df$padj<0.001, '***','n.s.')))
          test_df$enrichment_var <- colnames(test_df)[1]
          
          #option b
          test_df2 <- data.frame(q=vars_df[,2], 
                                 m=sum(vars_df[,2]), 
                                 n=sum(vars_df[,3]),
                                 k=vars_df[,4])
          p_hyper <- apply(test_df2, MARGIN = 1, function(x) 1-phyper(x[[1]]-1, x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
          padj <- p.adjust(p_hyper, method="BH")
          
          test_df2$p_hyper <- p_hyper
          test_df2$padj <- padj
          test_df2$Cluster <- vars_df$Cluster
          test_df2$Significance <- ifelse(test_df2$padj<0.05 & test_df2$padj>0.01, '*',
                                          ifelse(test_df2$padj<0.01 & test_df2$padj>0.001, '**',
                                                 ifelse(test_df2$padj<0.001, '***','n.s.')))
          
          test_df2$enrichment_var <- colnames(test_df2)[1]
          
          return(bind_rows(test_df, test_df2))
        }
#GO term analysis
        go_term_analysis <- function(.df = df, ontogeny = "BP", .sc = sc, organism = 'org.Mm.eg.db') {
          require(clusterProfiler)
          require(organism,character.only = TRUE)
          require(tidyverse)
          #keytypes(ontogeny)
          require(viridis)
          require(pheatmap)
          require(RaceID)
          
          back_genes <- rownames(as.matrix(.sc@ndata))[which(apply(as.matrix(.sc@ndata) > 0, 1, sum)>0)]
          back_genes <- gsub('_.*', '', back_genes)
          
          background <- bitr(back_genes, fromType = "SYMBOL",
                             toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                             OrgDb = organism)
          background <- background[!duplicated(background$ENTREZID),]
          
          #define empty data frame to collect data
          enrich_up <- data.frame(matrix(ncol = 10))
          colnames(enrich_up) <- c('ID','Description', 'GeneRatio', 'BgRatio' ,'pvalue', 'p.adjust', 'qvalue', 'geneID','Count' , 'Cluster')
          
          for (i in unique(.df$Cluster))  {
            
            tryCatch({
              gene <- .df$GENEID[.df$Cluster == i]
              gene.df <- bitr(gene, fromType = "SYMBOL",
                              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                              OrgDb = organism)
              
              
              ggo <- groupGO(gene     = gene.df[,3],
                             OrgDb    = organism,
                             ont      = ontogeny,
                             level    = 3,
                             readable = TRUE)
              
              
              ego <- enrichGO(gene          = gene.df[,3],
                              universe      = background[,3],
                              OrgDb         = organism,
                              minGSSize     = 1,
                              ont           = ontogeny,
                              pool          = TRUE,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)
              
              
              ego_simpl <- my_simplify(res = ego,
                                       semData = godata(ont = ego@ontology))
              ego_simpl <- ego_simpl[!duplicated(ego_simpl$geneID),]
              
              ego_simpl$Cluster <- rep(i, nrow(ego_simpl))
              
              enrich_up <- rbind(enrich_up, ego_simpl)
              
              
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
          }
          
          return(na.omit(enrich_up))
          
        }        

#build index file        
        index_analysis <- function(id_left_half='Pat_GM_1_1_', id_right_half='Pat_GM_1_1_',path_index_file,.sc=sc, path_bcd='data/384-Well-plate-layout.cvs') {
          micr_names <- colnames(.sc@ndata)
          bcd_plate <- read_csv(path_bcd)
          
          bcd_plate$ID <- ifelse(bcd_plate$Y <13, paste0(id_left_half, bcd_plate$Barcode), paste0(id_right_half, bcd_plate$Barcode))
          tsne <- cbind(.sc@tsne, micr_names)
          colnames(tsne)[3] <- 'ID'
          bcd_plate <- bcd_plate %>% left_join(tsne)
          
          #load index data
          index <- data.frame()
          for (i in path_index_file) {
            index <- index %>% 
              bind_rows(read_csv(i)[,-c(1:3)])
          }
          
          colnames(index)[1:2] <- c('Y', 'X')
          
          #Merge datasets
          index_full1 <- bcd_plate %>% left_join(index)
          colnames(index_full1)
          
          #Adjust the colnames
          a <- colnames(index_full1)
          a <- gsub(".*? (.+)", "\\1", a)
          a <- gsub(".*? (.+)", "\\1", a)
          a <- gsub(".*? (.+)", "\\1", a)
          a <- gsub("-", "_", a)
          
          a
          colnames(index_full1) <- a
          index_full1[,7:ncol(index_full1)] <- apply(index_full1[,7:ncol(index_full1)], MARGIN = 2, scale)
          return(index_full1)
        }
        
        go_term_analysis_seurat <- function(.df = df, ontogeny = "BP", .sc = all, organism = 'org.Mm.eg.db') {
          require(clusterProfiler)
          require(organism,character.only = TRUE)
          require(tidyverse)
          require(viridis)
          require(pheatmap)
          require(RaceID)
          
          back_genes <- rownames(.sc@assays$RNA@counts)[which(rowMeans(as.matrix(.sc@assays$RNA@counts)) > 0)]
          back_genes <- gsub('_.*', '', back_genes)
          
          background <- bitr(back_genes, fromType = "SYMBOL",
                             toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                             OrgDb = organism)
          background <- background[!duplicated(background$ENTREZID),]
          
          #define empty data frame to collect data
          enrich_up <- data.frame(matrix(ncol = 10))
          colnames(enrich_up) <- c('ID','Description', 'GeneRatio', 'BgRatio' ,'pvalue', 'p.adjust', 'qvalue', 'geneID','Count' , 'Cluster')
          
          for (i in unique(.df$cluster))  {
            
            tryCatch({
              gene <- .df$gene[.df$cluster == i]
              gene.df <- bitr(gene, fromType = "SYMBOL",
                              toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                              OrgDb = organism)
              
              
              ggo <- groupGO(gene     = gene.df[,3],
                             OrgDb    = organism,
                             ont      = ontogeny,
                             level    = 3,
                             readable = TRUE)
              
              
              ego <- enrichGO(gene          = gene.df[,3],
                              universe      = background[,3],
                              OrgDb         = organism,
                              minGSSize     = 1,
                              ont           = ontogeny,
                              pool          = TRUE,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)
              
              
              ego_simpl <- my_simplify(res = ego,
                                       semData = godata(ont = ego@ontology))
              ego_simpl <- ego_simpl[!duplicated(ego_simpl$geneID),]
              
              ego_simpl$Cluster <- rep(i, nrow(ego_simpl))
              
              enrich_up <- rbind(enrich_up, ego_simpl)
              
              
            }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
          }
          
          return(na.omit(enrich_up))
          
        }
        
        
        my_simplify <- function(res=ego,
                                measure = 'Wang',
                                semData = godata(ont = ego@ontology),
                                by="pvalue",
                                cutoff=0.7,
                                select_fun=min) {
          
          require(GOSemSim)
          
          sim <- mgoSim(res$ID, res$ID,
                        semData = semData,
                        measure=measure,
                        combine=NULL)
          
          ## to satisfy codetools for calling gather
          go1 <- go2 <- similarity <- NULL
          
          sim.df <- as.data.frame(sim)
          sim.df$go1 <- row.names(sim.df)
          sim.df <- gather(sim.df, go2, similarity, -go1)
          
          sim.df <- sim.df[!is.na(sim.df$similarity),]
          
          ## feature 'by' is attached to 'go1'
          sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
          sim.df$go2 <- as.character(sim.df$go2)
          
          ID <- res$ID
          
          GO_to_remove <- character()
          for (i in seq_along(ID)) {
            ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
            ## if length(ii) == 1, then go1 == go2
            if (length(ii) < 2)
              next
            
            sim_subset <- sim.df[ii,]
            
            jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))
            
            ## sim.df <- sim.df[-ii[-jj]]
            GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
          }
          
          enrich_go <- res[!res$ID %in% GO_to_remove, ]
        }
        
#plot expression seurat
        plot_expmap_seurat <- function(features, object=all, reduction = "umap", dims=c(1,2), point_size=1, logsc=FALSE, line_width=0, .retain_cl = retain_cl) {
          
          dims <- paste0(Key(object = object[[reduction]]), dims)
          data <- FetchData(object = object, vars = c(dims, "ident", features),  slot = "data")
          
          if (ncol(data) > 4) {
          data2 <- data.frame(data[,1:3], rowSums(data[, 4:ncol(data)]))
          } else {
            data2 <- data
          }
          
          l <- data2[[4]][which(data$ident %in% .retain_cl)]
          mi <- min(l)
          ma <- max(l)
          ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
          ColorLevels <- seq(mi, ma, length = length(ColorRamp))
          v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
          
          kk <- bind_cols(data.frame('l'=l), data[, dims][which(data$ident %in% .retain_cl),]) %>% arrange(l)
          colnames(kk)[2:3] <- c("UMAP_1", "UMAP_2")
          
          if(logsc) {
            plot <- ggplot(kk, aes(UMAP_1, UMAP_2, color = log(l+0.1))) +
              geom_point(size = point_size, pch = 19) +
              scale_color_gradientn('', colors = ColorRamp) +
              theme_void() +
              labs(title = paste(features, collapse = ',')) 
            return(plot)
          }
          else {
            plot <- ggplot(kk, aes(UMAP_1, UMAP_2, color = l)) +
              geom_point(size = point_size, pch = 19) +
              scale_color_gradientn('', colors = ColorRamp) +
              theme_void() +
              labs(title = paste(features, collapse = ','))
            return(plot)
          }
          
        }
        
        umap_plot_seurat <- function(data = df, FILL = Condition, fill_colors = colors_pat, point_size = 0.5, line_width = 0.25, point_shape=19, alpha_value = 1) {
          tsne_plot <- ggplot(data, aes(UMAP_1, UMAP_2, fill = FILL, color = FILL)) +
            geom_point(pch = point_shape, size = point_size, alpha = alpha_value) +
            theme(panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank(),
                  text=element_text(size=17),
                  axis.text.x=element_blank(), 
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  legend.title = element_blank(),
                  legend.key = element_blank()) +
            #scale_fill_manual(values = fill_colors) +
            scale_color_manual(values = fill_colors) 
          
          tsne_plot 
        }
        
        plot_continuous <- function(gene, .df, point_size=2.5, logsc=FALSE) {
          l <- .df[[gene]]
          mi <- min(l)
          ma <- max(l)
          ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
          ColorLevels <- seq(mi, ma, length = length(ColorRamp))
          v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
          
          kk <- bind_cols(data.frame('l'=l), .df[, c("UMAP_1","UMAP_2")]) %>% arrange(l)
          colnames(kk)[2:3] <- c("V1", "V2")
          if(logsc) {
            plot <- ggplot(kk, aes(V1, V2, color = log(l))) +
              geom_point(size = point_size, pch = 19) +
              scale_color_gradientn('', colors = ColorRamp) +
              theme_void() +
              labs(title = paste(gene, collapse = ',')) 
            return(plot)
          }
          else {
            plot <- ggplot(kk, aes(V1, V2, color = l)) +
              geom_point(size = point_size, pch = 19) +
              scale_color_gradientn('', colors = ColorRamp) +
              theme_void() +
              labs(title = paste(gene, collapse = ','))
            return(plot)
          }
          
        }
        
        hyper_test_n <- function(data = df, var1 = "Cluster", var2 = "Region") {
          require(tidyverse) 
          require(broom)
          
          .df <- data.frame()
          for (i in unique(data[[var2]])) {
            data2 <- data
            data2[[var2]] <- factor(ifelse(data2[[var2]] == i, i, paste0("non_",i)), levels = c(i, paste0("non_",i)))
            clusters <- as_tibble(table(data2$Cluster), .name_repair = 'unique')
            colnames(clusters) <- c(var1, 'cluster_size')
            vars <- as_tibble(table(data2[,var1], data2[,var2]), .name_repair = 'unique')
            colnames(vars) <- c(var1, var2, "freq_var2")
            vars_wide <- spread(vars, var2, freq_var2)
            
            vars_df <- vars_wide %>%
              left_join(clusters)
            
            
            #hypergeometric test
            #option a
            test_df<- data.frame(q=vars_df[,i], 
                                 m=sum(vars_df[,i]), 
                                 n=sum(vars_df[,paste0("non_",i)]),
                                 k=vars_df[,4])
            
            colnames(test_df)[1] <- "q"
            
            p_hyper <- apply(test_df, MARGIN = 1, function(x) 1-phyper(max(0,x[[1]]-1), x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
            
            test_df$p_hyper <- p_hyper
            test_df$Cluster <- vars_df$Cluster
            test_df$enrichment_var <- i
            .df <- .df %>%
              bind_rows(test_df[,c("q","m","n","cluster_size","p_hyper","Cluster","enrichment_var")])
          }
          
        
        .df$padj <- p.adjust(.df$p_hyper, method="bonferroni")
        .df$Significance <- ifelse(.df$padj<0.05 & .df$padj>0.01, '*',
                                   ifelse(.df$padj<0.01 & .df$padj>0.001, '**',
                                          ifelse(.df$padj<0.001, '***','n.s.')))
        
          return(.df)
        }
        
        test_binom <- function(data = metadata, var1 = "Cluster", var2 = "Cellstate") {
          require(broom)
          .df <- data.frame()
          for (i in unique(data[[var2]])) {
            data2 <- data
            data2[[var2]] <- factor(ifelse(data2[[var2]] == i, i, paste0("non_",i)), levels = c(i, paste0("non_",i)))
            clusters <- as_tibble(table(data2$Cluster), .name_repair = 'unique')
            colnames(clusters) <- c(var1, 'cluster_size')
            vars <- as_tibble(table(data2[,var1], data2[,var2]), .name_repair = 'unique')
            colnames(vars) <- c(var1, var2, "freq_var2")
            vars_wide <- spread(vars, var2, freq_var2)
            
            vars_df <- vars_wide %>%
              left_join(clusters)
            
            
            #binomgeometric test
            #option a
            test_df<- data.frame(q=vars_df[,i], 
                                 m=sum(vars_df[,i]), 
                                 n=sum(vars_df[,paste0("non_",i)]),
                                 k=vars_df[,4])
            colnames(test_df)[1] <- "q"
            p_binom <- apply(test_df, MARGIN = 1, function(x) {glance(binom.test(x = max(0,x[[1]]-1), n = x[[4]], p = x[[2]]/(x[[2]] + x[[3]]), alternative = "greater"))[["p.value"]]}) #probability to get q or more successes in populaton
            
            test_df$p_binom <- p_binom
            test_df$Cluster <- vars_df$Cluster
            test_df$enrichment_var <- i
            .df <- .df %>%
              bind_rows(test_df[,c("q","m","n","cluster_size","p_binom","Cluster","enrichment_var")])
          }
          
          .df$padj <- p.adjust(.df$p_binom, method="bonferroni")
          .df$Significance <- ifelse(.df$padj<0.05 & .df$padj>0.01, '*',
                                     ifelse(.df$padj<0.01 & .df$padj>0.001, '**',
                                            ifelse(.df$padj<0.001, '***','n.s.')))
          
          
          
          return(.df)
        }
        
        plot_index <- function(gene, .index=index_all, point_size=5, log=T) {
          l <- .index[[gene]] + 0.1
          mi <- min(l, na.rm = T)
          ma <- max(l, na.rm = T)
          ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
          ColorLevels <- seq(mi, ma, length = length(ColorRamp))
          v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
          
          kk <- bind_cols(data.frame('l'=l), .index[,c('UMAP_1', 'UMAP_2')]) %>% arrange(l)
          
          if (log){ 
            plot <- ggplot(na.omit(kk), aes(UMAP_1, UMAP_2, color = log(l))) +
              geom_point(size = point_size, pch = 19, stroke=0.25) +
              scale_color_gradientn('', colors = ColorRamp) +
              theme_void() +
              labs(title = paste(gene, collapse = ','))}
          
          if (!log){ 
            plot <- ggplot(na.omit(kk), aes(UMAP_1, UMAP_2, color = l)) +
              geom_point(size = point_size, pch = 19, stroke=0.25) +
              scale_color_gradientn('', colors = ColorRamp) +
              theme_void() +
              labs(title = paste(gene, collapse = ','))}
          return(plot)
        }

                
#export cumulative data
        export_cumulative_expression_data <- function(features, object=all, reduction = "umap", dims=c(1,2), point_size=1, logsc=FALSE, line_width=0, .retain_cl = retain_cl) {
          
          dims <- paste0(Key(object = object[[reduction]]), dims)
          data <- FetchData(object = object, vars = c(dims, "ident", features),  slot = "data")
          
          if (ncol(data) > 4) {
            data2 <- data.frame(data[,1:3], rowSums(data[, 4:ncol(data)]))
          } else {
            data2 <- data
          }
          
          l <- data2[[4]][which(data$ident %in% .retain_cl)]
          names(l) <- rownames(data2)
          mi <- min(l)
          ma <- max(l)
          ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
          ColorLevels <- seq(mi, ma, length = length(ColorRamp))
          v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
          
          kk <- bind_cols(ID = names(l), data.frame('l'=l), data[, dims][which(data$ident %in% .retain_cl),]) %>% arrange(l)
          colnames(kk)[2:4] <- c("Expression","UMAP_1", "UMAP_2")
          
          kk
          
        }
        
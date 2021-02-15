# adjusted from Josip Herman
library(readr)
library(data.table)
library(tidyverse)
library(tools)
library(assertthat)

#this code was adjusted from Josip Herman, MPI Freiburg
# Specify directories
file_paths <- list(file.path("data","human_counts"))
names(file_paths) <- file_paths


# Search for all files in the specified directories and extract files by a given extension
files_list            <- lapply(file_paths, list.files, recursive=T)
files_by_ext          <- lapply(files_list, function(x){x[endsWith(x, suffix=".coutt.csv")]} )

# Get complete paths to all files
all_file_paths        <- unlist(lapply(seq_along(files_by_ext), function(x) {  file.path(names(files_by_ext[x]), files_by_ext[[x]]) } ))

names(all_file_paths) <- gsub("data/human_counts/|.coutt.csv", "", all_file_paths)

# Calculate md5sums to check for duplicated
md5sums      <- lapply(all_file_paths, function(x) {md5sum(x)} )

# Check for duplicated data
assert_that(sum(duplicated(md5sums)) == 0)

# Check for duplicated names
assert_that(sum(duplicated(unname(unlist(lapply(strsplit(unlist(files_by_ext),split = "/"),tail,1))))) == 0)

####
#### LOADING
####
# Loading data using lapply
data_list   <- lapply(all_file_paths, function(x) {fread(x, header= T)} )

# Add dataset name prefix to all columns, Merge with remaining gene names
for (d in names(data_list)) { 
  colnames(data_list[[d]]) <- c("GENEID", paste(d, "_",1:192,sep="" ))
  
}

# Cbind list of data.tables and removing the GENEID column from data.tables
data_list_cbind <- reduce(data_list, full_join, by = "GENEID")
data_list_cbind[is.na(data_list_cbind)]   <- 0

# Remove ERCC and mitochondrial genes
prdata <- as.data.frame(data_list_cbind)
rownames(prdata) <- prdata$GENEID
prdata$GENEID    <- NULL

#remove Kcnq1ot1 and correlated genes
cs <- colSums(prdata)
prdata <- prdata[,cs > 499]

#exclude not detected genes
prdata <- prdata[rowSums(prdata) >0,]

cs <- cs[cs > 499]
f <- t(prdata["KCNQ1OT1__chr11",])/cs < .02

#define sigcor
sigcor <- function(x,y,cthr=.4){
  if ( min(var(x),var(y)) == 0 ) return(NA)
  fit <- lm(x ~ y)
  pv <- as.data.frame(summary(fit)[4])[2,4]
  y <- as.data.frame(summary(fit)[4])[2,1]
  if ( is.na(pv) | is.na(y) ) return( NA )
  z <- sign(y)*sqrt(summary(fit)$r.square)
  if ( is.na(z) ) return(NA)
  if ( pv < .01 & abs(z) >= cthr ) return(z) else return(NA)
}


h <- rep(TRUE,nrow(prdata))
for ( g in "KCNQ1OT1__chr11"){
  z <- apply(prdata,1,function(x,y) sigcor(x,y),y=t(prdata[g,]))
  h <- h & ( is.na(z) | z < .65 )
}
prdata <- prdata[h,]

prdata <- prdata[!duplicated(gsub("_.*", "", rownames(prdata))),]
rownames(prdata) <- gsub("_.*", "", rownames(prdata))

save(prdata, file = "data/prdata.RData")

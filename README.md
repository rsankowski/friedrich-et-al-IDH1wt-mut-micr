# friedrich-et-al-IDH1wt-mut-micr

This is the code to reproduce the findings of the paper by Friedrich, Sankowski, Bunse et al. accepted in Nat Cancer. This is a step by step code book to get to the analyses in this study. It is roughly in the order of the appearance of the plots in the paper.

If you have any questions, please raise an issue or feel free to contact me anytime under: roman.sankowski@uniklinik-freiburg.de.

The structure of this project is fairly straighforward. The data folder contains the data. If you download the data for this project from GEO you should put it here. The plots folder contains the plot output and the R folder contains the code.

To get your system up to speed please run the code in the 0_setup.R script. To run the following analyses you will need to download the counts files from the GEO repository associated with this paper: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166420. 
- The mouse data is under the sub-series GSE166218,
- the human data under GSE166418 and 
- the mouse bulk data under GSE166521. 

Instead of re-running the seurat analyses you can directly download the associated objects:
- seurat-integration-wt-rh-gbm-final.RData for the human data and
- seurat-integration-10x-wt-rh-gbm.RData for the mouse data. 

All you need to do is download these files and put them in the data subfolder of this directory.


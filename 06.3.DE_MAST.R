
## Script header: take ident.1 and ident.2 from command line ----------------------
# Check if the number of command-line arguments is correct
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  cat("Usage: Rscript script_name.R arg1 arg2\n")
  quit(status = 1)
}

# Extract command-line arguments
ident.1 <- as.character(commandArgs(trailingOnly = TRUE)[1])
ident.2 <- as.character(commandArgs(trailingOnly = TRUE)[2])

# Print the extracted arguments (optional)
cat("Ident.1:", ident.1, "\n")
cat("Ident.2:", ident.2, "\n")


## Initiate analysis ---------------------------------------------------------
# use the following line to set the library path to /data/zhangy68/R/4.3;
# this path is writtable and used by R in sinteractive mode
# this line is nessacery if running interactive web-based Rstudio
.libPaths('/data/zhangy68/R/4.3')

packages_to_load = c(
  "dplyr",
  "Seurat",
  "patchwork",
  "ggplot2",
  "ggrepel",
  "readxl",
  "future",
  "cowplot",
  "Azimuth",
  "dsb",
  "corrplot",
  "DT",
  "tidyr",
  "ggpubr",
  "scRNAtoolVis",
  "purrr",
  "limma",
  "wesanderson",
  "reshape",
  "edgeR",
  'ggfun',
  'ggrepel',
  'fgsea'
)
# Load packages without printing output
invisible(suppressPackageStartupMessages({
  lapply(packages_to_load, library, character.only = TRUE)
}))

# set variables 
setwd('./')
set.seed('123')
col21 <- rev(c("tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque",
               "greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4",
               "orchid2","seagreen3","purple4","dodgerblue2","red","gray27"))



## Load Seurat object --------------
pbmc = readRDS('srt.obj.merged.clusterCuratedRefined.classifiedSet.myeloidInfoAdded.metadataIdentifierAdded.RDS')

## Load covariates table and add covariates to seurat metadata ---------------------
cov = read.csv('covariates.csv')
cov$sex = gsub(' ','',cov$sex)
pbmc@meta.data = merge(pbmc@meta.data,cov,
                       by.x = 'orig.ident',
                       by.y = 'sample',all.x = T)
rownames(pbmc@meta.data) = pbmc@meta.data$Row.names

## Run the test ------------------------------------
DefaultAssay(pbmc) <- 'RNA'
Idents(pbmc) <- 'curated_clusters4_group_trisomy'

print(paste0('now testing ', ident.1, ' and ', ident.2, ' using ', DefaultAssay(pbmc), ' matrix'))
print(paste0(ident.1, ': ', nrow(subset(pbmc@meta.data,curated_clusters4_group_trisomy == ident.1)), ' cells; ',
             ident.2, ': ', nrow(subset(pbmc@meta.data,curated_clusters4_group_trisomy == ident.2)), ' cells'))
if (nrow(subset(pbmc@meta.data,curated_clusters4_group_trisomy == ident.1)) < 50 || nrow(subset(pbmc@meta.data,curated_clusters4_group_trisomy == ident.2)) < 50) {
  print('comparison skipped due to low cell number (< 50 cells)')
  next # skip the comparison if either group has cell number < 50
}
res = FindMarkers(pbmc, ident.1 = ident.1, ident.2 = ident.2, verbose = T,
                  min.cells.group = 50,
                  test.use = 'MAST',latent.vars = 'sex', logfc.threshold = 0.1)

res$ident.1 = ident.1
res$ident.1.cellNumbers = nrow(subset(pbmc@meta.data,curated_clusters4_group_trisomy == ident.1))
res$ident.2 = ident.2
res$ident.2.cellNumbers = nrow(subset(pbmc@meta.data,curated_clusters4_group_trisomy == ident.2))

filename = gsub(' ','',paste0('DEG_MAST.',ident.1,'_VS_',ident.2,'.csv'))
write.csv(res,file = filename,quote = F,row.names = T)


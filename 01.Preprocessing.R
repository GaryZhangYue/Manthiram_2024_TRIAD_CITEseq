# NCBR-323 clean codes

#####################################################################
###################### seurat_qc.R ##################################
#####################################################################

## load the packages -------------
packages_to_load = c('Seurat','dplyr','patchwork','hdf5r','umap','cowplot','ggplot2','SeuratData','Routliers','scDblFinder','SingleCellExperiment')
lapply(packages_to_load,library, character.only = T)


## Set up the function for doublet finding
doublets <-function(dfso){
  set.seed(123)
  dfso <- as.SingleCellExperiment(dfso)
  dfso <- scDblFinder(dfso)
  dfso <- as.Seurat(dfso)
  return(dfso)
}

## load the h5 dataset -------------
samplelist = c('FYJ', 'GPG', 'HMM', 'JG', 'MMA', 'NYP', 'SRN', 'VKM', 'VR')

for (sample in samplelist) {
  print(sample)
  #h5.file = paste0('batch2/multi-2nd-pass_finalMultiUsingMultiDemultiplexingAssignment_downloaded/',sample,'.sample_filtered_feature_bc_matrix.h5')
  h5.file = paste0('batch2/multi-3rd-pass_finalMultiUsingHTODemuxAssignment_downloaded/',sample,'/count//sample_filtered_feature_bc_matrix.h5')
  print(h5.file)
  # read h5 file
  h5.obj = Read10X_h5(h5.file)
  # create seurat object for RNA data
  srt.obj = CreateSeuratObject(counts = h5.obj$`Gene Expression`,project = sample)
  # add ADT data as another assay
  srt.obj[['ADT']] <- CreateAssayObject(counts = h5.obj$`Antibody Capture`,project = sample)
  # validate the object
  print(Assays(srt.obj))
  # assign object
  assign(sample,srt.obj)
}

## Perform QC for each sample
srt.obj.list = list(FYJ,GPG,HMM,JG,MMA,NYP,SRN,VKM,VR)
n=1
for (srt.obj in srt.obj.list) {
  DefaultAssay(srt.obj) <- "RNA"
  sample=samplelist[n]
  # We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, 
  # which calculates the percentage of counts originating from a set of features
  srt.obj[["percent.mt"]] <- PercentageFeatureSet(srt.obj, pattern = "^MT-")
  
  # Visualize QC metrics as a violin plot -------------
  pv_untrimmed = VlnPlot(srt.obj, features = c("nFeature_RNA", "nCount_RNA", 'percent.mt'), ncol = 3, pt.size = 0)
  pf_untrimmed = FeatureScatter(srt.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  pf_untrimmed2 = FeatureScatter(srt.obj, feature1 = "percent.mt", feature2 = "nFeature_RNA")
  pf_untrimmed_t = pf_untrimmed+pf_untrimmed2
  # Filter the unwanted cells using outliers.mad -------------
  nCount_outliers_mad.3 = outliers_mad(log2(srt.obj$nCount_RNA),threshold = 3)
  nCount_outliers_mad.5 = outliers_mad(log2(srt.obj$nCount_RNA),threshold = 5)
  nCount_out  = c(nCount_outliers_mad.3$LL_CI_MAD,nCount_outliers_mad.5$UL_CI_MAD)
  
  nFeature_outliers_mad.3 = outliers_mad(log2(srt.obj$nFeature_RNA),threshold = 3)
  #nFeature_outliers_mad.5 = outliers_mad(log2(srt.obj$nFeature_RNA),threshold = 5)
  nFeature_out = c(nFeature_outliers_mad.3$LL_CI_MAD, nFeature_outliers_mad.3$UL_CI_MAD)
  mt_out = outliers_mad(log2(srt.obj$percent.mt),threshold = 3)$UL_CI_MAD 
  
  cellsToRemove.Feature= colnames(srt.obj)[which(log2(srt.obj$nFeature_RNA)<nFeature_out[1] | log2(srt.obj$nFeature_RNA)>nFeature_out[2])] 
  cellsToRemove.Count = colnames(srt.obj)[which(log2(srt.obj$nCount_RNA)<nCount_out[1])] 
  cellsToRemove.Mito = colnames(srt.obj)[which(log2(srt.obj$percent.mt)>mt_out)] 
  srt.obj <- subset(srt.obj,cells=unique(c(cellsToRemove.Feature,cellsToRemove.Count,cellsToRemove.Mito)),invert=T)
  # Find doublets ---------------
  # Use scblfinder
  srt.dblfinder <- doublets(srt.obj) 
  # print out singlet and doublet number
  print(table(srt.dblfinder$scDblFinder.class))
  # save singlet data for downstream analysis
  srt.obj <- subset(srt.dblfinder, scDblFinder.class=="singlet")
  # Make the same plots to examine the trimmed data
  pf_trimmed = FeatureScatter(srt.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",)
  pf_trimmed2 = FeatureScatter(srt.obj, feature1 = "percent.mt", feature2 = "nFeature_RNA")
  pv_trimmed = VlnPlot(srt.obj, features = c("nFeature_RNA", "nCount_RNA", 'percent.mt'), ncol = 3, pt.size = 0)
  pf_trimmed_t = pf_trimmed+pf_trimmed2
  
  pt_pv = plot_grid(pv_untrimmed,pv_trimmed, nrow = 2)
  save_plot(plot = pt_pv, 
            paste0('batch2/seurat_analysis/',sample,'_violinplot.png'), base_height = 15, base_width = 10)
  pt_pf = plot_grid(pf_untrimmed_t,pf_trimmed_t, nrow = 2)
  save_plot(plot = pt_pf, 
            paste0('batch2/seurat_analysis/',sample,'_featurescatters.png'), base_height = 15, base_width = 10)
  
  assign(sample,srt.obj)
  n=n+1
}

## save output for SCTransform and integration on Biowulf
srt.obj.list = list(FYJ,GPG,HMM,JG,MMA,NYP,SRN,VKM,VR)
saveRDS(srt.obj.list, file = 'batch2/seurat_analysis/srt.obj.list_QCed_beforeSCT.rds')

#####################################################################
####################### seurat_integration.R ########################
#####################################################################

## session 4 --------------------------------------------------------
# read input
srt.obj.list = readRDS(file = "srt.obj.list_QCed_beforeSCT.rds")
# merge all srt obj into one to ease the analysis
nr = NULL
sample_list = c(NULL)
for (i in seq(1,9,1)) {
  print(srt.obj.list[[i]]@meta.data$orig.ident[1])
  sample_list = c(sample_list,as.character(srt.obj.list[[i]]@meta.data$orig.ident[1]))
  
  print(DefaultAssay(srt.obj.list[[i]]))
  DefaultAssay(srt.obj.list[[i]]) = 'RNA'
  print(DefaultAssay(srt.obj.list[[i]]))
  nr = sum(nr,ncol(srt.obj.list[[i]]))
  print(ncol(srt.obj.list[[i]]))
}

pbmc = merge(x = srt.obj.list[[1]], y = c(srt.obj.list[[2]], 
                                          srt.obj.list[[3]], 
                                          srt.obj.list[[4]], 
                                          srt.obj.list[[5]], 
                                          srt.obj.list[[6]], 
                                          srt.obj.list[[7]],
                                          srt.obj.list[[8]],
                                          srt.obj.list[[9]]),
             add.cell.ids = sample_list,
             project = "merged")
nr == ncol(pbmc)# sanity check
# assign control/trisomy8 groups
pbmc@meta.data$group = ifelse(pbmc@meta.data$orig.ident %in% c('FYJ', 'SRN', 'VR'), 'control', 'trisomy8')
# sanity check
head(pbmc@meta.data)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", 'percent.mt'), ncol = 3, pt.size = 0)

## standard preprocessing steps
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) # find features (genes) expressed most variably between cells

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


# scale data to prepare data for dimensional reduction
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)

# cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP)
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap",split.by = "orig.ident", ncol = 3)
DimPlot(pbmc, reduction = "umap",group.by = "orig.ident")

pbmc  = RunAzimuth(pbmc, reference = 'pbmcref')
(plot.dimplot.azimuth.l1.byGroup = DimPlot(pbmc, reduction = "umap",group.by = 'predicted.celltype.l1',label = TRUE, label.size = 5, split.by ='group') + NoLegend())
(plot.dimplot.azimuth.l1.bySample = DimPlot(pbmc, reduction = "umap",group.by = 'predicted.celltype.l1',label = TRUE, label.size = 5,split.by = "orig.ident", ncol = 3) + NoLegend())

saveRDS(pbmc,'srt.obj.merged.analyzed.RDS')

## session 5 -------------------------------------------
# ADT data analysis

# read data
pbmc = readRDS(file = "srt.obj.merged.analyzed.RDS")
# rename ADT features (remove "-TotalSeqC" suffix)
pbmc@assays$ADT@counts@Dimnames[[1]] -> adt.old.names
adt.new.names= gsub('-TotalSeqC','',adt.old.names)
pbmc@assays$ADT@counts@Dimnames[[1]] <- adt.new.names
pbmc@assays$ADT@data@Dimnames[[1]] <- adt.new.names

DefaultAssay(pbmc) <- 'ADT'
rownames(pbmc)

# ADT Isotypes controls (mouse and rat isotypes)
adt.IsotypeControl.list = c("IgG1-1",
                            "IgG2a-1",
                            "IgG2b-1",
                            "IgG1-2",
                            "IgG2a-2",
                            "IgG2b-2",
                            "IgG")
# extract raw ADT count
adt = GetAssayData(pbmc, slot = 'counts', assay = 'ADT')
# normalize ADTs
dsb.norm = ModelNegativeADTnorm(cell_protein_matrix = adt,
                                denoise.counts = TRUE, 
                                use.isotype.control = TRUE, 
                                isotype.control.name.vec = adt.IsotypeControl.list
)

# plot a few proteins
theme_set(theme_bw())
plist = list(geom_vline(xintercept = 0, color = 'red'), 
             geom_hline(yintercept = 0, color = 'red'), 
             geom_point(size = 0.2, alpha = 0.1))
d = as.data.frame(t(dsb.norm))
adt.t = as.data.frame(t(as.matrix(adt)))
# plot distributions
p1 = ggplot(d, aes(x = `CD19`, y = `CD3`)) + plist
p2 = ggplot(adt.t, aes(x = `CD19`, y = `CD3`)) + geom_point()
p1n2 = cowplot::plot_grid(p1,p2)
save_plot('dsb_geom.png', p1n2, ncol = 2, base_asp = 1.1)

pdf(paste("dsb_histgram.pdf"), width = 12, height = 10)
par(mfrow = c(4,2))
hist(adt["CD4", ], breaks = 45, col = 'red', main = "CD4", xlab = 'Raw ADT')
hist(d[, "CD4"], breaks = 45, col = '#009ACD80', main = "CD4", xlab = 'ModelNegativeADTnorm')
hist(adt["CD8",], breaks = 45, col = 'red', main = "CD8", xlab = 'Raw ADT')
hist(d[, "CD8"], breaks = 45, col = '#009ACD80', main = "CD8", xlab = 'ModelNegativeADTnorm')
hist(adt["CD3",], breaks = 45, col = 'red', main = "CD3", xlab = 'Raw ADT')
hist(d[, "CD3"], breaks = 45, col = '#009ACD80', main = "CD3", xlab = 'ModelNegativeADTnorm')
hist(adt["CD88",], breaks = 45, col = 'red', main = "CD88", xlab = 'Raw ADT')
hist(d[, "CD88"], breaks = 45, col = '#009ACD80', main = "CD88", xlab = 'ModelNegativeADTnorm')
dev.off()

pdf(paste("dsb_histgram2.pdf"), width = 12, height = 10)
par(mfrow = c(5,4))
for (i in seq(1,10,1)) {
  print(i)
  hist(adt[i, ], breaks = 45, col = 'red', main = rownames(adt)[i], xlab = 'Raw ADT')
  hist(d[, i], breaks = 45, col = '#009ACD80', main = names(d)[i], xlab = 'ModelNegativeADTnorm')
}
dev.off()

# save the object without dsb as pbmc.beforeDsb
pbmc.beforeDsb = pbmc

## add the normalized data to the Seurat object
pbmc = SetAssayData(pbmc, slot = 'data', 
                    assay = 'ADT', 
                    new.data = dsb.norm)

## Centering and scaling data matrix
VariableFeatures(pbmc) <- rownames(pbmc[["ADT"]])
pbmc <- ScaleData(pbmc, assay = "ADT")

## PCA
pbmc <- RunPCA(pbmc, features = rownames(pbmc), reduction.name = "adt_pca", reduction.key = "adt_pca_", verbose = FALSE)
DimPlot(pbmc, reduction = "adt_pca")
ElbowPlot(pbmc,reduction = "adt_pca")

# Run non-linear dimensional reduction (UMAP)
pbmc <- RunUMAP(pbmc, dims = 1:20,reduction = 'adt_pca',reduction.name = 'adt_umap')
DimPlot(pbmc, reduction = "adt_umap",split.by = "orig.ident", label = TRUE, label.size = 5, repel = T, ncol = 3) + NoLegend()
DimPlot(pbmc, reduction = "adt_umap",split.by = "group", label = TRUE, label.size = 5, repel = T, ncol = 3) + NoLegend()
DimPlot(pbmc, reduction = "adt_umap",group.by = "predicted.celltype.l1", 
        label = TRUE, label.size = 5, repel = T, ncol = 3) + NoLegend()


# process ADT for WNN # see the main dsb vignette for an alternate version
# run WNN 
pbmc <- FindMultiModalNeighbors(
  pbmc, reduction.list = list("pca", "adt_pca"), 
  dims.list = list(1:20, 1:20), modality.weight.name = "RNA.weight"
)

saveRDS(pbmc,'srt.obj.merged.ADTdsb.RDS')
saveRDS(pbmc.beforeDsb,'srt.obj.merged.ADTbeforeDsb.RDS')

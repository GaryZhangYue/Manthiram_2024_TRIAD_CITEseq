# Add a few columns to the seurat metadata for DE comparisons
pbmc = readRDS('srt.obj.merged.clusterCuratedRefined.classifiedSet.myeloidInfoAdded.RDS')
pbmc@meta.data$curated_clusters4_trisomy = paste0(pbmc@meta.data$curated_clusters4, '-', pbmc@meta.data$trisomy)
pbmc@meta.data$curated_clusters4_group = paste0(pbmc@meta.data$curated_clusters4, '-',pbmc@meta.data$group)
pbmc@meta.data$curated_clusters4_group_trisomy = paste0(pbmc@meta.data$curated_clusters4_group, '-',pbmc@meta.data$trisomy)


saveRDS(pbmc, "srt.obj.merged.clusterCuratedRefined.classifiedSet.myeloidInfoAdded.metadataIdentifierAdded.RDS")
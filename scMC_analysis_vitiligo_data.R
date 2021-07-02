rm(list=ls())
library(Seurat)
library(cowplot)
library(dplyr)
library(scMC)
library(cowplot)

setwd("./Google Drive/projects/Vitiligo/2021")
counts_data <- read.table(file = "Vitiligo_processed_count_data.txt",sep = "\t")
counts_data <- as.matrix(counts_data)
meta_data <- read.table(file = "Vitiligo_processed_metadata.txt",sep = "\t")
conditions <- unique(meta_data$condition)
a <- unique(conditions)
patients <- meta_data$patient
b <- unique(patients)
data_list <- list()
meta_list <- list()
# as PatientG normal data only has 33 cell
for (i in 1:length(b)) {
  if (i < length(b)) {
    index1 <- which(patients == b[[i]] & conditions == a[[1]])
    data_list[[i]] <- counts_data[,index1]
    meta_list[[i]] <- meta_data[index1,]
  }

  index2 <- which(patients == b[[i]] & conditions == a[[2]])
  data_list[[i+(length(b)-1)]] <- counts_data[,index2]
  meta_list[[i+(length(b)-1)]] <- meta_data[index2,]
}
sample.name <- c(paste0(b[1:5],"_",a[[1]]), paste0(b,"_",a[[2]]))

######### Part I: Setup the a list of Seurat objects, one per dataset ##############
future::plan("multiprocess", workers = 4)
options(future.rng.onMisuse="ignore")
object.list <- vector("list", length(sample.name))
names(object.list) <- sample.name
for (i in 1:length(object.list)) {
  # Initialize the Seurat object with the raw (non-normalized data) for each dataset
  object.list0 <- CreateSeuratObject(counts = data_list[[i]], min.cells = 0, min.features = 0, meta.data = meta_list[[i]])
  # calculate mitochondrial QC metrics
  #object.list0[["percent.mt"]] <- PercentageFeatureSet(object.list0, pattern = "^mt-")
  #object.list0$sample.name <- sample.name[i]
  object.list[[i]] <-  object.list0
  rm(object.list0)
}

# step1. pre-processing each dataset
##  QC and selecting cells for further analysis
#for (i in 1:length(object.list)) {
# Visualize QC metrics as a violin plot
#  gg <- VlnPlot(object.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001,cols=c("#a6cee3"))
#  print(gg)
#  Sys.sleep(2)
# plot1 <- FeatureScatter(object.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1,cols=c("black")
# plot2 <- FeatureScatter(object.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1,cols=c("black"))
# plot1 + plot2
#  object.list[[i]] <- subset(object.list[[i]], subset = nFeature_RNA < 7000 & nCount_RNA < 40000 & percent.mt < 10)
#}

## Normalizing, scaling the data and feature selection
for (i in 1:length(object.list)) {
  x <- NormalizeData(object.list[[i]], verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "mean.var.plot", verbose = FALSE, mean.cutoff = c(0.01, 5), dispersion.cutoff = c(0.25, Inf))
  # perform scaling on the previously identified variable features
  x <- ScaleData(x, verbose = FALSE)
  object.list[[i]] <- x
}

########### Part II: Perform an integrated analysis using scMC ###########
# step2. identify clusters with different resolution for each condition
# compute SNN
object.list <- identifyNeighbors(object.list)
# identify clusters
object.list <- identifyClusters(object.list, mode = "separate",resolution = 1)

## step3. detect cluster-specific cells with high confident
features.integration = identifyIntegrationFeatures(object.list)
object.list <- identifyConfidentCells(object.list, features.integration,quantile.cutoff = 0.5)

## step4. Identify marker genes associated with the putative cell clusters in each dataset
object.list <- identifyMarkers(object.list, test.use = "bimod")


### check cluster of each data
features = c('KRT15','KRT5','KRT14','TYMS','TOP2A','KRT1','KRT2','FLG','S100A8','CXCL10','PMEL','PTPRC','CD207','CD3D')
for (i in 1:length(object.list)) {
  object1 <- RunUMAP(object.list[[i]],reduction='pca',dims = 1:40)
  gg <- featurePlot(object1, features = features, show.legend = F, show.axes = F)
  gg <- patchwork::wrap_plots(plots = gg, ncol = 4)
  cowplot::save_plot(filename=paste0("overlayKnownMarkers_integration_data_",i, "_umap.pdf"), plot=gg, base_width = 5.5, base_height = 6)

  p <- dimPlot(object1, reduction = "umap", label = TRUE)
  cowplot::save_plot(filename=paste0("data_",i, "_umap.pdf"), plot=p, base_width = 4, base_height = 3)
}


## step 5. Learn technical variation between any two datasets
structured_mat <- learnTechnicalVariation(object.list, features.integration, similarity.cutoff = 0.65)

## step 7. Learn a shared embedding of cells across all datasets after removing technical variation
combined <- merge(x = object.list[[1]],y = object.list[2:length(x = object.list)])
VariableFeatures(combined) <- features.integration
combined <- integrateData(combined, structured_mat, lambda = 1)

########### Part III: Run the standard workflow for visualization and clustering ###########
nPC = 40
combined <- FindNeighbors(combined, reduction = "scMC", dims = 1:nPC)
combined <- FindClusters(combined, algorithm = 4, resolution = 1)
levels(Idents(combined))
combined <- RunUMAP(combined, reduction='scMC', dims = 1:nPC)

combined <- ScaleData(combined, feature = rownames(combined), verbose = FALSE)


## Visualization
dimPlot(combined, reduction = "umap", label = TRUE)
dimPlot(combined, reduction = "umap", split.by = "conditions", combine = T)

p1 <- dimPlot(combined, reduction = "umap", group.by = "patient", colors.ggplot = T)
p2 <- dimPlot(combined, reduction = "umap", label = TRUE)
gg <- patchwork::wrap_plots(p1, p2, ncol = 2)
gg

# subclustering of melena and immune
combined.sub <- subset(combined, idents = c('12','14'))
combined.sub <- FindNeighbors(combined.sub, reduction = "scMC", dims = 1:nPC)
combined.sub <- FindClusters(combined.sub, algorithm = 4, resolution = 0.15)
levels(Idents(combined.sub))

dimPlot(combined.sub, reduction = "umap", label = TRUE)

gg <- featurePlot(combined.sub, features = features, show.legend = F, show.axes = T)
gg <- patchwork::wrap_plots(plots = gg, ncol = 4)
gg

## rename subclusters
## Annotation
## Labeling the clusters by cell type ##
# rename all
new.cluster.ids <- c("S2G 1","Spinous","B2S 2", "Basal 1", "Basal 2","Basal 1","Granular", "B2S 1", "Cycling","S2G 2", "Stress 2","Immune","Stress 1", "Immune", "Cycling")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)


new.cluster.ids <- c("DC", "Melanocytes", "TC")
names(new.cluster.ids) <- levels(combined.sub)
combined.sub <- RenameIdents(combined.sub, new.cluster.ids)

combined$cluster <- as.character(Idents(combined))
combined$cluster[Cells(combined.sub)] <- as.character(Idents(combined.sub))
Idents(combined) <- combined$cluster
dimPlot(combined, reduction = "umap", label = TRUE)


combined1 <- subset(combined.sub, idents = c('Melanocytes'))
combined1 <- FindNeighbors(combined1, reduction = "scMC", dims = 1:nPC)
combined1 <- FindClusters(combined1, algorithm = 4, resolution = 0.15)
levels(Idents(combined1))

dimPlot(combined1, reduction = "umap", label = TRUE)

features <- c('PMEL','MLANA','CD74','IL32','PTPRC','CD207','CD3D','CXCR3','CXCR6','CXCL16')
gg <- StackedVlnPlot(combined1, features = features)
gg

FeaturePlot(combined1, features = features)
new.cluster.ids1 <- c("Melanocytes","DC","DC")
names(new.cluster.ids1) <- levels(combined1)
combined1 <- RenameIdents(combined1, new.cluster.ids1)

combined.sub$subcluster <- as.character(Idents(combined.sub))
combined.sub$subcluster[Cells(combined1)] <- as.character(Idents(combined1))
Idents(combined.sub) <- combined.sub$subcluster
dimPlot(combined.sub, reduction = "umap", label = TRUE)

combined$cluster <- as.character(Idents(combined))
combined$cluster[Cells(combined.sub)] <- as.character(Idents(combined.sub))
Idents(combined) <- combined$cluster
dimPlot(combined, reduction = "umap", label = TRUE)


p1 <- dimPlot(combined, reduction = "umap", group.by = "conditions", colors.ggplot = T)
p2 <- dimPlot(combined, reduction = "umap", label = F)
gg <- patchwork::wrap_plots(p1, p2, ncol = 2)
gg


### reorder
new.cluster.ord <- c("Basal 1", "Basal 2","Cycling", "B2S 1", "B2S 2", "Spinous", "S2G 1", "S2G 2", "Granular",   "Stress 1", "Stress 2","Melanocytes", "DC","TC")
combined@active.ident <- factor(combined@active.ident, levels = new.cluster.ord)

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(combined, features = top10$gene) + NoLegend() + theme(axis.text.y = element_text(size = 10))



features = c('KRT5','KRT14','KRT15','CYR61','HES1','MT1G','TYMS','MCM7','TOP2A','CENPF','KRT1','KRT10','ATF3','KRT2','NOTCH3','SLURP1','FLG','LOR','CCL2','S100A8','S100A9','CXCL10','PMEL','MLANA','CD74','IL32','PTPRC','CD207','CD3D',"DCT", "CSF1R", "KRT6B","KRT6A","CD3E", "CSF1R", "CD207", "CLEC10A")
gg <- featurePlot(combined, features = features, show.legend = F, show.axes = F)
gg <- patchwork::wrap_plots(plots = gg, ncol = 4)
gg

color.use <- scPalette(nlevels(Idents(combined)))
gg <- StackedVlnPlot(combined, features = features, colors.use = color.use)
gg

gg <- StackedVlnPlot(combined, features = features, colors.use = color.use,split.by = "conditions")
gg


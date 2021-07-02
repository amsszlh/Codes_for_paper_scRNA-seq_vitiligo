# CellChat analysis of cell-cell communication between skin cells from nonlesional and lesional vitiligo
rm(list = ls())
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# load single-cell data for analysis
setwd("./Google Drive/projects/Vitiligo/2021")
library(Seurat)
combined <- readRDS(file = "Vitiligo_seurat.rds")
dataACDEFG = list(data = GetAssayData(combined, assay = "RNA"), identity = data.frame(group = Idents(combined), conditions = combined@meta.data$conditions, patient = combined@meta.data$patient))
data.input0 = dataACDEFG$data
identity0 = dataACDEFG$identity
levels(identity0$group)

############## CellChat analysis of nonlesional skin data #####################
cell.use = rownames(identity0)[identity0$conditions == "normal"]
data.input = data.input0[, cell.use]
identity = subset(identity0[cell.use,], select = 'group')

color.use.all <- c('#b2df8a',"#33a02c", '#03663b','#cab2d6','#6a3d9a','#ed98ba','#e14891','#d9423d','#9c2c5b','#a6cee3','#1f78b4',
                   '#b15928','#e3d322','#eca130')
color.use <- color.use.all[levels(identity$group) %in% unique(identity$group)]

meta = data.frame(labels = identity$group, row.names = names(identity$group))
## Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, do.sparse = F)
levels(cellchat@idents)

# start cellular communication analysis
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

## identify over-expressed ligands/receptors and pairs
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

## compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, idents = c(10,11)) # remove the stress cells in nonlesional skin due to the extremely low number of cells compared to other keratenocytes
# compute communications on signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
## calculate the aggregated network by counting the number of links or summarizing the communication probability
cellchat <- aggregateNet(cellchat)
# network importance analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

## Create a directory
data.dir <- './nonlesional'
dir.create(data.dir)
setwd(data.dir)
saveRDS(cellchat, file = "cellchat_patientALL_nonlesionalVitiligo.rds")


############## CellChat analysis of vitiligo lesional skin data #####################
cell.use = rownames(identity0)[identity0$conditions == "vitiligo"]
data.input = data.input0[, cell.use]
identity = subset(identity0[cell.use,], select = 'group')

color.use.all <- c('#b2df8a',"#33a02c", '#03663b','#cab2d6','#6a3d9a','#ed98ba','#e14891','#d9423d','#9c2c5b','#a6cee3','#1f78b4',
                   '#b15928','#e3d322','#eca130')
color.use <- color.use.all[levels(identity$group) %in% unique(identity$group)]

meta = data.frame(labels = identity$group, row.names = names(identity$group))
## Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, do.sparse = F)
levels(cellchat@idents)

# start cellular communication analysis
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

## identify over-expressed ligands/receptors and pairs
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

## compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
# compute communications on signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
## calculate the aggregated network by counting the number of links or summarizing the communication probability
cellchat <- aggregateNet(cellchat)
# network importance analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

## Create a directory
data.dir <- './lesional'
dir.create(data.dir)
setwd(data.dir)
saveRDS(cellchat, file = "cellchat_patientALL_lesionalVitiligo.rds")



############## Comparison analysis of cell-cell communication between vitiligo nonlesional and lesional skin #####################
## Create a directory to save figures
data.dir <- './Vitiligo/2021/comparison'
dir.create(data.dir)
setwd(data.dir)

# Load CellChat object of each dataset and then merge together
normal <- readRDS("./nonlesional/cellchat_patientALL_nonlesionalVitiligo.rds")
vitiligo <- readRDS("./lesional/cellchat_patientALL_lesionalVitiligo.rds")

object.list = list(normal,vitiligo)
object.list.names <- c("NL","LS")
names(object.list) <- object.list.names
colors.use.conditions <- c("#1b9e77","#7570b3")
colors.use <- c('#b2df8a',"#33a02c", '#03663b','#cab2d6','#6a3d9a','#ed98ba','#e14891','#d9423d','#9c2c5b','#a6cee3','#1f78b4',
                '#b15928','#e3d322','#eca130')

# merge CellChat objects
sapply(object.list, function(x){length(levels(x@idents))})
cellchat.merged <- mergeCellChat(object.list, add.names = object.list.names)
cellchat.merged
save(object.list, file = "cellchat_object.list_normal_vitiligo.RData")
save(cellchat.merged, file = "cellchat_merged_normal_vitiligo.RData")

## Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat.merged, color.use = colors.use.conditions,show.legend = T, x.lab.rot = F, group = rep(factor(c("NL","LS"), levels = c("NL","LS")), 1))
cowplot::save_plot("comparisonNumInteractions_normal_vitiligo.pdf", gg1, base_height = 5, base_width = 2.5)

## Identify and visualize the conserved and context-specific signaling pathways
###By comparing the information flow/interaction strengh of each signaling pathway, we can identify signaling pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase, by change their information flow at one condition as compared to another condition.
### Compare the overall information flow of each signaling pathway
gg <- list()
for (i in 1) {
  gg[[i]] <- rankNet(cellchat.merged, mode = "comparison", stacked = T, comparison = c(2*i-1,2*i), do.stat = T, color.use = colors.use.conditions) + theme(legend.position="top")
}
g <- patchwork::wrap_plots(plots = gg)
g
cowplot::save_plot(filename=paste0("comparison_contributions_signalingPathways_normal_vitiligo.pdf"), plot=g, base_width = 2.5, base_height = 5)

## Compare the major sources and targets in 2D space
# Comparing the outgoing and incoming interaction strength in 2D space allows ready identification of the cell populations with significant changes in sending or receiving signals between different datasets.
num.link <- lapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(sapply(num.link, function(x) min(x))), max(sapply(num.link, function(x) max(x)))) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], group = NULL, color.use = color.use.all,title = names(object.list)[i], weight.MinMax = weight.MinMax, show.legend = T) +
    xlim(c(0,1)) +ylim(c(0, 2.6))
  #  xlim(c(0,0.55)) +ylim(c(0, 1))
}
g <- patchwork::wrap_plots(plots = gg, byrow = FALSE)
g
cowplot::save_plot(filename=paste0("signalingRole_scatter_comparison_normal_vitiligo.pdf"), plot=g, base_width = 8, base_height = 3.5)


## Identify signaling changes associated with one cell group
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, cell.use = "Stress 2", signaling.exclude = "MIF", color.use = c("grey10", colors.use.conditions))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat.merged, cell.use = "Melanocytes", signaling.exclude = c("MIF","PTN"), color.use = c("grey10", colors.use.conditions))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat.merged, cell.use = "TC", signaling.exclude = c("MIF"), color.use = c("grey10", colors.use.conditions))
g <- patchwork::wrap_plots(plots = list(gg1,gg2,gg3), byrow = FALSE)
g
cowplot::save_plot(filename=paste0("signalingChanges_Stress2_Mel_TC_DC_scatter_comparison_normal_vitiligo.pdf"), plot=g, base_width = 9, base_height = 6)


#  Identify the upgulated and down-regulated signaling ligand-receptor pairs
levels(cellchat.merged@idents[[1]])
gg <- netVisual_bubble(cellchat.merged, sources.use = c(10, 11), targets.use = c(12:14), comparison = c(1, 2),  angle.x = 45, remove.isolate = F,title.name = "Signaling from stressed keratinocyte to melanocytes and immune cells")
gg
pairLR.use = unique(subset(gg$data, !is.na(interaction_name))[,"interaction_name", drop = FALSE])
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use, sources.use = c(10, 11), targets.use = c(12), comparison = c(1, 2),  angle.x = 45, remove.isolate = F,title.name = "Signaling from stressed keratinocyte to melanocytes and immune cells")
gg1
cowplot::save_plot("bubblePlot_Signaling_stressedKeratinToMelanocyteImmune.pdf", gg, base_height = 3, base_width = 5)
gg <- netVisual_bubble(cellchat.merged, targets.use = c(12), comparison = c(1, 2),max.dataset = 1, angle.x = 45, remove.isolate = F,title.name = "Decreased signaling to Melanocytes (NL vs. LS)", signaling = c("WNT","BMP","LIGHT"))
gg
cowplot::save_plot("bubblePlot_toMelanocytes_normal_vitiligo_decreasing.pdf", gg, base_height = 4, base_width = 5.5)


# Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
cairo_pdf(file = "circle_WNT_CXCL.pdf", width = 4, height = 4)
netVisual_aggregate(object.list[[1]], signaling = "WNT", layout = "circle", signaling.name = paste("WNT", names(object.list)[1]), color.use = color.use)
netVisual_aggregate(object.list[[2]], signaling = "CXCL", layout = "circle", signaling.name = paste("CXCL", names(object.list)[2]), color.use = color.use)
dev.off()


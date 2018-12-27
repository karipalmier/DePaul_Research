
library(Seurat)
library(dplyr)
library(Matrix)

dataPath = "C:\\DePaulCoursework\\Research\\Known_Examples\\Satija_Example_PBMC_4K\\hg19"

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = dataPath)

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
print(dense.size)

sparse.size <- object.size(x = pbmc.data)
print(sparse.size)

print(dense.size/sparse.size)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                           project = "10X_PBMC")


# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell. We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts. The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)


# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

# Normalize the data 
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Detection of variable genes across the single cells
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

print(length(x = pbmc@var.genes))

# Scaling the data and removing unwanted sources of variation
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

# Perform linear dimensional reduction
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = pbmc, pcs.use = 1:2)

PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)

PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)

JackStrawPlot(object = pbmc, PCs = 1:12)

PCElbowPlot(object = pbmc)

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)

# Run Non-linear dimensional reduction (tSNE)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc)

save(pbmc, file = "C:\\DePaulCoursework\\Research\\SatijaLabExamples\\pbmc_tutorial.Robj")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), 
                                min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)

VlnPlot(object = pbmc, features.plot = c("MS4A1", "CD79A"))

# you can plot raw UMI counts as well
VlnPlot(object = pbmc, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)

FeaturePlot(object = pbmc, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", 
                                             "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

#################################################################################################################
# Kari Palmier 2/7/2018 - Added this code to extract the heatmap data then print it in JSON format to a text file

# Get the path of the current based on if RStudio is run or if run from command line
cmdArgs <- commandArgs(trailingOnly = FALSE)

# If run from RStudio, use the rstudio api to get the path of the current R script
if (cmdArgs[1] == "RStudio"){
  mainPath = dirname(rstudioapi::getSourceEditorContext()$path)
  mainPath = gsub("/", "\\\\", mainPath)
  mainPath = paste(mainPath, '\\', sep="")
  
# Else if running from command line or R console, extract path info from command arguments
} else {
  
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  
  # If running from command line using RScript command
  if (length(match) > 0) {
    # Rscript
    fullPath = normalizePath(sub(needle, "", cmdArgs[match]))

  # Else running from R Console
  } else {
    # 'source'd via R console
    fullPath = normalizePath(sys.frames()[[1]]$ofile)
  }
  
  currentRFile = 'SajitaLabExample.R'
  mainPath = gsub(currentRFile, "", fullPath)
}

# Sourse the GetHeatmapData function using the full path (with directory above)
# *** Note that this assumes the GetHeatmapData.R file is in the same directory as the current R script
source(paste(mainPath, "GetHeatmapData.R", sep=""))

# Call GetHeatmapData.  The output will be JSON.
jsonData = GetHeatmapData(object = pbmc, genes.use = top10$gene)

# Write the output JSON variable contents to a text file
# *** Note 1: The path will need to be changed to match where the system should actually write the file
# *** Note 2: To write to standard output, just run the print command (remove all file interface & sink commands)
outFileDir = gsub("Code", "R_Output", mainPath)
outFileName = paste(outFileDir, "HeatmapData.txt", sep="")

outFile = file(outFileName, open="wt")
sink(file=outFile, append=TRUE)

print(jsonData)

sink()
close(outFile)
closeAllConnections()

#################################################################################################################



#current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
#new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
#                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
#pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
#TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)

# First lets stash our identities for later
#pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the
# SNN, and can simply put: pbmc <- FindClusters(pbmc,resolution = 0.8)
#pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
#                     resolution = 0.8, print.output = FALSE)

# Demonstration of how to plot two tSNE plots side by side, and how to color
# points based on different criteria
#plot1 <- TSNEPlot(object = pbmc, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
#plot2 <- TSNEPlot(object = pbmc, do.return = TRUE, group.by = "ClusterNames_0.6", 
#                  no.legend = TRUE, do.label = TRUE)
#plot_grid(plot1, plot2)

# Find discriminating markers
#tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)

# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we
# can see that CCR7 is upregulated in C0, strongly indicating that we can
# differentiate memory from naive CD4 cells.  cols.use demarcates the color
# palette from low to high expression
#FeaturePlot(object = pbmc, features.plot = c("S100A4", "CCR7"), cols.use = c("green", 
#                                                                             "blue"))
#pbmc <- SetAllIdent(object = pbmc, id = "ClusterNames_0.6")
#save(pbmc, file = "C:\\DePaulCoursework\\Research\\SatijaLabExamples\\pbmc3k_final.Rda")



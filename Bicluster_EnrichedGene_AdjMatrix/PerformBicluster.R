#####################################################################################################################
#
#  PerformBicluster.R
#
#  This script imports and filters data for mitochondrial and number of genes, then perform PCA, 
#  Seurat biclustering (using PCA components), and finally extract markers based on differential expressions.
#
#  Refer the the tutorial at https://satijalab.org/seurat/pbmc3k_tutorial.html for details on Seurat functions
#
#  Input Arguments:
#     filename = data file path  
#     colnames_row_1 = column names present (TRUE/FALSE) 
#     rownames_col_1 = row names present (TRUE/FALSE)
#     cluster_res = biclustering resolution to use
#     cluster_dims_used = number of PCA dimensions to use for biclustering
#     markers_min_pct = markers minimum percent to use in differential expression 
#     markers_test_use = markers test to use in differential expression
#     markers_thresh_use = markers threshold to use in differential expression
#     ngene_up_limit = number of gene upper limit for fitering
#     ngene_low_limit = number of gene lower limit for fitering
#     perc_mito_up_limit = percent mitochondrial upper limit for filtering
#     perc_mito_low_limit = percent mitochondrial lower limit for filtering
#     jpg_path = path to save heatmap image - if is empty, image is not saved
#     filter_colsums_5000 = flag to filter out any cell columns that sum to less than 5000
#
#  Outputs:
#     JSON string of output_data R list (contents below):
#        gene_cluster_data: list of following
#                             gene = list of genes from clustering 
#                             cluster = cluster assignments for each gene
#        heatmap_data: matrix of gene symbols (rows), cell id (columns), and expression levels (cell values)
#
#  Date         Name          Comment
#  7/2/2018    Kari Palmier   File Created
#  7/11/2018   Kari Palmier   File updated to use command arguments, to run through Python, to save off heatmap,
#                             and to pass out both heatmap and gene/cluster info
#  8/29/2018   Kari Palmier   File commented
#  9/20/2018   Kari Palmier   Added data filtering limit script inputs
#  1/9/2019    Kari Palmier   Removed number of PC components from inputs - hard-coded to 20
#                             This parameter had no effect on clustering results
#
#
#####################################################################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(jsonlite)

cmdArgs = commandArgs(trailingOnly=FALSE)
needle <- "--file="
match <- grep(needle, cmdArgs)
fullPath = normalizePath(sub(needle, "", cmdArgs[match]))
currentRFile = 'PerformBicluster.R'
mainPath = gsub(currentRFile, "", fullPath)

# mainPath = "C:\\DePaulCoursework\\Research\\Code\\"
   
source(paste(mainPath, "GetHeatmapListData.R", sep=""))

num_args_expected = 14

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
num_args_in = length(args)

if (num_args_in < num_args_expected){
 msg = paste("Missing input arguments. Expected = ", num_args_expected, ". Received = ", num_args_in, ".", sep = "")
 stop(msg)
}

# Get inputs
filename = args[1]
colnames_row_1 = args[2]
rownames_col_1 = args[3]
cluster_res = args[4]
cluster_dims_used = args[5]
markers_min_pct = args[6]
markers_test_use = args[7]
markers_thresh_use = args[8]
ngene_up_limit = args[9]
ngene_low_limit = args[10]
perc_mito_up_limit= args[11]
perc_mito_low_limit = args[12]
jpg_path = args[13]
filter_colsums_5000 = args[14]

# Test inputs
# filename =
#   'C:\\DePaulCoursework\\Research\\PBMC_Examples\\Satija_Example_PBMC_4K\\hg19.csv'
# colnames_row_1 = 'TRUE'
# rownames_col_1 = 'TRUE'
# cluster_res = '0.5'
# cluster_dims_used = '10'
# markers_min_pct = '0.25'
# markers_test_use = 'roc'
# markers_thresh_use = '0.25'
# ngene_up_limit = '5000'
# ngene_low_limit = '1'
# perc_mito_up_limit = '0.04'
# perc_mito_low_limit = '-Inf'
# jpg_path =
#   'C:\\DePaulCoursework\\Research\\PBMC_Examples\\Satija_Example_PBMC_4K\\201896_18_8_6_Cluster_Results\\'
# filter_colsums_5000 = 'FALSE'

# Convert inputs
cluster_res = as.numeric(cluster_res)
cluster_dims_used = as.numeric(cluster_dims_used)
markers_min_pct = as.numeric(markers_min_pct)
markers_thresh_use = as.numeric(markers_thresh_use)
ngene_up_limit = as.numeric(ngene_up_limit)
ngene_low_limit = as.numeric(ngene_low_limit)
perc_mito_up_limit = as.numeric(perc_mito_up_limit)
perc_mito_low_limit = as.numeric(perc_mito_low_limit)

# Verify csv colnames_row_1 parameter is valid
if ((colnames_row_1 != "TRUE") && (colnames_row_1 != "FALSE")){
  msg = "Incorrect is_csv parameter value."
  stop(msg)
}

# Verify csv rownames_row_1 parameter is valid
if ((rownames_col_1 != "TRUE") && (rownames_col_1 != "FALSE")){
  msg = "Incorrect first_row_names parameter value."
  stop(msg)
}

# Verify csv markers_test_use parameter is valid
if ((markers_test_use != "roc") && (markers_test_use != "bimodal") && (markers_test_use != "wilcox")){
  msg = "Incorrect markers_test_use parameter value."
  stop(msg)
}

# Read input data based on whether there are row and column titles
if ((colnames_row_1 == "TRUE") && (rownames_col_1 == "TRUE")){
  input_data = read.table(filename, sep=",", header=TRUE, row.names=1)
} else if (colnames_row_1 == "TRUE"){
  input_data = read.table(filename, sep=",", header=TRUE)
} else if (rownames_col_1 == "TRUE"){
  input_data = read.table(filename, sep=",", row.names=1)
} else{
  input_data = read.table(filename, sep=",")
}

# If filter_colsums_5000 is TRUE, filter out columns with sums less than 4999
if (filter_colsums_5000 == "TRUE"){
  total_filtered_data<- input_data[, colSums(input_data) > 4999]
  total_filtered_data<- total_filtered_data[rowSums(total_filtered_data) > 0,]
} else{
  total_filtered_data = input_data
}

# Set number of PCA components to use
pca_num_pcs = 20

# convert input data to Seurat object
pbmci <- CreateSeuratObject(total_filtered_data, min.cells = 1, min.genes = 1)

# Get mitochondrial information
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmci@data), value = TRUE)
percent.mito <- colSums(pbmci@raw.data[mito.genes, ]) / colSums(pbmci@raw.data)

# Add mitochondrial information to Seurat object
pbmci <- AddMetaData(pbmci, percent.mito, "percent.mito")
pbmc <- AddMetaData(object = pbmci, metadata = percent.mito, col.name = "percent.mito")

# Filter cells by number of genes and percent mitochondrial
pbmci <- FilterCells(object = pbmci, subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(ngene_low_limit, perc_mito_low_limit), 
                    high.thresholds = c(ngene_up_limit, perc_mito_up_limit))

# Normalize the data by total expression, then multiply by 10,000 and log-transform result
pbmci <- NormalizeData(object = pbmci, normalization.method = "LogNormalize",
                      scale.factor = 10000, display.progress = FALSE)

# Detect highly variable genes for later usage
pbmci <- FindVariableGenes(object = pbmci, mean.function = ExpMean, dispersion.function = LogVMR,
                           x.low.cutoff = 0.0125, y.cutoff = 0.0125, display.progress = FALSE, do.plot = FALSE)

# Remove unwanted sources of variation through regression undesired signals out
pbmci <- ScaleData(object = pbmci, vars.to.regress = c("nUMI", "percent.mito"), display.progress = FALSE)

# Perform primary componenet analysis linear dimension reduction
pbmci <- RunPCA(object = pbmci, pc.genes = pbmci@var.genes, do.print = FALSE, pcs.compute = pca_num_pcs)

# Perform biclustering
pbmci <- FindClusters(object = pbmci, reduction.type = "pca", dims.use = 1:cluster_dims_used,
                     resolution = cluster_res, print.output = FALSE, save.SNN = TRUE, force.recalc = TRUE)

# Perform non-linear dimensional reduction
pbmci <- RunTSNE(object = pbmci, dims.use = 1:cluster_dims_used, do.fast = TRUE)

# Find differentially expressed genes (cluster biomarkers)
pbmciroc_beta.markers <- FindAllMarkers(pbmci, only.pos=TRUE, min.pct = markers_min_pct, 
                                        thresh.use = markers_thresh_use, test.use = markers_test_use, 
                                        do.print = FALSE, print.bar = FALSE)

# Create heatmap file name and save if jpg path was provided
if (!(jpg_path == "")){
  date_str = format(Sys.time(), format="%Y%m%d_%H_%M_%s")
  savename = paste(jpg_path, date_str, "_res_", cluster_res, "_test_", markers_test_use, "_ngene_",  
                   ngene_up_limit, "_", ngene_low_limit, "_percmito_",  perc_mito_up_limit, "_",   
                   perc_mito_low_limit, "_cdims_", cluster_dims_used, "_Heatmap.jpg", 
                   sep = '')
  jpeg(file = savename)
}

DoHeatmap(object = pbmci, genes.use = pbmciroc_beta.markers$gene, slim.col.label = TRUE,
          remove.key = TRUE)

if (!(jpg_path == "")){
  garbage = dev.off()
}

# Extract heatmap gene/cell/expression data
heatmap_data = GetHeatmapListData(object = pbmci, genes.use = pbmciroc_beta.markers$gene)

# Create list of gene name and cluster assignments
gene_cluster_data = list(gene = pbmciroc_beta.markers$gene,
                         cluster = pbmciroc_beta.markers$cluster)

# Create combined list for all outputs
output_data = list(heatmap_data = heatmap_data, gene_cluster_data = gene_cluster_data)

# Convert combined output list to JSON to pass to Python
jsonData = jsonlite::toJSON(output_data, pretty = TRUE)

# Print combined output JSON to Python
print(jsonData)
 

#####################################################################################################################
#
#  Convert10xToCSV.R
#
#  This script converts data that has been saved into the 10X format into a standard 2D dataframe, then saves
#  the dataframe to a CSV file.
#
#
#  Date         Name          Comment
#  9/6/2018     Kari Palmier  File Created
#
#####################################################################################################################

library(Seurat)
library(dplyr)
library(Matrix)

#data_path = "C:\\DePaulCoursework\\Research\\SatijaLabExamples\\filtered_gene_bc_matrices\\hg19"
data_path = choose.dir(default = "", "Select Folder Containing 10X Data:")

# Formulate starting path for output file saving based on input file path
search_results = gregexpr("\\\\", data_path)
num_slashes = length(search_results[[1]])
slash_ndx = search_results[[1]][num_slashes]
start_ndx = slash_ndx + 1
last_ndx = nchar(data_path)
out_name = substr(data_path, start_ndx, last_ndx)
out_base_path = substr(data_path, 1, slash_ndx)
min_genes = 200

out_path = choose.dir(default = out_base_path, "Select Folder To Save CSV Data To:")

# Load the 10X dataset
obj.data <- Read10X(data.dir = data_path)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
obj <- CreateSeuratObject(raw.data = obj.data, min.cells = 3, min.genes = min_genes, 
                           project = "10X_PBMC")

# Convert dgTMatrix to regular matrix
data_df = as.data.frame(as.matrix(obj@raw.data))

# Form output file name based on input data directory and output path chosen
out_file = paste(out_path, "\\", out_name, ".csv", sep = "")

#write.csv(data_df, file = out_file, row.names=TRUE, col.names=TRUE)
write.csv(data_df, file = out_file)
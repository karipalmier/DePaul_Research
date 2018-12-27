#####################################################################################################################
#
#  CreateAdjToGraphml.R
#
#
#  Date         Name          Comment
#  9/21/2018    Kari Palmier  File Created
#
#
####################################################################################################################

library("intergraph")
library("igraph")

df_filename = "C:\\DePaulCoursework\\Research\\PBMC_Examples\\Bryan_WebLink_PBMC_4K\\2018921_0_9_53_Bicluster_Results\\2018921_3_46_6_Bicluster_Filtered_97_Adj_Matrix.csv"

graphml_ndx = regexpr(".graphml", df_filename)
base_name = substring(df_filename, 1, (graphml_ndx[1] - 1))
graphml_filename = paste(base_name, ".graphml", sep = "")

bicluster_df = read.csv(df_filename, header = T, row.names = 1, sep = ",", stringsAsFactors = FALSE)

bicluster_matrix = data.matrix(bicluster_df)

bicluster_gr = graph_from_adjacency_matrix(bicluster_matrix,  mode = "undirected", weighted = TRUE)

V(bicluster_gr)$wdegree = graph.strength(bicluster_gr)

# Save updated files
write_graph(bicluster_gr, graphml_filename, format = "graphml")




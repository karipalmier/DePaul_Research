#####################################################################################################################
#
#  CreateGraphmlFile.R
#
#
#  Date         Name          Comment
#  7/15/2018    Kari Palmier  File Created
#
#
####################################################################################################################

library("intergraph")
library("igraph")

# Get command line arguments 
args = commandArgs(trailingOnly=TRUE)

# Get inputs
df_filename = args[1]
graphml_base_filename = args[2]

#df_filename = "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_1_22_Cluster_Results\\2018715_21_44_53_Bicluster_Distance_Matrix.csv"
#graphml_base_filename = "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_1_22_Cluster_Results\\2018715_21_44_53_Bicluster_Graph.graphml"

graphml_ndx = regexpr(".graphml", graphml_base_filename)
base_name = substring(graphml_base_filename, 1, (graphml_ndx[1] - 1))
graphml_filename25 = paste(base_name, "_Filter_Low25.graphml", sep = "")
graphml_filename50 = paste(base_name, "_Filter_Low50.graphml", sep = "")
graphml_filename75 = paste(base_name, "_Filter_Low75.graphml", sep = "")

bicluster_df = read.csv(df_filename, header = T, row.names = 1, sep = ",", stringsAsFactors = FALSE)

bicluster_matrix = data.matrix(bicluster_df)

bicluster_gr = graph_from_adjacency_matrix(bicluster_matrix,  mode = "undirected", weighted = TRUE)

# edge_w = E(bicluster_gr)$weight
# min_e = min(edge_w)
# max_e = max(edge_w)
# diff_e = max_e - min_e
# e_25_thresh = ceiling(diff_e / 4)
# e_50_thresh = ceiling(diff_e / 2)
# e_75_thresh = ceiling((diff_e / 4) * 3)

bicluster_filter = delete.edges(bicluster_gr, which(E(bicluster_gr)$weight == 1))
bicluster_filter = delete.vertices(bicluster_filter, which(degree(bicluster_filter) == 0))
V(bicluster_filter)$wdegree = graph.strength(bicluster_filter)

edge_w = E(bicluster_filter)$weight
sorted_edge_w = sort(edge_w)
e_25_thresh = quantile(sorted_edge_w, 0.25)
e_50_thresh = quantile(sorted_edge_w, 0.50)
e_75_thresh = quantile(sorted_edge_w, 0.75)

bicluster_filter_low25 = delete.edges(bicluster_gr, which(E(bicluster_gr)$weight < e_25_thresh))
bicluster_filter_low25 = delete.vertices(bicluster_filter_low25, which(degree(bicluster_filter_low25) == 0))
V(bicluster_filter_low25)$wdegree = graph.strength(bicluster_filter_low25)

bicluster_filter_low50 = delete.edges(bicluster_gr, which(E(bicluster_gr)$weight < e_50_thresh))
bicluster_filter_low50 = delete.vertices(bicluster_filter_low50, which(degree(bicluster_filter_low50) == 0))
V(bicluster_filter_low50)$wdegree = graph.strength(bicluster_filter_low50)

bicluster_filter_low75 = delete.edges(bicluster_gr, which(E(bicluster_gr)$weight < e_75_thresh))
bicluster_filter_low75 = delete.vertices(bicluster_filter_low75, which(degree(bicluster_filter_low75) == 0))
V(bicluster_filter_low75)$wdegree = graph.strength(bicluster_filter_low75)

# Save updated files
write_graph(bicluster_filter_low25, graphml_filename25, format = "graphml")
write_graph(bicluster_filter_low50, graphml_filename50, format = "graphml")
write_graph(bicluster_filter_low75, graphml_filename75, format = "graphml")




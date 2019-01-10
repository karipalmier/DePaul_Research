#####################################################################################################################
#
#  CreateFilteredAdjacency.R
#
#  This script converts a given distance matrix to an igraph object, filters the graph by the percentage passed in,
#  calculates weighted degree and betweenness centraility and several modularity values, converts the filtered graph
#  to an adjacency matrix and saves it, then returns the centraility and modularity information per node to Python.
#
#  Input Arguments:
#     df_filename = file containing distance matrix data frame creating in Python
#     adj_filename = file to same output filtered adjacency matrix to
#     filter_prop = proportion of edge weights to be filtered out (for example 0.95 filteres out lowest 95% of edges)
#     save_graphml = flag for whether or not to save graphml file
#     save_hist = flag for whether or not to save histograms to jpg files
#
#  Outputs:
#     node_attrs: list of following
#                   names = list of node gene names 
#                   wdegree_cent = list of node weighted degree centrailites 
#                   between_cent = list of node betweenness centrailities 
#                   louvain_mod = list of node louvain modularity communities 
#                   leade_mod = list of node leading eigenvector modularity communities
#                   fastg_mod = list of node fast greedy modularity communities
#                   ebetween_mod = list of node edge betweenness modularity communities
#
#  Date         Name          Comment
#  8/29/2018    Kari Palmier  File Created
#  10/10/2018   Kari Palmier  Added modularity and edge weight histograms
#  1/9/2019     Kari Palmier  Updated to return full graph modularites and centrailities (not components separately)
#
#
#####################################################################################################################

library("intergraph")
library("igraph")
library("ggplot2")

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)

# Get inputs
df_filename = args[1]
adj_filename = args[2]
filter_prop = args[3]
save_graphml = args[4]
save_hist = args[5]

# df_filename = "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_PBMC_4K\\20181010_22_56_36_Bicluster_Results\\20181011_10_59_34_Bicluster_Distance_Matrix.csv"
# adj_filename = "C:\\Temp\\20181011_10_59_34_Bicluster_Filtered_99_Adj_Matrix.csv"
# filter_prop = "99"
# save_graphml = "TRUE"
# save_hist = "TRUE"

# Convert inputs
filter_prop = as.numeric(filter_prop)

if (filter_prop > 1){
  filter_prop = filter_prop / 100
}


#################### Create Filtered Graphml Network Object #########################################################

bicluster_df = read.csv(df_filename, header = T, row.names = 1, sep = ",", stringsAsFactors = FALSE)

bicluster_matrix = data.matrix(bicluster_df)

bicluster_gr = graph_from_adjacency_matrix(bicluster_matrix,  mode = "undirected", weighted = TRUE)

bicluster_filter = delete.edges(bicluster_gr, which(E(bicluster_gr)$weight == 1))
bicluster_filter = delete.vertices(bicluster_filter, which(degree(bicluster_filter) == 0))
V(bicluster_filter)$wdegree = graph.strength(bicluster_filter)

edge_w = E(bicluster_filter)$weight
sorted_edge_w = sort(edge_w)
e_thresh = quantile(sorted_edge_w, filter_prop)

bicluster_filter_prop = delete.edges(bicluster_gr, which(E(bicluster_gr)$weight <= e_thresh))
bicluster_filter_prop = delete.vertices(bicluster_filter_prop, which(degree(bicluster_filter_prop) == 0))
V(bicluster_filter_prop)$wdegree = graph.strength(bicluster_filter_prop)
V(bicluster_filter_prop)$betweenness = betweenness(bicluster_filter_prop, normalized=TRUE)

min_edge = min(E(bicluster_filter_prop)$weight)
max_edge = max(E(bicluster_filter_prop)$weight)


############################ Create Edge Weight Histograms ##########################################################

graphml_ndx = regexpr("_Adj_Matrix.csv", adj_filename)
base_name = substring(adj_filename, 1, (graphml_ndx[1] - 1))

if (save_hist == "TRUE"){
  jpg_all_path = paste(base_name, "_Hist_All_Edges.jpg", sep = "")
  jpeg(file = jpg_all_path)
}

g_all = ggplot(data = data.frame(weights = E(bicluster_gr)$weight), aes(x=weights)) +
  geom_histogram() +
  ggtitle("All Edge Weights Histogram") + 
  labs(x = "Edge Weights", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

if (save_hist == "TRUE"){
  garbage = dev.off()
}


if (save_hist == "TRUE"){
  jpg_filtered_path = paste(base_name, "_Hist_Edges.jpg", sep = "")
  jpeg(file = jpg_filtered_path)
}

g_filtered = ggplot(data = data.frame(weights = E(bicluster_filter_prop)$weight), aes(x=weights)) +
  geom_histogram() +
  ggtitle(paste(filter_prop, "Filtered Edge Weights Histogram")) + 
  labs(x = "Edge Weights", y = "Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

if (save_hist == "TRUE"){
  garbage = dev.off()
}


#################### Create Adjacency Matrix and DataFrame #########################################################

num_nodes = vcount(bicluster_filter_prop)
filtered_adj_filename = paste(base_name, "_", num_nodes, "_Nodes_Adj_Matrix.csv")
bicluster_filter_mat = as.matrix(get.adjacency(bicluster_filter_prop, type = "both", names = TRUE, attr = "weight"))
bicluster_filter_df = as.data.frame(bicluster_filter_mat)

write.csv(bicluster_filter_df, file = filtered_adj_filename)

if (save_graphml == "TRUE"){
  graphml_filename = paste(base_name, "_EdgeWs_", min_edge, "to", max_edge, ".graphml", sep = "")
  write_graph(bicluster_filter_prop, graphml_filename, format = "graphml")
  
}


############################ Find Modularity Communities ############################################################

cluster_filter_names = V(bicluster_filter_prop)$name
cluster_filter_wdeg = V(bicluster_filter_prop)$wdegree
cluster_filter_betw = V(bicluster_filter_prop)$betweenness

cluster_filter_louvain = rep(0, length(cluster_filter_names))
cluster_filter_leade = rep(0, length(cluster_filter_names))
cluster_filter_fastg = rep(0, length(cluster_filter_names))
cluster_filter_ebetween = rep(0, length(cluster_filter_names))

bicluster_decomp = decompose(bicluster_filter_prop)

num_comps = length(bicluster_decomp)
lv_offset = 0
le_offset = 0
fg_offset = 0
bt_offset = 0
for (i in 1:num_comps){
  
  bicluster_temp = bicluster_decomp[[i]]
  temp_nodes = V(bicluster_temp)$name
  temp_num_nodes = vcount(bicluster_temp)
  
  comm_lv = cluster_louvain(bicluster_temp)
  V(bicluster_temp)$comm_lv = comm_lv$membership
  
  comm_le = cluster_leading_eigen(bicluster_temp)
  V(bicluster_temp)$comm_le = comm_le$membership
  
  comm_fg = cluster_fast_greedy(bicluster_temp)
  V(bicluster_temp)$comm_fg = comm_fg$membership
  
  comm_bt = cluster_edge_betweenness(bicluster_temp, weights = NULL)
  V(bicluster_temp)$comm_bt = comm_bt$membership
  
  for (n in 1:temp_num_nodes){
    node_ndx = cluster_filter_names == temp_nodes[n]
    cluster_filter_louvain[node_ndx] = V(bicluster_temp)$comm_lv[n] + lv_offset
    cluster_filter_leade[node_ndx] = V(bicluster_temp)$comm_le[n] + le_offset
    cluster_filter_fastg[node_ndx] = V(bicluster_temp)$comm_fg[n] + fg_offset
    cluster_filter_ebetween[node_ndx] = V(bicluster_temp)$comm_bt[n] + bt_offset
    
  }
  
  lv_offset = max(V(bicluster_temp)$comm_lv)
  le_offset = max(V(bicluster_temp)$comm_le)
  fg_offset = max(V(bicluster_temp)$comm_fg)
  bt_offset = max(V(bicluster_temp)$comm_bt)
  
}

# Create list of node attributes
node_attrs = list(names = cluster_filter_names,
                  wdegree_cent = cluster_filter_wdeg,
                  between_cent = cluster_filter_betw,
                  louvain_mod = cluster_filter_louvain,
                  leade_mod = cluster_filter_leade,
                  fastg_mod = cluster_filter_fastg,
                  ebetween_mod = cluster_filter_ebetween)

# Convert combined output list to JSON to pass to Python
jsonData = jsonlite::toJSON(node_attrs, pretty = TRUE)

# Print combined output JSON to Python
print(jsonData)








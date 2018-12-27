#####################################################################################################################
#
#  CreateFilteredAdjacency.R
#
#
#  Input Arguments:
#     df_filename = file containing distance matrix data frame creating in Python
#     adj_filename = file to same output filtered adjacency matrix to
#     filter_prop = proportion of edge weights to be filtered out (for example 0.95 filteres out lowest 95% of edges)
#
#  Outputs:
#     None
#
#  Date         Name          Comment
#  8/29/2018    Kari Palmier  File Created
#  10/10/2018   Kari Palmier  Added modularity and edge weight histograms
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
save_components = args[6]

# df_filename = "C:\\DePaulCoursework\\Research\\PBMC_Examples\\Bryan_WebLink_PBMC_4K\\20181010_22_56_36_Bicluster_Results\\20181011_10_59_34_Bicluster_Distance_Matrix.csv"
# adj_filename = "C:\\DePaulCoursework\\Research\\PBMC_Examples\\Bryan_WebLink_PBMC_4K\\20181010_22_56_36_Bicluster_Results\\20181011_10_59_34_Bicluster_Filtered_99_Adj_Matrix.csv"
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
print(g_all)

if (save_hist == "TRUE"){
  dev.off()
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
print(g_filtered)

if (save_hist == "TRUE"){
  dev.off()
}


#################### Create Adjacency Matrix and DataFrame #########################################################

num_nodes = vcount(bicluster_filter_prop)
filtered_adj_filename = paste(base_name, "_AllComp_", num_nodes, "_Nodes_Adj_Matrix.csv")
bicluster_filter_mat = as.matrix(get.adjacency(bicluster_filter_prop, type = "both", names = TRUE, attr = "weight"))
bicluster_filter_df = as.data.frame(bicluster_filter_mat)

write.csv(bicluster_filter_df, file = filtered_adj_filename)

if (save_graphml == "TRUE"){
  graphml_filename = paste(base_name, "_AllComp_EdgeWs_", min_edge, "to", max_edge, ".graphml", sep = "")
  
  write_graph(bicluster_filter_prop, graphml_filename, format = "graphml")
  
}


############################ Find Modularity Communities ############################################################

if (save_hist == "TRUE"){
  bicluster_decomp = decompose(bicluster_filter_prop)

  num_comps = length(bicluster_decomp)
  for (i in 1:num_comps){
    
    bicluster_temp = bicluster_decomp[[i]]
    temp_num_nodes = vcount(bicluster_temp)
    
    min_edge_temp = min(E(bicluster_temp)$weight)
    max_edge_temp = max(E(bicluster_temp)$weight)
    
    comm_sg = cluster_spinglass(bicluster_temp)
    V(bicluster_temp)$comm_sg = comm_sg$membership
    
    comm_lv = cluster_louvain(bicluster_temp)
    V(bicluster_temp)$comm_lv = comm_lv$membership
    
    comm_le = cluster_leading_eigen(bicluster_temp)
    V(bicluster_temp)$comm_le = comm_le$membership
    
    comm_fg = cluster_fast_greedy(bicluster_temp)
    V(bicluster_temp)$comm_fg = comm_fg$membership
    
    comm_bt = cluster_edge_betweenness(bicluster_temp, weights = NULL)
    V(bicluster_temp)$comm_bt = comm_bt$membership
    
    comm_wt3 = cluster_walktrap(bicluster_temp, steps = 3)
    V(bicluster_temp)$comm_wt5 = comm_wt3$membership
    
    comm_wt5 = cluster_walktrap(bicluster_temp, steps = 5)
    V(bicluster_temp)$comm_wt10 = comm_wt5$membership
    
    comm_wt7 = cluster_walktrap(bicluster_temp, steps = 7)
    V(bicluster_temp)$comm_wt10 = comm_wt7$membership
  
  
    if (save_hist == "TRUE"){
      jpg_temp_path = paste(base_name, "_Hist_Edges_Comp", i, ".jpg", sep = "")
      jpeg(file = jpg_temp_path)
    }
    
    g_temp = ggplot(data = data.frame(weights = E(bicluster_temp)$weight), aes(x=weights)) +
      geom_histogram() +
      ggtitle(paste(filter_prop, "Component", i, "Filtered Edge Weights Histogram")) + 
      labs(x = "Edge Weights", y = "Count") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(g_temp)
    
    if (save_hist == "TRUE"){
      dev.off()
    }
    
    temp_adj_filename = paste(base_name, "_Adj_Matrix_Comp", i, ".csv", sep = "")
    bicluster_mat_temp = as.matrix(get.adjacency(bicluster_temp, type = "both", names = TRUE, attr = "weight"))
    bicluster_df_temp = as.data.frame(bicluster_mat_temp)
    
    write.csv(bicluster_df_temp, file = temp_adj_filename)
    
    if (save_graphml == "TRUE"){
      graphml_filename = paste(base_name, "_Comp", i, "_", temp_num_nodes, "_Nodes_EdgeWs_", min_edge_temp, "to", max_edge_temp, ".graphml", sep = "")
      
      write_graph(bicluster_temp, graphml_filename, format = "graphml")
      
    }
    
  }
}






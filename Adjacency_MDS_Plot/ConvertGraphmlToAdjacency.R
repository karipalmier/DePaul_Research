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

graph25_filename = "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_44_53_Bicluster_Graph_Filter_Low25.graphml"
graph50_filename = "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_44_53_Bicluster_Graph_Filter_Low50.graphml"
graph75_filename = "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_44_53_Bicluster_Graph_Filter_Low75.graphml"

bicluster_graph25 = read.graph(graph25_filename, format = "graphml")
bicluster_graph50 = read.graph(graph50_filename, format = "graphml")
bicluster_graph75 = read.graph(graph75_filename, format = "graphml")

bicluster_mat25 = as.matrix(get.adjacency(bicluster_graph25, type = "both", names = TRUE, attr = "weight"))
bicluster_df25 = as.data.frame(bicluster_mat25)

bicluster_mat50 = as.matrix(get.adjacency(bicluster_graph50, type = "both", names = TRUE, attr = "weight"))
bicluster_df50 = as.data.frame(bicluster_mat50)

bicluster_mat75 = as.matrix(get.adjacency(bicluster_graph75, type = "both", names = TRUE, attr = "weight"))
bicluster_df75 = as.data.frame(bicluster_mat75)

adj25_filename = "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_44_53_Bicluster_Adjacency_Low25.csv"
adj50_filename = "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_44_53_Bicluster_Adjacency_Low50.csv"
adj75_filename = "C:\\DePaulCoursework\\Research\\BryanExample\\2018715_21_44_53_Bicluster_Adjacency_Low75.csv"

write.csv(bicluster_df25, file = adj25_filename)
write.csv(bicluster_df50, file = adj50_filename)
write.csv(bicluster_df75, file = adj75_filename)



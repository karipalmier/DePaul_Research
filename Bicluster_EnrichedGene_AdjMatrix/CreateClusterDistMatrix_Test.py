# -*- coding: utf-8 -*-
"""
CreateClusterDistMatrix_Test.py

This function calls the create_cluster_wrapper function from the
CreateClusterDistMatrix.py file, which creates the distance matrix for multiple
clusterings, then calls the create_adj_matrixfunction from the same py file to
create the filtered adjacency matrix.  
This is done to test out the use cases of the Python and R scripts in the 
Django framework.

Date      Name              Comment
8/29/18   Kari Palmier      Created
1/9/19    Kari Palmier      Updated to handle PC component and graph component
                            saving removals
                            
"""

from CreateClusterDistMatrix import create_cluster_wrapper
from CreateClusterDistMatrix import create_adj_matrix

# Define command and arguments
command = 'Rscript'

#dir_path = os.path.dirname(os.path.realpath(__file__))
R_dir_path = "C:\\DePaulCoursework\\Research\\Code\\"

#cluster_data_path = \
#  "C:\\DePaulCoursework\\Research\\Known_Examples\\Satija_Example_PBMC_4K\\hg19.csv"
#cluster_data_path = \
#  "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_PBMC_4K\\GRCh38.csv"
#cluster_data_path = \
#  "C:\\DePaulCoursework\\Research\\BryanExample\\Bryan_Example_Sample1_infected.csv"
#cluster_data_path = \
#  "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_PBMC_8K\\GRCh38.csv"
cluster_data_path = \
  "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_Mouse_Neurons\\mm10.csv"

# Enriched significanct level
sig_level = '0.05'
#sig_level = '0.00000000001'

test_types = ['roc']
#test_types = ['roc', 'bimodal']

#res_levels = [0.6, 1.2]
res_levels = [0.6, 0.8, 1.0, 1.2, 1.4]

#ngene_up_limit = [5000]
#ngene_low_limit = [1]
#perc_mito_up_limit = [0.04]
#perc_mito_low_limit = ["-Inf"]
ngene_up_limit = [2500]
ngene_low_limit = [200]
#ngene_up_limit = [2500, 5000]
#ngene_low_limit = [200, 1]
perc_mito_up_limit = [0.05, 0.04]
#perc_mito_up_limit = [0.05, 0.04, 0.6]
perc_mito_low_limit = ["-Inf"]

cluster_dims_used = [15]
#cluster_dims_used = [10, 15]
markers_min_pct = 0.25
markers_thresh_use = 0.25

#filter_colsums_5000 = 'TRUE'
filter_colsums_5000 = 'FALSE'
print_log = True
save_cluster_data = True
log_to_file = True

output_rule_path, output_dist_path = create_cluster_wrapper(R_dir_path, 
                                        cluster_data_path, test_types, res_levels, 
                                        sig_level, ngene_up_limit, ngene_low_limit, 
                                        perc_mito_up_limit, perc_mito_low_limit,  
                                        cluster_dims_used, markers_min_pct, 
                                        markers_thresh_use, 
                                        filter_colsums_5000, print_log, 
                                        save_cluster_data, log_to_file)

save_graphml = "TRUE"
save_hist = "TRUE"

output_dist_path = "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_Mouse_Neurons\\201919_18_9_4_Bicluster_Results\\201919_19_0_55_Bicluster_Distance_Matrix.csv"

filter_perc = 97
adj_filename, node_filename = create_adj_matrix(R_dir_path, filter_perc, output_dist_path, save_graphml, 
                                                save_hist)

filter_perc = 99
adj_filename, node_filename = create_adj_matrix(R_dir_path, filter_perc, output_dist_path, save_graphml, 
                                                save_hist)

filter_perc = 99.5
adj_filename, node_filename = create_adj_matrix(R_dir_path, filter_perc, output_dist_path, save_graphml, 
                                                save_hist)


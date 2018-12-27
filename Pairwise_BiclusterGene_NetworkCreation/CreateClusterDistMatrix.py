# -*- coding: utf-8 -*-
"""
CreateClusterDistMatrix.py

The functions in this file are used to perform Seurat biclustering given test types and  
resolutions, perform hypergeometric enrichment testing to get genes enriched per cluster 
in each clustering, create a distance matrix of number of times a gene pair is enriched 
in the same cluster, convert the distance matrix to a graph object and filter low value edges, 
then create an adjacency matrix from the filtered graph object.

Date      Name              Comment
7/11/18   Kari Palmier      Created
8/29/18   Kari Palmier      Commented Code
10/4/18   Kari Palmier      Added in additional paramters for clustering grid search

"""

import subprocess
import json
import pandas as pd
import numpy as np
import time
import datetime
import sys
import os

# Function to generate the current date and time string used in creating file names
def get_current_time():
    curr_datetime = datetime.datetime.now()
    curr_date = str(curr_datetime.year) + str(curr_datetime.month) + str(curr_datetime.day)
    curr_time = str(curr_datetime.hour) + "_" + str(curr_datetime.minute) + "_" + \
    str(curr_datetime.second)
    curr_datetime = curr_date + "_" + curr_time
    return curr_datetime

# Function to perform clustering, find enriched genes, and create distance matrix
# This function calls the PerformBicluster and GetEnrichedGO.R scripts.  
def create_cluster_dist_matrix(R_dir_path = "C:\\", sig_level = 0.05, 
                               cluster_data_path = "C:\\", cluster_out_path = "C:\\", 
                               test_types = ["roc"], res_levels = [0.05], 
                               ngene_up_limit = [5000], ngene_low_limit = [1], 
                               perc_mito_up_limit = [0.05], perc_mito_low_limit = ["-Inf"], 
                               pca_num_pcs = [20], cluster_dims_used = [15], 
                               markers_min_pct = 0.25, markers_thresh_use = 0.25, 
                               filter_colsums_5000 = "TRUE", print_log = False, 
                               save_cluster_data = False):
    
    if print_log:
        start_script_time = time.time()
        curr_datetime = get_current_time()
        old_stdout = sys.stdout
        log_filename = cluster_out_path + curr_datetime + "_Bicluster_Log.txt"
        log_file = open(log_filename, "w")
        sys.stdout = log_file
        print("Script Start Time:", curr_datetime)
               
    cluster_script_path = R_dir_path + 'PerformBicluster.R'
    enriched_script_path = R_dir_path + 'GetEnrichedGO.R'

    cluster_jpg_path = cluster_out_path + "Heatmpas\\"
    os.mkdir(cluster_jpg_path)
    cluster_inter_path = cluster_out_path + "Intermediate\\"
    os.mkdir(cluster_inter_path)   
                             
    rule_dict = {}
    
    num_clusterings = len(test_types) * len(res_levels) * len(ngene_up_limit) * len(ngene_low_limit) * \
    len(perc_mito_up_limit) * len(perc_mito_low_limit) * len(pca_num_pcs) * len(cluster_dims_used)
    
    if print_log:
        print("\n")
        print("Clustering Parameters:")
        print("Significance Level:", sig_level)
        print("Test Type:", test_types)
        print("Resolution Levels:", res_levels)
        print("Number of Genes Upper Limit:", ngene_up_limit)
        print("Number of Genes Lower Limit:", ngene_low_limit)
        print("Percent Mitochondrial Upper Limit:", perc_mito_up_limit)
        print("Percent Mitochondrial Lower Limit:", perc_mito_low_limit)
        print("Number of PCA Primary Components:", pca_num_pcs)
        print("Number of CLuster Dimensions:", cluster_dims_used)
        print("Marker Minimum Percent:", markers_min_pct)
        print("Marker Threshold:", markers_thresh_use)
        print("Filter Col Sums Below 5000:", filter_colsums_5000)
        print("Total number of clusterings:", num_clusterings)
        print("\n")
        
    
    # Loop over possible parameters
    num_iter = 0
    clustering_count = 0
    for num_pca_pcs in pca_num_pcs:
        for num_cluster_dims in cluster_dims_used:
            for ngene_ul in ngene_up_limit:
                for ngene_ll in ngene_low_limit:
                    for perc_mito_ul in perc_mito_up_limit:
                         for perc_mito_ll in perc_mito_low_limit:
                            for test_type in test_types:                                    
                                for res_level in res_levels:
                                    if print_log:
                                        print("==========================================================")
                                        print("Current Number of PCA PCs:", num_pca_pcs)
                                        print("Current Number of Cluster Dimensions:", num_cluster_dims)
                                        print("Current NGene Upper Limit:", ngene_ul)
                                        print("Current NGene Lower Limit:", ngene_ll)
                                        print("Current Perc Mito Upper Limit:", perc_mito_ul)
                                        print("Current Perc Mito Lower Limit:", perc_mito_ll)
                                        print("Current Test Type:", test_type)
                                        print("Current Resolution:", res_level)
                                        start_loop_time = time.time()
                                                                     
                                    cluster_arg_list = [cluster_data_path, 'TRUE', 'TRUE', str(res_level),  
                                                        str(num_pca_pcs), str(num_cluster_dims),  
                                                        str(markers_min_pct), test_type,
                                                        str(markers_thresh_use), str(ngene_ul), 
                                                        str(ngene_ll),  str(perc_mito_ul), 
                                                        str(perc_mito_ll), cluster_jpg_path, 
                                                        filter_colsums_5000]
                                    
                                    # Build subprocess command
                                    cluster_cmd = ["Rscript", cluster_script_path] + cluster_arg_list
                                    
                                    # Call R PerformBicluster function
                                    # check_output will run the command and store to result
                                    cluster_json = subprocess.check_output(cluster_cmd, 
                                                                           universal_newlines = True)
                                    cluster_dict = json.loads(cluster_json)
                                    
                                    # Set argument to R script - note this will need to be passed from 
                                    # another layer in the final code
                                    gene_cluster_df = pd.DataFrame(cluster_dict['gene_cluster_data'])
                                    gene_cluster_df['cluster'] = gene_cluster_df['cluster'].astype('int')
                                    cluster_nums = list(gene_cluster_df['cluster'].unique())
                                    
                                    # Save gene cluster assignments for current test type and  
                                    # resolution values
                                    if save_cluster_data:                                
                                        curr_datetime_str = get_current_time()
                                        output_cluster_path = cluster_inter_path + curr_datetime_str + \
                                        "_res_" + str(res_level) + "_test_" + test_type + \
                                        "_ngene_" + str(ngene_ul) + "_" + str(ngene_ll) + "_mito_" + \
                                        str(perc_mito_ul) + "_" + str(perc_mito_ll) + "_pcas_" + \
                                        str(num_pca_pcs) + "_cdims_" + str(num_cluster_dims) + ".csv"
                                        gene_cluster_df.to_csv(output_cluster_path, header = True, 
                                                               index = True)
                                    
                                    if print_log:
                                        end_cluster_time = time.time()
                                        cluster_duration_s = end_cluster_time - start_loop_time
                                        cluster_duration_min = cluster_duration_s / 60
                                        print("\n")
                                        print("Cluster Duration (sec):", cluster_duration_s)
                                        print("Cluster Duration (min):", cluster_duration_min)
                                        
                                    # Find enriched genes per cluster and create/update pairwise gene count 
                                    # data frame
                                    clustering_count += 1
                                    for cluster in cluster_nums:
                                        if print_log:
                                            start_dist_time = time.time()
                                        
                                        if print_log:
                                            print("\n")
                                            print("Cluster:", cluster)
                                            
                                        # Get data and list of genes belonging to current cluster
                                        cluster_ndx = gene_cluster_df['cluster'] == cluster
                                        temp_df = gene_cluster_df[cluster_ndx]
                                        cluster_gene_list = list(temp_df['gene'])
                                    
                                        # Build subprocess command
                                        enriched_cmd = ["Rscript", enriched_script_path] + [sig_level] + \
                                        cluster_gene_list
                                        
                                        # Call GetEnrichedGO.R function
                                        # check_output will run the command and store to result
                                        enriched_json = subprocess.check_output(enriched_cmd, 
                                                                                universal_newlines=True)
                                        enriched_dict = json.loads(enriched_json)
                                        
                                        # Get list of enriched GO genes for current cluster
                                        num_enriched_genes = enriched_dict['num_enriched_genes'][0]
                                        curr_genes = list(enriched_dict['enriched_genes'])
                                    
                                         # Build dictionary with list of genes per cluster as each value.
                                        # For use in association rule mining.
                                        dict_str = str(clustering_count) + '_' + str(cluster + 1)
                                        rule_dict[dict_str] = curr_genes
                                        
                                       # If the first iteration, create arrays and dataframes to be 
                                        # populated
                                        if num_iter == 0:   
                                            temp_array = np.zeros((num_enriched_genes, num_enriched_genes), 
                                                                  dtype = int)
                                            dist_df = pd.DataFrame(temp_array)
                                            dist_df.columns = curr_genes
                                            dist_df.index = curr_genes
                                            
                                        # Else find the list of genes not present in the current dataframe, 
                                        # then add columns and rows to the dataframe for each
                                        else:
                                            dist_genes = list(dist_df.columns.values)
                                            
                                            unmatched_genes = list(set(curr_genes) - set(dist_genes))
                                            
                                            if len(unmatched_genes) > 0:
                                                for gene in unmatched_genes:
                                                    # Get current dataframe row and column sizes
                                                    dist_rows = dist_df.shape[0]
                                                    dist_cols = dist_df.shape[1]
                                                    dist_index = list(dist_df.index)
                                                    
                                                    # Create temp column and row of all zeros 
                                                    temp_col = np.zeros((dist_rows, 1), dtype = int)
                                                    num_new_cols = dist_cols + 1
                                                    temp_row = np.zeros((1, num_new_cols), dtype = int)
                                                    
                                                    # Add temp column and row to existing dataframe 
                                                    dist_df[gene] = pd.DataFrame(temp_col, 
                                                           index = dist_index)
                                                    new_cols = list(dist_df.columns)
                                                    new_row_df = pd.DataFrame(temp_row, columns = new_cols)
                                                    new_row_df.index = [gene]
                                                    dist_df = pd.concat([dist_df, new_row_df])
                                                    
                                        # Increment each cell corresponding to a gene pair by 1 (unless  
                                        # same gene is row and column)
                                        for row_gene in curr_genes:
                                            for col_gene in curr_genes:
                                                if row_gene == col_gene:
                                                    continue
                                                
                                                dist_df.loc[row_gene, col_gene] = \
                                                dist_df.loc[row_gene, col_gene] + 1
                                        
                                        num_iter += 1
                                        
                                        if print_log:
                                            end_dist_time = time.time()
                                            dist_duration_s = end_dist_time - start_dist_time
                                            dist_duration_min = dist_duration_s / 60
                                            print("Cluster Matrix Duration (sec):", dist_duration_s)
                                            print("Cluster Matrix Duration (min):", dist_duration_min)
                                            print("\n")
                                                             
                                    if print_log:
                                        end_loop_time = time.time()
                                        loop_duration_s = end_loop_time - start_loop_time
                                        loop_duration_min = loop_duration_s / 60
                                        print("\n")
                                        print("Loop Duration (sec):", loop_duration_s)
                                        print("Loop Duration (min):", loop_duration_min)
                                        print("\n")
        
    # Get current time stamp to be used in distance matrix file naming
    curr_datetime = get_current_time()
    
    # Save distance matrix
    output_dist_path = cluster_out_path + curr_datetime + "_Bicluster_Distance_Matrix.csv"
    dist_df.to_csv(output_dist_path, header = True, index = True)
       
    # Save association rule dictionary
    output_rule_path = cluster_out_path + curr_datetime + "_Bicluster_Rules_Dict.txt"
    dict_file = open(output_rule_path, 'w')
    rule_json = json.dumps(rule_dict)
    dict_file.write(rule_json)
    dict_file.close()
       
    if print_log:
        print("\n")
        print("Distance Matrix File Path:")
        print(output_dist_path)
        print("\n")
        print("Rule Dictionary File Path:")
        print(output_rule_path)
        
        min_edge_w = min(dist_df.min())
        max_edge_w = max(dist_df.max())
        print("\n")
        print("Min Overall Pair Count:", min_edge_w)
        print("Max Overall Pair Count:", max_edge_w)
    
        end_script_time = time.time()
        script_duration_s = end_script_time - start_script_time
        script_duration_min = script_duration_s / 60  
        print("\n")
        print("Script Duration (sec):", script_duration_s)
        print("Script Duration (min):", script_duration_min)

        # Assign standard output back to terminal
        curr_datetime = get_current_time()
        print("\n")
        print("Script Stop Time:", curr_datetime)
        sys.stdout = old_stdout
        log_file.close()

    # Remove temporary distance matrix file
    #os.remove(output_dist_path)
    
    return output_rule_path, output_dist_path


# Call function to create filtered adjacency matrix
def create_adj_matrix(R_dir_path = "C:\\", filter_perc = 97, output_dist_path = "C:\\", 
                      save_graphml = "FALSE", save_hist = "FALSE", save_components = "FALSE"):

    adj_script_path = R_dir_path + 'CreateFilteredAdjacency.R'

    # If filter is a proportion, convert to a percentage
    if filter_perc < 1:
        filter_perc = filter_perc * 100
        
    # Create filter value string to be used to create output adjacency file name
    filter_perc_str = str(filter_perc).replace(".", "p")
    
    # Create adjacency file name
    end_ndx = output_dist_path.find("Distance_Matrix")
    base_name = output_dist_path[0:end_ndx]
    filter_prop = filter_perc / 100
    adj_filename = base_name + "Filtered_" + filter_perc_str + "_Adj_Matrix.csv"  
    
    # Build subprocess command
    adj_cmd = ["Rscript", adj_script_path, output_dist_path, adj_filename, str(filter_prop), 
               save_graphml, save_hist, save_components]
    
    # Call CreateFilteredAdjacency.R function
    subprocess.check_output(adj_cmd, universal_newlines = True)            

    return adj_filename
    
# Call wrapper function which calls function to create distance matrix and creates text log
# This function is only used when interested in looking at log data (debugging)
def create_cluster_wrapper(R_dir_path = "C\\", cluster_data_path = "C\\", test_types = ["roc"],
                           res_levels = [0.05], sig_level = 0.05, 
                           ngene_up_limit = [5000], ngene_low_limit = [1], 
                           perc_mito_up_limit = [0.05], perc_mito_low_limit = ["-Inf"], 
                           pca_num_pcs = 20, cluster_dims_used = 15, markers_min_pct = 0.25, 
                           markers_thresh_use = 0.25, filter_colsums_5000 = "TRUE", 
                           print_log = False, save_cluster_data = False):
    
    # Create folder to use to save distance matrix and log files
    curr_datetime = get_current_time()
    slash_ndx = cluster_data_path.rfind('\\')
    base_out_path = cluster_data_path[0:slash_ndx + 1]
    cluster_out_path = base_out_path + curr_datetime + "_Bicluster_Results\\"
    os.mkdir(cluster_out_path)

    # Call function to create distance matrix
    output_rule_path, output_dist_path = create_cluster_dist_matrix(R_dir_path, sig_level, 
                                         cluster_data_path, cluster_out_path, test_types, res_levels, 
                                         ngene_up_limit, ngene_low_limit, perc_mito_up_limit, 
                                         perc_mito_low_limit, pca_num_pcs, cluster_dims_used, 
                                         markers_min_pct, markers_thresh_use, filter_colsums_5000,
                                         print_log, save_cluster_data)

    # Return output rule and distance to Django framework to be used in visualization 
    # generation
    return output_rule_path, output_dist_path

        

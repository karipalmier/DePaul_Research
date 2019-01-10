# -*- coding: utf-8 -*-
"""
CompileResults.py

The functions in this file are used to compile statistics and elapsed time information for 
a previously run biclustering set.  Note that the intermediate csv data needs to be in its own
folder because all of the files in the folder get looped over (expects all of the same csv formats). 

Date      Name              Comment
10/22/18  Kari Palmier      Created

"""

import os
import numpy as np
import pandas as pd

def get_cluster_statistics(data_folder):
    
    dir_list = os.listdir(data_folder)
    
    column_names = ["resolution", "test_type", "num_genes_max", "num_genes_min", "perc_mito_max",
                    "perc_mito_min", "num_pca_pcs", "num_cluster_dims", "tot_num_genes", 
                    "num_uniq_genes", "num_notuniq_genes", "num_clusters", "min_cluster_size", 
                    "max_cluster_size", "mean_cluster_size", "median_cluster_size", 
                    "q1_cluster_size", "q3_cluster_size", "cluster_size_std", "cluster_size_var"]
    stat_df = pd.DataFrame(columns = column_names)
    for x in range(len(dir_list)):
        
        file = dir_list[x]
        name_lst = file.strip().split("_")
        res = float(name_lst[name_lst.index("res") + 1])
        test = name_lst[name_lst.index("test") + 1]
        ngene_max = float(name_lst[name_lst.index("ngene") + 1])
        ngene_min = float(name_lst[name_lst.index("ngene") + 2])
        mito_max = float(name_lst[name_lst.index("percmito") + 1])
        mito_min = float(name_lst[name_lst.index("percmito") + 2])
#        pcas = float(name_lst[name_lst.index("pcas") + 1])
        pcas = 20
        cdims = float(name_lst[name_lst.index("cdims") + 1].split('.')[0])
               
        file_name = data_folder + "\\" + file
        temp_df = pd.read_csv(file_name, delimiter = ',', header = 0, index_col = 0)
        num_tot_genes = len(temp_df["gene"])
        num_uniq_genes = len(pd.unique(temp_df["gene"]))
        num_notuniq_genes = num_tot_genes - num_uniq_genes
        cluster_nums = pd.unique(temp_df["cluster"])
        num_clusters = len(cluster_nums)
        
        cluster_sizes = []
        cluster_sizes = [len(temp_df[temp_df["cluster"] == i]) for i in cluster_nums]
        cluster_sizes = np.array(cluster_sizes)
        
        min_cluster_size = np.min(cluster_sizes)
        max_cluster_size = np.max(cluster_sizes)
        mean_cluster_size = np.mean(cluster_sizes)
        median_cluster_size = np.median(cluster_sizes)
        sd_cluster_size = np.std(cluster_sizes)
        var_cluster_size = np.var(cluster_sizes)
        q3_cluster_size = np.percentile(cluster_sizes, 0.75)
        q1_cluster_size = np.percentile(cluster_sizes, 0.25)
        
        temp_row = [res, test, ngene_max, ngene_min, mito_max, mito_min, pcas, cdims, num_tot_genes, 
                    num_uniq_genes, num_notuniq_genes, num_clusters, min_cluster_size, 
                    max_cluster_size, mean_cluster_size, median_cluster_size, q1_cluster_size, 
                    q3_cluster_size, sd_cluster_size, var_cluster_size]
        temp_row = np.reshape(np.array(temp_row), (1, len(temp_row)))           
        new_row_df = pd.DataFrame(temp_row, columns = column_names)
        new_row_df.index = [x]
        stat_df = pd.concat([stat_df, new_row_df])
        
    return stat_df

def get_elapsed_times(log_file):
    file = open(log_file, "r")

    column_names = ["num_pca_pcs", "num_cluster_dims", "num_genes_max", "num_genes_min", "perc_mito_max",
                    "perc_mito_min", "test_type", "resolution", "elapsed_time_sec", "elapsed_time_min"]
    time_df = pd.DataFrame(columns = column_names)
    entry_ndx = 0
    while True:
        line = file.readline()
        if not line:
            break
        
        if line.find("Current Number of PCA PCs") >= 0:
            pcas = float(line.strip().split(':')[1])
        if line.find("Current Number of Cluster Dimensions") >= 0:
            cdims = float(line.strip().split(':')[1])
        if line.find("Current NGene Upper Limit") >= 0:
            ngene_max = float(line.strip().split(':')[1])
        if line.find("Current NGene Lower Limit") >= 0:
            ngene_min = float(line.strip().split(':')[1])
        if line.find("Current Perc Mito Upper Limit") >= 0:
            mito_max = float(line.strip().split(':')[1])
        if line.find("Current Perc Mito Lower Limit") >= 0:
            mito_min = float(line.strip().split(':')[1])
        if line.find("Current Test Type") >= 0:
            test = line.strip().split(':')[1].strip()
        if line.find("Current Resolution") >= 0:
            res = float(line.strip().split(':')[1])
        if line.find("Cluster Duration (sec)") >= 0:
            elapsed_time_s = float(line.strip().split(':')[1])
        if line.find("Cluster Duration (min)") >= 0:
            elapsed_time_m = float(line.strip().split(':')[1])
            temp_row = [pcas, cdims, ngene_max, ngene_min, mito_max, mito_min, test, res, 
                        elapsed_time_s, elapsed_time_m]
            temp_row = np.reshape(np.array(temp_row), (1, len(temp_row)))           
            new_row_df = pd.DataFrame(temp_row, columns = column_names)
            new_row_df.index = [entry_ndx]
            time_df = pd.concat([time_df, new_row_df])
            entry_ndx += 1
    
    return time_df

def combine_results_times(stat_df, time_df):
    
    elapsed_time_sec_lst = []
    elapsed_time_min_lst = []
    for i in range(len(stat_df)):
        
         match_ndx = (time_df["resolution"] == stat_df.iloc[i]["resolution"]) & \
             (time_df["test_type"] == stat_df.iloc[i]["test_type"]) & \
             (time_df["num_genes_max"] == stat_df.iloc[i]["num_genes_max"]) & \
             (time_df["num_genes_min"] == stat_df.iloc[i]["num_genes_min"]) & \
             (time_df["perc_mito_max"] == stat_df.iloc[i]["perc_mito_max"]) & \
             (time_df["perc_mito_min"] == stat_df.iloc[i]["perc_mito_min"]) & \
             (time_df["num_pca_pcs"] == stat_df.iloc[i]["num_pca_pcs"]) & \
             (time_df["num_cluster_dims"] == stat_df.iloc[i]["num_cluster_dims"])
   
         elapsed_time_sec_lst.append(float(time_df.loc[match_ndx]["elapsed_time_sec"]))
         elapsed_time_min_lst.append(float(time_df.loc[match_ndx]["elapsed_time_min"]))
    
    stat_df["elapsed_time_sec"] = elapsed_time_sec_lst
    stat_df["elapsed_time_min"] = elapsed_time_min_lst
    
    return stat_df
    
def get_summary_results(results_df):
    
    results_df["tot_num_genes"] = results_df.tot_num_genes.astype(int)
    min_tot_num_genes = np.min(results_df["tot_num_genes"])
    max_tot_num_genes = np.max(results_df["tot_num_genes"])
    mean_tot_num_genes = np.mean(results_df["tot_num_genes"])
    median_tot_num_genes = np.median(results_df["tot_num_genes"])
    q3_tot_num_genes = np.percentile(results_df["tot_num_genes"], 0.75)
    q1_tot_num_genes = np.percentile(results_df["tot_num_genes"], 0.25)
    sd_tot_num_genes = np.std(results_df["tot_num_genes"])
    var_tot_num_genes = np.var(results_df["tot_num_genes"])
   
    results_df["num_uniq_genes"] = results_df.num_uniq_genes.astype(int)
    min_uniq_num_genes = np.min(results_df["num_uniq_genes"])
    max_uniq_num_genes = np.max(results_df["num_uniq_genes"])
    mean_uniq_num_genes = np.mean(results_df["num_uniq_genes"])
    median_uniq_num_genes = np.median(results_df["num_uniq_genes"])
    q3_uniq_num_genes = np.percentile(results_df["num_uniq_genes"], 0.75)
    q1_uniq_num_genes = np.percentile(results_df["num_uniq_genes"], 0.25)
    sd_uniq_num_genes = np.std(results_df["num_uniq_genes"])
    var_uniq_num_genes = np.var(results_df["num_uniq_genes"])
    
    results_df["num_clusters"] = results_df.num_clusters.astype(int)
    min_num_clusters = np.min(results_df["num_clusters"])
    max_num_clusters = np.max(results_df["num_clusters"])
    mean_num_clusters = np.mean(results_df["num_clusters"])
    median_num_clusters = np.median(results_df["num_clusters"])
    q3_num_clusters = np.percentile(results_df["num_clusters"], 0.75)
    q1_num_clusters = np.percentile(results_df["num_clusters"], 0.25)
    sd_num_clusters = np.std(results_df["num_clusters"])
    var_num_clusters = np.var(results_df["num_clusters"])
   
    min_elapsed_time_min = np.min(results_df["elapsed_time_min"])
    max_elapsed_time_min = np.max(results_df["elapsed_time_min"])
    mean_elapsed_time_min = np.mean(results_df["elapsed_time_min"])
    median_elapsed_time_min = np.median(results_df["elapsed_time_min"])
    q3_elapsed_time_min = np.percentile(results_df["elapsed_time_min"], 0.75)
    q1_elapsed_time_min = np.percentile(results_df["elapsed_time_min"], 0.25)
    sd_elapsed_time_min = np.std(results_df["elapsed_time_min"])
    var_elapsed_time_min = np.var(results_df["elapsed_time_min"])
   
    temp_row = [min_tot_num_genes, max_tot_num_genes, mean_tot_num_genes, median_tot_num_genes, 
                q3_tot_num_genes, q1_tot_num_genes, sd_tot_num_genes, var_tot_num_genes, 
                min_uniq_num_genes, max_uniq_num_genes, mean_uniq_num_genes, median_uniq_num_genes,
                q3_uniq_num_genes, q1_uniq_num_genes, sd_uniq_num_genes, var_uniq_num_genes, 
                min_num_clusters, max_num_clusters, mean_num_clusters, median_num_clusters, 
                q3_num_clusters, q1_num_clusters, sd_num_clusters, var_num_clusters, 
                min_elapsed_time_min, max_elapsed_time_min, mean_elapsed_time_min, median_elapsed_time_min, 
                q3_elapsed_time_min, q1_elapsed_time_min, sd_elapsed_time_min, var_elapsed_time_min]
    
    column_names = ["min_tot_num_genes", "max_tot_num_genes", "mean_tot_num_genes", "median_tot_num_genes", 
                "q3_tot_num_genes", "q1_tot_num_genes", "sd_tot_num_genes", "var_tot_num_genes", 
                "min_uniq_num_genes", "max_uniq_num_genes", "mean_uniq_num_genes", "median_uniq_num_genes",
                "q3_uniq_num_genes", "q1_uniq_num_genes", "sd_uniq_num_genes", "var_uniq_num_genes", 
                "min_num_clusters", "max_num_clusters", "mean_num_clusters", "median_num_clusters", 
                "q3_num_clusters", "q1_num_clusters", "sd_num_clusters", "var_num_clusters", 
                "min_elapsed_time_min", "max_elapsed_time_min", "mean_elapsed_time_min", 
                "median_elapsed_time_min",  "q3_elapsed_time_min", "q1_elapsed_time_min", 
                "sd_elapsed_time_min", "var_elapsed_time_min"]

    temp_row = np.reshape(np.array(temp_row), (1, len(temp_row)))           
    summary_df = pd.DataFrame(temp_row, columns = column_names)
 
    return summary_df
    
    
    
    
    
    
    

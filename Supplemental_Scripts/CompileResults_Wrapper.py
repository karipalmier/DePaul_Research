# -*- coding: utf-8 -*-
"""
CompileResults_Wrapper.py

This file is used to call the CompileData.py file with different clustering result folders and files.
Note that the folder containing the intermediate csv data is hard-coded and needs to be changed to match the
actual folder name (changed folder from Temp to Intermediate during development). 

Date      Name              Comment
10/22/18  Kari Palmier      Created

"""


from CompileResults import get_cluster_statistics
from CompileResults import get_elapsed_times
from CompileResults import combine_results_times
from CompileResults import get_summary_results

base_path = "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_PBMC_4K\\20181010_22_56_36_Bicluster_Results\\"
#base_path = "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_PBMC_4K\\2018103_21_1_17_Bicluster_Results\\"
#base_path = "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_Mouse_Neurons\\20181020_11_34_42_Bicluster_Results\\"
#base_path = "C:\\DePaulCoursework\\Research\\Known_Examples\\Bryan_WebLink_Mouse_Neurons\\20181018_15_52_16_Bicluster_Results\\"

data_path = base_path + "Temp\\"
stat_df = get_cluster_statistics(data_path)

log_file = base_path + "20181010_22_56_36_Bicluster_Log.txt"
#log_file = base_path + "2018103_21_1_17_Bicluster_Log.txt"
#log_file = base_path + "20181020_11_34_42_Bicluster_Log.txt"
#log_file = base_path + "20181018_15_52_16_Bicluster_Log.txt"
time_df = get_elapsed_times(log_file)

results_df = combine_results_times(stat_df, time_df)
summary_df = get_summary_results(results_df)

# Create adjacency file name
end_ndx = log_file.find("Log.txt")
base_name = log_file[0:end_ndx]

results_path = base_name + "Results.csv"
results_df.to_csv(results_path, header = True, index = True)

summary_path = base_name + "Results_Summary.csv"
summary_df.to_csv(summary_path, header = True, index = True)



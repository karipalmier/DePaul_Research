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

base_path = "C:\\DePaulCoursework\\Research\\BryanExample\\20180627_Sample1_Infected\\2019127_14_48_8_Bicluster_Results\\"

data_path = base_path + "Temp\\"
stat_df = get_cluster_statistics(data_path)

log_file = base_path + "2019127_14_48_8_Bicluster_Log.txt"

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



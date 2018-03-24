# -*- coding: utf-8 -*-
"""

Reg_FD_Test.py

This code tests the functionality of the Create_Reg_FD_Data function.

Date      Name              Comment
2/28/18   Kari Palmier      Created

"""
import pandas as pd
import random

from Create_Reg_FD_Data import Create_Reg_FD_Data

reg_file = "C:\\DePaulCoursework\\Research\\human\\human.source"

# Example 1 - A small set of genes chosen to fulfill all group and value types
#gene_lst = ["AATF", "AC007092.1", "BTK", "C18orf54", "CALR3", "CCP110", 
#            "GM2A", "IL12RB2", "LRRC27", "NDRG2", "RPL7A", "POLR2A", "POLR2C"]
#new_json_str = Create_Reg_FD_Data(gene_lst, reg_file)
#print(new_json_str)

# Example 2 - A more realistic size list of genes
#             Genes are chosen randomly from the total list of genes in the
#             data csv file provided.  The num_genes determines how many genes
#             are chose
data_file = "C:\\DePaulCoursework\\Research\\Data\\Sample1_Infected_5000UMI.csv"
data_df = pd.read_csv(data_file, delimiter=",", index_col=0, header=0)
all_gene_lst = list(data_df.index)

num_genes = 100
gene_lst = random.sample(all_gene_lst, num_genes)
new_json_str = Create_Reg_FD_Data(gene_lst, reg_file)
print(new_json_str)


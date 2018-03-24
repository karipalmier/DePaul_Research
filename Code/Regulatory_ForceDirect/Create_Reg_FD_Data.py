# -*- coding: utf-8 -*-
"""
CreateRegForceDirect.py

This code generates an output JSON string containing the information required
to create a force direct diagram of the regulatory network of a list of genes.
The json string contains a dictionary entry for the force direct links and one
for the nodes.  Note that the entries of the nodes and links dictionaries have  
the entries in random orders.  This could not be changed due to how Python 
handles dicts.

inputs:
gene_lst = list of genes to be analyzed
reg_file = path of human.source file containing regulatory information

NOTE: An option is to pass in the dataframe of the reg information instead of
      the file name.  This will save file loading time.  Didn't do this 
      because didn't know how the calling app backbone code is laid out.

Date      Name              Comment
2/28/18   Kari Palmier      Created

"""

import pandas as pd
import numpy as np
import json

#reg_file = "C:\\DePaulCoursework\\Research\\human\\human.source"
#gene_lst = ["AATF", "AC007092.1", "BTK", "C18orf54", "CALR3", "CCP110", 
#            "GM2A", "IL12RB2", "LRRC27", "NDRG2", "RPL7A", "POLR2A", "POLR2C"]

def Create_Reg_FD_Data(gene_lst, reg_file):

    # read the regulatory info from the human.source file into a dataframe
    all_reg_df = pd.read_csv(reg_file, delimiter="\t", index_col=None, header=None)
    all_reg_df.columns = ["regulator", "reg_num", "target", "target_num"]
    
    #initialize dataframe and set to be populated in loop
    reg_df = pd.DataFrame(columns = ['source', 'target', 'value'])
    gene_set = set()
    
    # loop over number of genes in list passed in
    for i in range(len(gene_lst)):
        
        # initialize temporary dataframes
        new_reg_df = pd.DataFrame(columns = ['source', 'target', 'value'])
        new_target_df = pd.DataFrame(columns = ['source', 'target', 'value'])
        
        # get gene name and add it to the set (will be list of unique genes)
        gene = gene_lst[i]
       
        # populate a dataframe with entries that gene is the regulator
        reg_ndx = all_reg_df['regulator'] == gene
        if any(reg_ndx):
            tmp_reg_df = all_reg_df[reg_ndx]
            new_reg_df = pd.concat([tmp_reg_df['regulator'], tmp_reg_df['target']], axis=1, 
                                keys=['source', 'target'])
            new_reg_df['value'] = pd.Series(np.ones(tmp_reg_df.shape[0]), 
                      index=tmp_reg_df.index, dtype='int32')
            gene_set.add(gene)
            gene_set.update(list(new_reg_df['target']))
                        
        # populate a dataframe with entries that gene is the target
        target_ndx = all_reg_df['target'] == gene
        if any(target_ndx):
            tmp_target_df = all_reg_df[target_ndx]
            new_target_df = pd.concat([tmp_target_df['regulator'], tmp_target_df['target']], axis=1, 
                                keys=['source', 'target'])
            new_target_df['value'] = pd.Series(np.ones(tmp_target_df.shape[0]), 
                         index=tmp_target_df.index, dtype='int32')
            gene_set.add(gene)
            gene_set.update(list(new_target_df['source']))
                 
        # if the genes that the gene regulates are also regulators in the target
        # dataframe, set the counts to 2 (gene pair is present in both directions)
        if new_reg_df.shape[0] > 0 and new_target_df.shape[0] > 0:
            for j in range(len(new_reg_df['target'])):
                match_ndx = new_target_df['source'] == new_reg_df.iloc[j, 1]
                if any(match_ndx):
                    new_reg_df.iloc[j, 2] = 2
                    tmp_ndx = new_target_df[match_ndx].index
                    new_target_df.at[tmp_ndx, 'value'] = 2
    
        # add temporary dataframes to the main one
        reg_df = reg_df.append(new_reg_df)
        reg_df = reg_df.append(new_target_df)
        
    # initialize a list to be populated with group values
    gene_group = []
    total_gene_lst = list(gene_set)
    
    # loop over the number of unique genes
    for x in range(len(total_gene_lst)):
        tmp_gene = total_gene_lst[x]
        
        # if current gene is in the orig list passed in, set group to 1
        if tmp_gene in gene_lst:
            gene_group.append(1)
        # if current gene is both a regulator and target, set group to 4
        elif any(reg_df['source'] == tmp_gene) and any(reg_df['target'] == tmp_gene):
            gene_group.append(4)
        # if current gene is a regulator, set group to 2
        elif any(reg_df['source'] == tmp_gene):
            gene_group.append(2)
        # if current gene is a target (is regulated by something) set group to 3
        elif any(reg_df['target'] == tmp_gene):
            gene_group.append(3)
    
    # create dataframe of genes and groups
    node_df = pd.concat([pd.Series(total_gene_lst), pd.Series(gene_group)], axis=1, 
                        keys=['id', 'group'])
    
    # create list of dict entries of reg_df
    link_dict = reg_df.to_dict('index')
    link_lst = list(link_dict.values())
    
    # create list of dict entries of node_df
    node_dict = node_df.to_dict('index')
    node_lst = list(node_dict.values())
    
    # create dict with both links and nodes
    final_dict = {'links': link_lst}
    final_dict['nodes'] = node_lst
    
    # convert final dict to json string
    json_str = json.dumps(final_dict)
    return json_str


       


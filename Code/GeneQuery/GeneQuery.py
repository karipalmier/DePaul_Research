# -*- coding: utf-8 -*-
"""
GeneQuery.py

This code calls the R GetGeneInfo script with a string name.  The R script
returns JSON text containing all of the gene information.  This JSON text 
is then converted and parsed.

Date      Name              Comment
2/6/18    Kari Palmier      Created

"""

import subprocess
import json
import os 
# Define command and arguments
command = 'Rscript'

dir_path = os.path.dirname(os.path.realpath(__file__))
script_path = dir_path + '\\GetGeneInfo.R'

# Set argument to R script - note this will need to be passed from another
# layer in the final code
gene_str = ['POLR2A']
#gene_str = ['MTRNR2L1']
#gene_str = ['TFIIA']
#gene_str = ['MT-RNR2']

# Build subprocess command
cmd = [command, script_path] + gene_str

# check_output will run the command and store to result
gene_json = subprocess.check_output(cmd, universal_newlines=True)
gene_dict = json.loads(gene_json)

# Parse out each gene entry and save to variables
genomic_pos = gene_dict['genomic_pos']
genomic_pos_hg19 = gene_dict['genomic_pos_hg19']
type_of_gene = gene_dict['type_of_gene']
pathway = gene_dict['pathway']
generif = gene_dict['generif']
other_names = gene_dict['other_names']
name = gene_dict['name']
summary = gene_dict['summary']
wikipedia = gene_dict['wikipedia']
go = gene_dict['go']

# If pathyway entries are present, parse out individual pathway entries
try:
    biocarta_pathway = gene_dict['pathway']['biocarta']
except:
    biocarta_pathway = {}
    
try:
    kegg_pathway = gene_dict['pathway']['kegg']
except:
    kegg_pathway = {}

try:
    reactome_pathway = gene_dict['pathway']['reactome']
except:
    reactome_pathway = {}

try:
    wiki_pathway = gene_dict['pathway']['wikipathways']
except:
    wiki_pathway = {}
 
# If genomic oncology are entries are present, parse out individual GO entries
try:
    go_BP = gene_dict['go']['BP']
except:
    go_BP = {}
   
try:
    go_CC = gene_dict['go']['CC']
except:
    go_CC = {}
   
try:
    go_MF = gene_dict['go']['MF']
except:
    go_MF = {}
        
# Below are print statements to test out the functionality of the code
# These are left in the code so others can see how to access entries
"""
print('Genomic Position HG19:')
print(genomic_pos_hg19)
print('\n')
print('Kegg Pathway Entries:')
print(kegg_pathway)
print('\n')
print('Generif First Entry:')
print(generif[0])
print('\n')
print('GO CC Second Entry:')
print(go_CC[1])
print('\n')
"""


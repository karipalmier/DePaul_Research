# -*- coding: utf-8 -*-
"""
GetEnrichedGO.py

This function calls the GetEnrichedGO.R script.  This is done to test out the
use case of the R script in the Django framework.

Date      Name              Comment
5/1/18    Kari Palmier      Created

"""

import subprocess
import json

# Define command and arguments
command = 'Rscript'

#dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = "C:\\DePaulCoursework\\Research\\Code"
script_path = dir_path + '\\GetEnrichedGO.R'

# Set argument to R script - note this will need to be passed from another
# layer in the final code
#gene_list = ["POLR2A", "LDHB", "NOSIP", "AQP3", "CD79A", "CD79B"]
gene_list = ["0.05", "LDHB", "NOSIP", "AQP3", "CD79A", "CD79B"]

# Build subprocess command
cmd = [command, script_path] + gene_list

# check_output will run the command and store to result
enriched_json = subprocess.check_output(cmd, universal_newlines=True)
enriched_dict = json.loads(enriched_json)


# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 14:55:18 2023

@author: alexamin
"""
import os
import re
import numpy as np
import glob
import matplotlib.pyplot as plt

# Define the directory and pattern
directory = 'workshop-sandbox/ngr/U1341B-16H-1-W_SHLF1110331_2009730151613'
pattern = '*_SHLF*.spe'

# Form the full pattern and get the list of files
full_pattern = os.path.join(directory, pattern)
spe_files = glob.glob(full_pattern)

# Initialize ndata, data and data_info
ndata = 1024
data = np.zeros((ndata, len(spe_files)))
data_info = np.zeros((4, len(spe_files)))
t_data = np.zeros(len(spe_files))

# Iterate over the spe_files and extract the required information and data
for j, file_path in enumerate(spe_files):
    file = os.path.basename(file_path)  # Get only the file name without the directory
    with open(file_path, 'r') as fid:
        print(file)  # equivalent to spe_files(j).name in MATLAB
        
        # Extract information from file name using regular expression
        info = [int(x) for x in re.findall(r'\d+', file)]
        data_info[:, j] = np.array(info)[[2, 3, 6, 7]]
        
        lines = fid.readlines()
        for i, tline in enumerate(lines):
            if tline.strip().startswith('$MEAS_TIM:'):
                # Extract t_data from the next line
                tline = lines[i + 1].strip()
                t_data[j] = float(re.search(r"[-+]?\d*\.\d+|\d+", tline).group(0))
                
            elif tline.strip().startswith('$DATA:'):
                # Extract data from the following ndata lines
                for k in range(ndata):
                    tline = lines[i + k + 2].strip()
                    data[k, j] = float(tline)
                break  # Exit the loop once all data are extracted

plt.Figure, plt.plot(data)

print(data)
print(data_info)
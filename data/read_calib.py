# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 14:50:44 2023

@author: alexamin
"""
import numpy as np
import os
import glob
import re

# Specify the directory containing the calibration files
directory = 'workshop-sandbox/ngr/U1341B-16H-1-W_SHLF1110331_2009730151613'
# Create a list of calibration files in the specified directory
cal_files = glob.glob(os.path.join(directory, '*_calibration*.spe'))

# Initialize the ENER_FIT array
ener_fit = np.zeros((2, len(cal_files)))

# Iterate over the calibration files and extract the ENER_FIT values
for j, file in enumerate(cal_files):
    with open(file, 'r') as fid:
        lines = fid.readlines()
        for i, tline in enumerate(lines):
            # Search for the line starting with '$ENER_FIT:'
            if tline.strip().startswith('$ENER_FIT:'):
                # The values are on the next line
                tline = lines[i + 1].strip()
                temp = [float(val) for val in re.findall(r"[-+]?\d*\.\d+|\d+", tline)]
                ener_fit[:, j] = temp
                break  # Exit the inner loop once the values are found
print(ener_fit)
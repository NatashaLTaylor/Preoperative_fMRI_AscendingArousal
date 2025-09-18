#%% Multiple subjects

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 10:43:27 2022

@author: medelero, natas, JoshB
"""

from nilearn.interfaces import fmriprep
from nilearn.maskers import NiftiLabelsMasker
from nilearn.image import load_img
# from nilearn import datasets, plotting
from nilearn.plotting import plot_roi
from scipy.io import savemat
import nibabel as nib
# import matplotlib.pyplot as plt
import os

#%% Define Functions

def parcellate(fname, atlas, custom_atlas, confounds_simple, MNI_res2):
    
    img = load_img(fname) 
        
    if confounds_simple==True:
        confounds, sample_mask = fmriprep.load_confounds(
            fname, strategy=["motion", "wm_csf"], # You can add your favorite confounds!
            motion="derivatives", wm_csf="basic",demean=True)
    else:
        confounds, sample_mask = fmriprep.load_confounds(fname)

    if custom_atlas==True:
        if MNI_res2==True:
            masker = NiftiLabelsMasker(labels_img=atlas, standardize=True, verbose=5,
                                       low_pass=0.1, high_pass=0.01, t_r=2.6)
        else:
            masker = NiftiLabelsMasker(labels_img=atlas, standardize=True, verbose=5,
                                       resampling_target="labels", low_pass=0.1,
                                       high_pass=0.01, t_r=2.6)
    else:
        if MNI_res2==True:
            masker = NiftiLabelsMasker(labels_img=atlas['maps'], standardize=True,
                                       verbose=5, low_pass=0.1, high_pass=0.01,
                                       t_r=2.6)
        else:
            masker = NiftiLabelsMasker(labels_img=atlas['maps'], standardize=True,
                                       verbose=5, resampling_target="labels",
                                       low_pass=0.1, high_pass=0.01, t_r=2.6)
        
    time_series = masker.fit_transform(img, confounds=confounds)
    
    return time_series

       
#%% Run Example

# 1. Prepare Data

# Load list of subj (you can also use the .tsv list and open it with pandas)

path_to_dataset = '/scratch/IPOD_B3/V1_MRI/fmriprep_output'
outpath = '/scratch/IPOD_B3/V1_MRI/timeseries/'

# We will use html for listing subjects

subjs = [pos[:+8] for pos in os.listdir(path_to_dataset) if pos.startswith('sub-')]
print(subjs)


# List sessions

sessions = ['ses-1']

task = 'rest'

# 2. Prepare Atlas

# Load in online atlas
# atlas = datasets.fetch_atlas_schaefer_2018()
# labels = atlas['labels']

# Load in custom atlas
#atlas = nib.load('/scratch/IPOD_B3/denoise/voltron_400.nii.gz') # path to custom atlas .nii.gz
atlas = nib.load('/scratch/IPOD_B3/atlas/voltron_1000.nii.gz') # path to custom atlas .nii.gz

# Get atlas image information
#print(atlas)
# Visualise ROI's of atlas
#plot_roi(atlas)


# 3. Loop and save .mat for each subj across conds (this can take a while)

for s in subjs:
    # Loop across sessions and extract time-series for each session of each subject
    for ses in sessions:
        try:
            fname = (path_to_dataset + '/' + s  + '/'  + ses + '/func/' + s + '_' + ses +
                     '_task-' + task + '_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz') # This may vary depending on dataset names
            ts = parcellate(fname, atlas, custom_atlas=True, confounds_simple=True,
                            MNI_res2=True) 
        # And save! Saves as time x ROI
            savemat(outpath + s + '_' + ses + '_timeseries.mat', {'ts': ts})
        except:
            pass
  
    
# Check time-series
# plt.plot(ts, marker='o')


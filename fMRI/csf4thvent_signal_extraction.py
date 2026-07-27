## Script to extract mean time series from 4th ventricle masks for multiple subjects based on their individual extracted csf 4th ventricle masks and corresponding BOLD images.
import os
import glob
import numpy as np
import nibabel as nib
from nilearn.maskers import NiftiMasker
from nilearn.image import load_img
from nilearn.interfaces import fmriprep


def parcellate(fname, atlas, custom_atlas, confounds_simple, MNI_res2):
    img = load_img(fname)

    if confounds_simple is True:
        confounds, sample_mask = fmriprep.load_confounds(
            fname,
            #strategy=["motion", "wm_csf"],
            strategy=["motion"],
            motion="derivatives",
            wm_csf="basic",
            demean=True,
        )
    else:
        confounds, sample_mask = fmriprep.load_confounds(fname)

    # In this script, custom_atlas is the subject-specific CSF 4th ventricle mask.
    if custom_atlas is True:
        masker = NiftiMasker(
            mask_img=atlas,
            standardize=True,
            detrend=True,
            low_pass=0.1,
            high_pass=0.01,
            t_r=2.6,
        )
    else:
        masker = NiftiMasker(
            mask_img=atlas["maps"],
            standardize=True,
            detrend=True,
            low_pass=0.1,
            high_pass=0.01,
            t_r=2.6,
        )

    time_series = masker.fit_transform(img, confounds=confounds, sample_mask=sample_mask)
    return time_series

# ----------------------------
# PATHS (EDIT THESE)
# ----------------------------
mask_dir = "/Volumes/PRJ-IPOD_B3/fMRI_data/freesurfer_extract/outputs/csf_4thvent_mask"
fmriprep_dir = "/Volumes/PRJ-IPOD_B3/fMRI_data/V1/fmriprep_output"
output_dir = "/Users/ntaylor/Desktop/AASDelirium_timeseries_csf4thvent"

# Set to a specific subject (e.g., "sub-001") to run one subject only.
# Keep as None to run all subjects.
target_sub_id = None

os.makedirs(output_dir, exist_ok=True)

# ----------------------------
# GET SUBJECT LIST FROM FMRIPREP OUTPUT
# ----------------------------
subject_dirs = sorted(glob.glob(os.path.join(fmriprep_dir, "sub-*")))

if target_sub_id is not None:
    subject_dirs = [
        subject_dir
        for subject_dir in subject_dirs
        if os.path.basename(subject_dir) == target_sub_id
    ]
    if len(subject_dirs) == 0:
        raise ValueError(f"Requested subject not found in fmriprep_dir: {target_sub_id}")

# ----------------------------
# LOOP THROUGH SUBJECTS
# ----------------------------
for subject_dir in subject_dirs:
    sub_id = os.path.basename(subject_dir)

    # Find corresponding 4th ventricle mask for this subject
    #sub-001_csf_4thvent_mni152.nii.gz
    mask_pattern = os.path.join(mask_dir, f"{sub_id}_csf_4thvent_mni152.nii.gz")
    mask_files = glob.glob(mask_pattern)
    if len(mask_files) == 0:
        print(f"No mask file found for {sub_id}")
        continue
    mask_file = mask_files[0]

    # Find corresponding BOLD file
    bold_dir = os.path.join(fmriprep_dir, f"{sub_id}/ses-1/func")
    #sub-001_ses-1_task-rest_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz
    bold_pattern = os.path.join(bold_dir, f"{sub_id}*_ses-1_task-rest_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz")
    bold_files = glob.glob(bold_pattern)

    if len(bold_files) == 0:
        print(f"No BOLD file found for {sub_id}")
        continue
    bold_file = bold_files[0]
    
    print(f"Processing {sub_id}")
    
    # ----------------------------
    # LOAD AND VERIFY MASK
    # ----------------------------
    mask_img = nib.load(mask_file)
    mask_data = mask_img.get_fdata()
    unique_vals = np.unique(mask_data)
    
    # Check if mask is binarized (only contains 0 and 1)
    is_binary = np.all(np.isin(unique_vals, [0, 1]))
    if not is_binary:
        print(f"WARNING: Mask for {sub_id} is not binarized. Unique values: {unique_vals}")
    else:
        print(f"  Mask verified: binarized (unique values: {unique_vals})")
    
    # ----------------------------
    # EXTRACT TIME SERIES
    # ----------------------------
    # Output shape: (n_timepoints, n_voxels_in_mask)
    time_series = parcellate(
        fname=bold_file,
        atlas=mask_img,
        custom_atlas=True,
        confounds_simple=True,
        MNI_res2=True,
    )
    
    # Average across voxels → single CSF signal
    mean_ts = np.mean(time_series, axis=1)
    
    # ----------------------------
    # SAVE OUTPUT
    # ----------------------------
    out_file = os.path.join(output_dir, f"{sub_id}_4thvent_timeseries.txt")
    np.savetxt(out_file, mean_ts)
    
    print(f"Saved: {out_file}")

print("Done.")
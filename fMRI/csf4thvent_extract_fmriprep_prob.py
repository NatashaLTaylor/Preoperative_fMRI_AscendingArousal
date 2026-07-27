## Script to extract mean time series from 4th Ventricle & Lateral Ventricle (Harvard-Oxford Atlas) for multiple subjects based on their individual csf fmriprep-probability atlas and corresponding BOLD images.
import os
import glob
import numpy as np
import nibabel as nib
from nilearn.maskers import NiftiMasker
from nilearn.image import load_img, math_img
from nilearn.interfaces import fmriprep
from nilearn import datasets
# ----------------------------
# PATHS (EDIT THESE)
# ----------------------------
fmriprep_dir = "/Volumes/PRJ-IPOD_B3/fMRI_data/V1/fmriprep_output"
output_dir = "/Users/ntaylor/Desktop/AASDelirium_timeseries_csf4thvent"
mask_dir = "/Volumes/PRJ-IPOD_B3/fMRI_data/freesurfer_extract/outputs/csf_4thvent_mask"

target_sub_id = None  # Set to a specific subject (e.g., "sub-001") to run one subject only. Keep as None to run all subjects.
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

    csf_dir= os.path.join(subject_dir, "ses-1/anat")
    # Find corresponding 4th ventricle mask for this subject
    csf_prob = os.path.join(csf_dir, f"{sub_id}_ses-1_space-MNI152NLin6Asym_res-2_label-CSF_probseg.nii.gz")
    csf_files = glob.glob(csf_prob)
    if len(csf_files) == 0:
        print(f"No CSF file found for {sub_id}")
        continue
    csf_file = csf_files[0]
    csf_img = nib.load(csf_file)

    # Find corresponding 4th ventricle mask for this subject
    # sub-001_csf_4thvent_mni152.nii.gz
    mask_pattern = os.path.join(mask_dir, f"{sub_id}_csf_4thvent_mni152.nii.gz")
    mask_files = glob.glob(mask_pattern)
    if len(mask_files) == 0:
        print(f"No mask file found for {sub_id}")
        continue
    mask_file = mask_files[0]
    mask_img = nib.load(mask_file)
    #-----------------------------
    # GENERATE CSF MASK
    #------------------------------
    # only csf_img > 90% probability
    csf_mask = math_img("img > 0.9", img=csf_img)
    
    # Intersect CSF probability with 4th ventricle mask defined from freesurfer
    csf_prob_4thvent = math_img("img1 * (img2 > 0)", img1=csf_mask, img2=mask_img)
    # print out of number of voxels in the final mask
    csf_prob_4thvent_data = csf_prob_4thvent.get_fdata()
    num_voxels = np.sum(csf_prob_4thvent_data > 0)
    print(f"{sub_id}: Number of voxels in CSF 4th ventricle mask: {num_voxels}")
    
    # Intersect CSF probability with 4th ventricle region
    # (csf_prob_4thvent is already defined above)
    nib.save(csf_prob_4thvent, os.path.join(output_dir, f"{sub_id}_csf_4thvent_mask.nii.gz"))
    
    #atlas = datasets.fetch_atlas_harvard_oxford('sub-maxprob-thr0-2mm')
    #atlas_labels = atlas["labels"]
    #atlas_maps = atlas["maps"] 
    #print(atlas.labels)
    #select for the fourth ventricles
    #left_lat_vent = atlas_maps == atlas_labels.index('Left Lateral Ventricle') + 1  # +1 because atlas labels are 1-indexed
    #right_lat_vent = atlas_maps == atlas_labels.index('Right Lateral Ventricle') + 2
    #combine lat vent into one mask
    #lat_vent_mask = left_lat_vent | right_lat_vent
 
    #-----------------------------
    # TIMESERIES EXTRACTION - BOLD SIGNAL
    #------------------------------
    
    # Find corresponding BOLD file
    bold_dir = os.path.join(fmriprep_dir, f"{sub_id}/ses-1/func")
    #sub-001_ses-1_task-rest_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz
    bold_pattern = os.path.join(bold_dir, f"{sub_id}*_ses-1_task-rest_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz")
    bold_files = glob.glob(bold_pattern)

    if len(bold_files) == 0:
        print(f"No BOLD file found for {sub_id}")
        continue
    bold_file = bold_files[0]
    
    out_file = os.path.join(output_dir, f"{sub_id}_ses-1_csf4thvent_timeseries.mat")
    if os.path.exists(out_file):
        print(f"Skipping {sub_id}: output already exists ({out_file})")
        continue

    print(f"Processing {sub_id}")

    # Extract timeseries using the subject-specific CSF 4th ventricle mask
    img = load_img(bold_file)
    confounds, _ = fmriprep.load_confounds(
        bold_file,
        strategy=["motion"],
        motion="derivatives",
        demean=True,
    )
    masker = NiftiMasker(
        mask_img=csf_prob_4thvent,
        standardize=True,
        detrend=True,
        low_pass=0.1,
        high_pass=0.01,
        t_r=2.6,
    )
    ts = masker.fit_transform(img, confounds=confounds)

    # Average across voxels -> single timeseries vector
    mean_ts = np.mean(ts, axis=1)

    from scipy.io import savemat
    savemat(out_file, {"ts": mean_ts})
    print(f"Saved: {out_file}")

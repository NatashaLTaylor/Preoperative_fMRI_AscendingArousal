#!/bin/bash
#==============================================================================
# Extract GM volumes for AAS nuclei in NATIVE T1w space
#==============================================================================

#------------------------------------------------------------------------------
# PATHS
#------------------------------------------------------------------------------
DATA_PATH=/Volumes/PRJ-IPOD_B3/fMRI_data/V1/fmriprep_output
OUTPUT_PATH=~/onedrive_usyd/Postdoc_Rob/Analysis/LC_volume/1mm_native
ATLAS_PATH=~/onedrive_usyd/Postdoc_Rob/code
ATLAS="voltron_400.nii.gz"

LABELS_FILE=${OUTPUT_PATH}/labels_aas_nuclei.txt
SUBJECT_LIST=${OUTPUT_PATH}/subject_ids.txt
FAILED_LOG=${OUTPUT_PATH}/failed_subjects_aasnuclei.csv

#------------------------------------------------------------------------------
# SETUP
#------------------------------------------------------------------------------
cd "$OUTPUT_PATH" || exit 1

# Initialise conda so antsApplyTransforms is available
eval "$(conda shell.bash hook)"
conda activate ANTs

# Verify ANTs is accessible
if ! command -v antsApplyTransforms &> /dev/null; then
    echo "ERROR: antsApplyTransforms not found. Check conda environment."
    exit 1
fi

# Initialise failure log
if [ ! -f "$FAILED_LOG" ]; then
  echo "SubjectID,Reason" > "$FAILED_LOG"
fi

#------------------------------------------------------------------------------
# MAIN LOOP
#------------------------------------------------------------------------------
while IFS= read -r SUBJECT_ID; do
  # Strip hidden characters (Windows line endings, trailing spaces)
  SUBJECT_ID=$(echo "$SUBJECT_ID" | tr -d '\r\n ')

  echo "=================================================================="
  echo "Processing: $SUBJECT_ID"
  echo "=================================================================="

  #--- Define paths for this subject ---
  ANAT_PATH=${DATA_PATH}/${SUBJECT_ID}/ses-1/anat
  NATIVE_T1W=${ANAT_PATH}/${SUBJECT_ID}_ses-1_desc-preproc_T1w.nii.gz
  GM_IMAGE=${ANAT_PATH}/${SUBJECT_ID}_ses-1_label-GM_probseg.nii.gz
  INVERSE_XFM=${ANAT_PATH}/${SUBJECT_ID}_ses-1_from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5
  ATLAS_NATIVE=${OUTPUT_PATH}/${SUBJECT_ID}_atlas_native.nii.gz
  OUTPUT_CSV=${OUTPUT_PATH}/${SUBJECT_ID}_GM_aas_nuclei_volumes.csv

  #--- Skip if already processed ---
  if [ -f "$OUTPUT_CSV" ]; then
    echo "  Output CSV already exists. Skipping..."
    continue
  fi

  #--- Check all required files exist ---
  MISSING=0
  for FILE in "$NATIVE_T1W" "$GM_IMAGE" "$INVERSE_XFM"; do
    if [ ! -f "$FILE" ]; then
      echo "  MISSING: $FILE"
      MISSING=1
    fi
  done
  if [ $MISSING -eq 1 ]; then
    echo "  Skipping $SUBJECT_ID due to missing files."
    echo "${SUBJECT_ID},missing input files" >> "$FAILED_LOG"
    continue
  fi

  #--- Step 1: Warp atlas to native space ---
  echo "  Warping atlas to native space..."
  antsApplyTransforms \
    -d 3 \
    -i "${ATLAS_PATH}/${ATLAS}" \
    -r "$NATIVE_T1W" \
    -t "$INVERSE_XFM" \
    -n NearestNeighbor \
    -o "$ATLAS_NATIVE"

  if [ $? -ne 0 ]; then
    echo "  ERROR: antsApplyTransforms failed for $SUBJECT_ID"
    echo "${SUBJECT_ID},antsApplyTransforms failed" >> "$FAILED_LOG"
    continue
  fi

  #--- Step 2: Get voxel dimensions from native GM image ---
  VOXEL_VOL=$(fslval "$GM_IMAGE" pixdim1)
  DIM2=$(fslval "$GM_IMAGE" pixdim2)
  DIM3=$(fslval "$GM_IMAGE" pixdim3)
  VOXEL_VOL=$(echo "$VOXEL_VOL * $DIM2 * $DIM3" | bc -l)
  echo "  Voxel volume: $VOXEL_VOL mm3"

  #--- Step 3: Write CSV header ---
  echo "Label,Region,MeanGM_NonZero,MeanGM_AllVoxels,MaskVoxels,ProbGM_Volume_mm3,BinaryMask_Volume_mm3" > "$OUTPUT_CSV"

  #--- Step 4: Loop through each label ---
  while IFS=, read -r LABEL REGION; do
    LABEL_NUMBER=$(echo "$LABEL" | awk '{print $1}')
    TMP_MASK=/tmp/tmp_mask_${SUBJECT_ID}_${LABEL_NUMBER}.nii.gz

    # Create binary mask for this label from the WARPED atlas
    fslmaths "$ATLAS_NATIVE" -thr "$LABEL_NUMBER" -uthr "$LABEL_NUMBER" -bin "$TMP_MASK"

    # Mean GM probability across non-zero voxels only
    STATS=$(fslstats "$GM_IMAGE" -k "$TMP_MASK" -M -V)
    MEAN_NONZERO=$(echo $STATS | awk '{print $1}')
    VOXELS=$(echo $STATS | awk '{print $2}')
    VOL=$(echo $STATS | awk '{print $3}')

    # Mean GM probability across ALL voxels in the mask (including zeros)
    MEAN_ALL=$(fslstats "$GM_IMAGE" -k "$TMP_MASK" -m)

    # Number of voxels in the binary mask
    NUM_MASK_VOXELS=$(fslstats "$TMP_MASK" -V | awk '{print $1}')

    # Probabilistic GM volume = mean_all * num_mask_voxels * voxel_volume
    PROB_GM_VOL=$(echo "$MEAN_ALL * $NUM_MASK_VOXELS * $VOXEL_VOL" | bc -l)

    # Binary mask volume (atlas region volume regardless of GM)
    BINARY_VOL=$(echo "$NUM_MASK_VOXELS * $VOXEL_VOL" | bc -l)

    # Append to CSV
    echo "$LABEL_NUMBER,$REGION,$MEAN_NONZERO,$MEAN_ALL,$NUM_MASK_VOXELS,$PROB_GM_VOL,$BINARY_VOL" >> "$OUTPUT_CSV"

    # Clean up temp mask
    rm -f "$TMP_MASK"

  done < "$LABELS_FILE"

  echo "  Done. Output: $OUTPUT_CSV"

done < "$SUBJECT_LIST"

echo "=================================================================="
echo "ALL SUBJECTS COMPLETE"
echo "=================================================================="

#!/bin/bash

# Set paths
INPUT_DIR="../../../processed_data/cellbender"
OUTPUT_DIR_BASE="../../../processed_data/cellbender"
SAMPLES_FILE="${INPUT_DIR}/samples.txt"

echo "-----------printing samples file------------"
cat "$SAMPLES_FILE"
echo "--------------------------------------------"

# Loop over each sample name
while read -r sample_name; do
  echo "Converting $sample_name"

  # choose FPR 0.0 which removes no cells. empty droplets will be removed manually downstream.
  input_h5="${INPUT_DIR}/${sample_name}/cellbender_output_FPR_0.0_filtered.h5" 
  output_dir="${OUTPUT_DIR_BASE}/${sample_name}"
  output_h5="${output_dir}/cellbender_output_FPR_0.0_filtered_seurat.h5"

  # Check if the input file exists
  if [[ ! -f "$input_h5" ]]; then
    echo "Warning: Input file not found for $sample_name at $input_h5. Skipping..."
    continue
  fi

  ptrepack --complevel 5 ${input_h5}:/matrix ${output_h5}:/matrix

    
done < "$SAMPLES_FILE"
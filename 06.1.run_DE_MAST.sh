#!/bin/bash

# Specify the path to your .txt file
TXT_FILE="DE.idents.txt"

# Check if the file exists
if [ -e "$TXT_FILE" ]; then
  # Read the file line by line using a for loop
  while IFS=$'\t' read -r col1 col2; do
    # Process each line by assigning values to variables and feeding them to cal.R
    echo "Processing values: ident.1=$col1, ident.2=$col2"
    
    # Example: Run cal.R with the extracted values
    sbatch DE_MAST_sbatch.csh "$col1" "$col2"
  done < "$TXT_FILE"
else
  echo "Error: File $TXT_FILE not found."
fi


#!/bin/sh

#SBATCH
#SBATCH --job-name=DEMAST
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=norm
#SBATCH --mem=200g
#SBATCH --gres=lscratch:500

ml R/4.3.0

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 arg1 arg2"
  exit 1
fi

# Assign command-line arguments to variables
arg1="$1"
arg2="$2"

echo "$arg1 and $arg2 are inputs to R scripts"
echo "starting DE_MAST.R"
Rscript DE_MAST.R "$arg1" "$arg2"

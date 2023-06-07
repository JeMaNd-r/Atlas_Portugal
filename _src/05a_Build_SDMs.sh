#!/bin/bash

#SBATCH --job-name=05a_Build_SDMs
#SBATCH --output=/work/%u/%x-%A-%a.log
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-type=BEGIN,END,FAIL

# module load foss/2020b R/4.0.4-2
module load foss/2022b R/4.2.2

taxon_name="$1"
species_csv="_intermediates/SDM_$taxon_name.csv"

output_dir="/work/$USER/$SLURM_JOB_NAMEs-$SLURM_JOB_ID"
mkdir -p "$output_dir"


Rscript _src/05a_Build_SDMs.R "$taxon_name" "$species_csv" "$output_dir"

# cd ~/Atlas_Portugal
# sbatch -a 1-$(xsv count _intermediates/SDM_earthworms.csv) _src/05a_Build_SDMs.sh earthworms 




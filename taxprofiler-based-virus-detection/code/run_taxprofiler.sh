#!/bin/bash
#SBATCH --job-name=TaxProfiler
#SBATCH --partition=bahl_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=TaxProfiler_%j.out
#SBATCH --error=TaxProfiler_%j.err

cd "$SLURM_SUBMIT_DIR"

# Check if required files exist
echo "Checking configuration files..."
if [ ! -f "samplesheet.csv" ]; then
    echo "Error: samplesheet.csv file does not exist"
    exit 1
fi

if [ ! -f "databases.csv" ]; then
    echo "Error: databases.csv file does not exist"
    exit 1
fi

echo "Configuration file check complete"

# Load conda environment (only needed for Nextflow itself)
module load Miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nextflow_env

# Set Singularity cache directory to avoid warnings
export NXF_SINGULARITY_CACHEDIR="/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/singularity_cache"

# Set Singularity bind paths
export SINGULARITY_BIND="/scratch/sp96859/Meta-genome-data-analysis/databases:/databases,\
/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/data:/data,\
/scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/results:/results"

# Run Taxprofiler with explicit revision
nextflow run nf-core/taxprofiler -r 1.2.4 \
  -profile singularity \
  --input /scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/samplesheet.csv \
  --databases /scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/databases.csv \
  --outdir /scratch/sp96859/Meta-genome-data-analysis/Apptainer/taxprofiler/results \
  --run_kraken2


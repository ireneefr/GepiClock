#!/bin/bash
#SBATCH --job-name=umap_tcga
#SBATCH --mail-type=ALL
#SBATCH --mail-user=irene.fernandez@fht.org
#SBATCH --partition=cpuq
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=umap_tcga%j.out.log
#SBATCH --error=umap_tcga%j.err.log
#SBATCH --mem=250G

######################################################################################################################
### Set the environment
###################################################################################################################### 
module load R/4.1.0

Rscript /group/iorio/Irene/epiclock_dev/src/run_umap_tcga.R

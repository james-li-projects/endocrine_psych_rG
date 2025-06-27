#!/bin/bash
#SBATCH --job-name=compute_rg_only_hormone
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/xtrait/logs/compute_rg_only_hormone_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/xtrait/logs/compute_rg_only_hormone_%A.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=25gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load miniconda3/23.1.0

source activate ldsc 

thormone1=${ARGS1}
thormone2=${ARGS2}

# munge sumstats
/gpfs/data/huo-lab/BCAC/james.li/conda_env/ldsc/bin/python /gpfs/data/huo-lab/BCAC/james.li/preconfluence/software/ldsc/ldsc.py \
    --rg /gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_thormone/${thormone1},/gpfs/data/huo-lab/BCAC/james.li/xtrait/input/munge_thormone/${thormone2} \
    --ref-ld-chr /gpfs/data/huo-lab/BCAC/james.li/thormone/EUR_hm3_corrected_ldsc/ \
    --w-ld-chr /gpfs/data/huo-lab/BCAC/james.li/thormone/EUR_hm3_corrected_ldsc/ \
    --out /gpfs/data/huo-lab/BCAC/james.li/xtrait/output/rG_only_hormone/${thormone1}_${thormone2}

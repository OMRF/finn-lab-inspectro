#!/usr/bin/bash -l 
#SBATCH --cpus-per-task 4
#SBATCH --mem 64G
#SBATCH --mail-type END,FAIL

# check if in slurm first
if [ -z "${SLURM_JOB_ID}" ]; then
    echo "Please run this using sbatch"
    exit
fi

ml python
input=$1
command="import pandas; pandas.read_parquet('$input').to_csv('$input.csv', index=False)"
python -c "$command"

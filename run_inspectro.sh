ml inspectro
ml slurm
sbatch --cpus-per-task 20 --mem 248G --wrap 'snakemake --cores all --snakefile Snakefile'

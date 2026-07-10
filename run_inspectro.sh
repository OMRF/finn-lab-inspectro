ml inspectro
ml slurm
sbatch --job-name inspectro --cpus-per-task 20 --mem 248G --wrap 'snakemake --cores all --snakefile Snakefile'

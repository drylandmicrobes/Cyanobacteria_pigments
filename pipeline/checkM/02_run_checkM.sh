#!/usr/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 4 --out logs/checkm_public.%a.%A.log

module load checkm
module unload miniconda2/4.4.10
module load anaconda3/4.5.4

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then            
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

checkm taxonomy_wf phylum Cyanobacteria -t 8 -x .pep.fasta --genes Public_genomes/pep/ checkM_public/checkM/ -f checkM_taxonomy_wf_results_public --tab_table


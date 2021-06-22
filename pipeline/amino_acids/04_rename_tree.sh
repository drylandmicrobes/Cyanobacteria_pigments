#!/usr/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 4 --out logs/phmmer.tre_rename.%a.%A.log
module unload miniconda2
module load perl

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
	echo "need to provide a number by --array or cmdline"
	exit
    fi
fi

GENEID=geneIDs_pHMMer
DAT=carot_genes.txt
IFS=,

tail -n +2 $DAT | sed -n ${N}p | while read GENE

do
    perl scripts/rename_tree_nodes.pl genetree_pHMMer/${GENE}.tre prefixes_pHMMer/prefix.${GENE}.tab > genetree_pHMMer/${GENE}_long.tre
done

#!/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 4 --out logs/phmmer.%a.%A.log

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

DAT=lib/ncbi_accessions.csv
DBFOLDER=$(realpath databases)
IFS=,
tail -n +2 $DAT | sed -n ${N}p | while read SPECIES STRAIN ACCESSION ASMNAME BIOPROJECT SRA HABITAT LOCATION
do
    OUTNAME=$(echo "${SPECIES}_$STRAIN" | perl -p -e 's/\s+$//; s/[\s\/\;]/_/g')
    echo "Processing $OUTNAME"


cd /pHMMer_run && mv $OUTNAME* $OUTNAME 

#!/usr/bin/bash
#SBATCH -p short --mem 16gb -N 1 -n 4 --out logs/autometa_run.%a.%A.log


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
AUTOMETAOUT=autometa_runs
COVERAGE=coverage
DAT=lib/ncbi_accessions.csv
IFS=,
tail -n +2 $DAT | sed -n ${N}p | while read SPECIES STRAIN ACCESSION ASMNAME BIOPROJECT SRA HABITAT LOCATION
do
    OUTNAME=$(echo "${SPECIES}_$STRAIN" | perl -p -e 's/\s+$//; s/[\s\/\;]/_/g')
    echo "Processing $OUTNAME"
#    GENOMEFILE=genomes/${OUTNAME}.dna.fasta
    COVTAB=$COVERAGE/$OUTNAME.coverage.tab

    if [ -d $AUTOMETAOUT/$OUTNAME ]; then
	module load autometa/1.0.2
	module load git
	time run_autometa.py -k bacteria -a $AUTOMETAOUT/$OUTNAME/Bacteria.fasta --ML_recruitment \
		--processors $CPU --length_cutoff 1500 --taxonomy_table $AUTOMETAOUT/$OUTNAME/taxonomy.tab -o $AUTOMETAOUT/$OUTNAME/Bacteria_run -v $COVTAB
    fi
done

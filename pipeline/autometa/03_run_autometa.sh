#!/usr/bin/bash
#SBATCH -p batch --mem 16gb -N 1 -n 24 --out autometa_run.%a.log



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
IFS=,
tail -n +2 $DAT | sed -n ${N}p | while read SPECIES STRAIN ACCESSION ASMNAME BIOPROJECT SRA HABITAT LOCATION
do
    OUTNAME=$(echo "${SPECIES}_$STRAIN" | perl -p -e 's/\s+$//; s/[\s\/\;]/_/g')
    echo "Processing $OUTNAME"
    OUT=coverage
#    GENOMEFILE=genomes/${OUTNAME}.dna.fasta
    COVTAB=$(realpath $OUT/$OUTNAME.coverage.tab)

    if [ -d $OUTNAME.autometa ]; then
	pushd $OUTNAME.autometa

	module load autometa/1.0.2
	source activate autometa
	module load git

	run_autometa.py -k bacteria -a Bacteria.fasta --ML_recruitment \
	    --processors $CPU --length_cutoff 1500 --taxonomy_table taxonomy.tab -o Bacteria_run -v $COVTAB
    fi
done

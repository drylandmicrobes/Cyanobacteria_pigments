#!/usr/bin/bash
#SBATCH -p batch,intel --mem 32gb -N 1 -n 24 --out logs/autometa_taxonomy.%a.%A.log

# see module load below

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi


N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
fi

if [ -z $N ]; then
    echo "need to provide a number by --array or cmdline"
    exit
fi
DATABASES=databases
if [ ! -d $DATABASES ]; then
    ln -s /srv/projects/db/autometa/1.0.2 $DATABASES
fi
DAT=lib/ncbi_accessions.csv
IFS=,
tail -n +2 $DAT | sed -n ${N}p | while read SPECIES STRAIN ACCESSION ASMNAME BIOPROJECT SRA HABITAT LOCATION
do
    OUTNAME=$(echo "${SPECIES}_$STRAIN" | perl -p -e 's/\s+$//; s/[\s\/\;]/_/g')
    echo $OUTNAME
    GENOMEFILE=genomes/${OUTNAME}.dna.fasta
    COVTAB=coverage/$OUTNAME.coverage.tab
    AUTOMETAOUT=autometa_runs
    mkdir -p $AUTOMETAOUT
    if [ ! -f $COVTAB ]; then
	bash pipeline/autometa/00_make_cov.sh $N
    fi
    if [ ! -d $AUTOMETAOUT/$OUTNAME ]; then
	module load autometa/1.0.2
	make_taxonomy_table.py -a $GENOMEFILE -p $CPU -o $AUTOMETAOUT/$OUTNAME --cov_table $COVTAB
    fi
done

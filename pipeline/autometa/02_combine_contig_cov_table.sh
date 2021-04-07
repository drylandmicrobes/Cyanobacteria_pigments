#!/usr/bin/bash
#SBATCH -N 1 -n 1 -p short --out logs/contig_cov_table.%a.log

module load autometa/1.0.2
source activate autometa
module load git

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
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
    echo $OUTNAME
    OUT=coverage
    GENOMEFILE=genomes/${OUTNAME}.dna.fasta
    COVTAB=$OUT/$OUTNAME.coverage.tab
    CONTIGCOV=$OUT/$OUTNAME.contig_cov_table.tab
    if [ ! -f $CONTIGCOV ]; then
	make_contig_table.py -a $BASE -c $COVTAB -o $$CONTIGCOV
    fi
done

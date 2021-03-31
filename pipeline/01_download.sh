#!/usr/bin/bash
#SBATCH -p short --out logs/download.log
module load aspera
DAT=lib/ncbi_accessions.csv
OUT=source/NCBI_ASM
STAGE=genomes
mkdir -p $STAGE $OUT
IFS=,
tail -n +2 $DAT | while read SPECIES STRAIN ACCESSION ASMNAME BIOPROJECT SRA HABITAT LOCATION
do
    OUTNAME=$(echo "${SPECIES}_$STRAIN" | perl -p -e 's/\s+$//; s/[\s\/\;]/_/g')
    PRE=$(echo $ACCESSION | cut -d_ -f1 )
    ONE=$(echo $ACCESSION | cut -d_ -f2 | awk '{print substr($1,1,3)}')
    TWO=$(echo $ACCESSION | cut -d_ -f2 | awk '{print substr($1,4,3)}')
    THREE=$(echo $ACCESSION | cut -d_ -f2 | awk '{print substr($1,7,3)}')
    ASMNAME=$(echo $ASMNAME | perl -p -e 's/ /_/g')
    echo "anonftp@ftp.ncbi.nlm.nih.gov:/genomes/all/$PRE/$ONE/$TWO/$THREE/${ACCESSION}_${ASMNAME}/"
    if [ ! -d $OUT/${ACCESSION}_$ASMNAME ]; then
      ascp -k1 -Tdr -l400M -i $ASPERAKEY --overwrite=diff anonftp@ftp.ncbi.nlm.nih.gov:/genomes/all/$PRE/$ONE/$TWO/$THREE/${ACCESSION}_$ASMNAME ./$OUT/
    fi
    if [ ! -s $STAGE/${OUTNAME}.dna.fasta ]; then
      pigz -dc $OUT/${ACCESSION}_${ASMNAME}/${ACCESSION}_${ASMNAME}_genomic.fna.gz > $STAGE/${OUTNAME}.dna.fasta
    fi
done

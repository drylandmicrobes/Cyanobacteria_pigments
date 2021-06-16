#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 --mem 2gb --out logs/download.log

#CYANOBASELIST=https://genome.microbedb.jp/cyanobase.csv
#mkdir -p lib
# kind of slow so not re-downloading unless the file is removed or empty
#if [ ! -s lib/cyanobase.csv ]; then
#	curl -o lib/cyanobase.csv -C - -L $CYANOBASELIST
#fi
INLIST=lib/Cyanobase_genomes.csv
NCBI=lib/Cyanobase_genomes.ncbi.json
NCBICSV=$(echo $NCBI | perl -p -e 's/\.json/.csv/')
BINDIR=bin
mkdir -p $BINDIR

if [ ! -f $BINDIR/datasets ]; then
	curl -L -C - -o $BINDIR/dataformat https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/dataformat
	curl -L -C - -o $BINDIR/datasets https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets

	chmod +x $BINDIR/dataformat $BINDIR/datasets
fi


if [ ! -s $NCBI ]; then
	cut -d, -f2 $INLIST | grep -v Assembly > assemblies.$$
	$BINDIR/datasets summary genome accession  --inputfile assemblies.$$ > $NCBI
	unlink assemblies.$$
fi

if [ ! -s $NCBICSV ]; then
	./scripts/assembly_json_process.py --infile $NCBI --outfile $NCBICSV
	perl -i -p -e 's/\r\n/\n/g; s/\r/\n/g;' $NCBICSV
fi

module load aspera
# one file for now but there may be more
DAT=(lib/Cyanobase_genomes.ncbi.csv)
OUT=../source/NCBI_ASM
STAGE=../genomes
mkdir -p $STAGE $OUT
IFS=,
for file in ${DAT[@]};
do
  echo "file is $file"
  tail -n +2 $file | while read -r ACCESSION SPECIES STRAIN TAXID BIOPROJECT ASMLEN N50 ASMNAME
  do
	  echo "acc=$ACCESSION sp=$SPECIES strain=$STRAIN name=$ASMNAME"
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
done

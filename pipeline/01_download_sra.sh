#!/usr/bin/bash
#SBATCH -p batch,intel -N 1 -n 2 --mem 2gb --time 48:00:00 --out logs/download_sra.log
module load aspera
module load sratoolkit
DAT=lib/ncbi_accessions.csv
OUT=source/NCBI_SRA
mkdir -p $OUT
IFS=,

tail -n +2 $DAT | cut -f6 -d, | while read SRA
do
#	URL=$(echo $SRA | perl scripts/sra2url.pl --ftp )
#	echo "ascp -QT -l 1000m -P33001 -i "$ASPERAKEY" $URL $OUT/"
	if [ ! -f $OUT/${SRA}_2.fastq.gz ]; then
		fastq-dump --gzip --split-e -O $OUT $SRA
	fi
#	for DIRECTION in 1 2;
#	do
#		if [ ! -f $OUT/${SRA}_${DIRECTION}.fastq.gz ];
#			ascp -k1 -Tdr --overwrite=diff -QT -l 1000m -P33001 -i "$ASPERAKEY" $URL/${SRA}_${DIRECTION}.fastq.gz $OUT/${SRA}_${DIRECTION}.fastq.gz
#			if [ ! -f $OUT/${SRA}_${DIRECTION}.fastq.gz ];

				#curl -C- -L -o $OUT/${SRA}_${DIRECTION}.fastq.gz $URL/${SRA}_${DIRECTION}.fastq.gz
#			fi
#	done
	break
done

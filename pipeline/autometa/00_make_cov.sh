#!/usr/bin/bash
#SBATCH -N 1 -n 16 -p batch --mem 16gb --out logs/make_cov.%a.log

module load bwa
module load samtools/1.11
module load bedtools

N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi


if [ -z $N ]; then
 N=$1
fi

if [ -z $N ]; then
 echo "need to provide a number by --array or cmdline"
 exit
fi

DAT=lib/ncbi_accessions.csv
IFS=,
tail -n +2 $DAT | sed -n ${N}p | while read SPECIES STRAIN ACCESSION ASMNAME BIOPROJECT SRA HABITAT LOCATION
do
    OUTNAME=$(echo "${SPECIES}_$STRAIN" | perl -p -e 's/\s+$//; s/[\s\/\;]/_/g')
    echo "Processing $OUTNAME"
    OUT=coverage
    mkdir -p $OUT
    GENOMEFILE=genomes/${OUTNAME}.dna.fasta
    # small speedup would be to write this to /scratch instead of the current directory
    BAM=/scratch/$OUTNAME.remap.bam
    COV=$OUTNAME.cov
    COVTAB=$OUT/$OUTNAME.coverage.tab
    FWD=source/NCBI_SRA/${SRA}_1.fastq.gz
    REV=source/NCBI_SRA/${SRA}_2.fastq.gz
    if [ ! -f $COVTAB ]; then
	if [ ! -f $GENOMEFILE.bwt ]; then
	    bwa index $GENOMEFILE
	fi
	bwa mem -t $CPU $GENOMEFILE $FWD $REV | samtools sort --threads $CPU -T /scratch -O bam -o $BAM -
	samtools index $BAM
	
	# can replace this also with samtools faidx and a cut cmd
	#fasta_length_table.pl $ASSEMBLY > $BASE.genome.lengths
	genomeCoverageBed -ibam $BAM  > $OUT/$OUTNAME.genome_cov.bed


	module load autometa/1.0.2
	source activate autometa

	contig_coverage_from_bedtools.pl $OUT/$OUTNAME.genome_cov.bed > $COVTAB
	
	rm -f $BAM $BAM.bai
    fi
done

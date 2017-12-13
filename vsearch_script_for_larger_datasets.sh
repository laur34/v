#! /bin/bash

## VSEARCH-based pipeline for Illumina paired-end metabarcoding reads
## 12.2017
## Run by typing bash vsearch_script_for_larger_datasets.sh

##Define Variables (Change to your primer sequences, these are MiniLGC, threads in computer, and directories):
#PRIMER_F="TCATGGWACWGGWTGAACWGTWTAYCCYCC"
#PRIMER_R_RC="TGRTTYTTTGGTCACCCTGAAGTTTACTCC"
#PRIMER_R="GGAGTAAACTTCAGGGTGACCAAARAAYCA"
#PRIMER_F_RC="GGRGGRTAWACWGTTCAWCCWGTWCCATGA"

PRIMER_F="GGWACWGGWTGAACWGTWTAYCCYCC"
PRIMER_R_RC="TGRTTYTTTGGTCACCCTGAAGTTTACTCC"
PRIMER_R="TAAACTTCAGGGTGACCAAARAAYCA"
PRIMER_F_RC="GGRGGRTAWACWGTTCAWCCWGTWCCATGA"


MINLEN=300

THREADS=8
VSEARCH=$(which vsearch2)

#DATADIR="/media/laur/INTENSO/AIMSEQ10082017 - RIMA/RAW DATA/RAW"
#OUTDIR=/media/laur/8b3e5647-ec62-46b4-acf8-3edd161f78f5/rima
################# Processing of FASTQ files separately ##############################################
set -e
echo
echo Checking one FASTQ file
$VSEARCH --fastq_chars "$(ls -1 *.fastq | head -1)"

echo
echo Paired-end merging

for i in *_R1_*.fastq; do
	$VSEARCH --fastq_mergepairs "$i" --reverse "${i/_R1_/_R2_}" --fastqout "${i/R1_001.fastq/merged.fq}" --relabel "$i"
done

#Because relabel option does not work in this command, rename the merged fastqs to contain the sample names:

echo =================================================================================
echo
echo Relabeling merged fastq files
for f in *_merged.fq; do
	s=$(cut -d_ -f1-3 <<< "$f")
	echo Processing sample "$s"
	sed -i "s/^@M/@${s}|/" "$f"
done


echo
echo
echo Removing primers
for f in *merged.fq; do
	s=$(cut -d_ -f1 <<< "$f")
	echo ========================================================================
	echo Processing sample $s
	echo
	cutadapt -g ${PRIMER_F} -o "${s}_trimmed1.fq" "$f"
	cutadapt -a ${PRIMER_R_RC} -o "${s}_trimmed2.fq" "${s}_trimmed1.fq"
	cutadapt -g ${PRIMER_R} -o "${s}_trimmed3.fq" "${s}_trimmed2.fq"
	cutadapt -a ${PRIMER_F_RC} -o "${s}_trimmed4.fq" "${s}_trimmed3.fq"
done

##Quality filtering and sequence lengths, and dereplication in one loop
##note: sizeout option important for subsequent de novo chimera removal step
echo

for f in *_trimmed4.fq; do
	s=$(cut -d_ -f1 <<< "$f")
	echo keeping sequences above min length, with quality scores above 1
	echo processing sample $s
	$VSEARCH --threads $THREADS --fastq_filter "$f" --fastq_minlen 300 --fastq_maxee 1 --fastaout "${s}.qf.fa"
	echo
	echo Dereplication
	$VSEARCH --threads $THREADS --derep_fulllength "${s}.qf.fa" --sizeout --relabel Uniq --relabel "${s}." --output "${s}_uniques.fa"
done


echo Sum of unique sequences in all samples: $(cat *_uniques.fa | grep -c "^>")

# At this point there should be one uniques fasta file for each sample, quality filtered and dereplicated                              
echo
echo
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo Concatenating all samples and processing together
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##Note: fastq/a files from earlier being in the directory could get added, so they should not be in the working directory.
#Concatenate merged fastq's into one file:
cat *_uniques.fa > all.fa

echo
echo Dereplicate combined file and discard singletons.
$VSEARCH --derep_fulllength all.fa --minuniquesize 2 --sizein --sizeout --uc all.derep.uc --output all.derep.fa

echo
echo Total unique non-singleton sequences: $(grep -c "^>" all.derep.fa)

#Cluster before chimera detection, so uchime won't take as long:
echo
echo Cluster at 98% into OTUs, relabeling with OTU_
$VSEARCH --threads $THREADS --cluster_size all.derep.fa --id 0.98 --iddef 0 --sizein --sizeout --relabel OTU_ --centroids otusch.fa 

echo
echo OTUs before chimera detection: $(grep -c "^>" otusch.fa)


echo
echo De novo chimera detection
$VSEARCH --threads $THREADS --uchime_denovo otusch.fa --sizein --sizeout --nonchimeras otus.fa

echo
echo OTUs after chimera detection: $(grep -c "^>" otus.fa)


echo
echo Formatting raw sequences for OTU table. . . .
##Converting merged fastq's to fasta's:

for f in *_merged.fq; do
	s=$(cut -d. -f1 <<< "$f")
	$VSEARCH --threads $THREADS --fastq_filter ${s}.fq --fastaout ${s}.fa
done

## Substitute hyphens with underscores so VSEARCH doesn't chop the headers there:
for f in *_merged.fa; do
	sed -i "s/-/_/g" "$f"
done

## Concatenate merged fasta files into one:
cat *_merged.fa > allmerged.fa

echo
echo Creating OTU table
$VSEARCH --usearch_global allmerged.fa -db otus.fa --id 0.97 --otutabout otu_table.txt

echo Done.
date


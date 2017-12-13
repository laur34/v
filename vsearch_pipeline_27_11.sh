#! /bin/bash

## VSEARCH-based pipeline for Illumina paired-end metabarcoding reads
## 11.2017 LH, VB.

#All fastq files to use must be in working directory, and maybe unzipped, depending on vsearch version.

set -e

##Define Variables (Change to your primer sequences, these are MiniLGC):
PRIMER_F="TCATGGWACWGGWTGAACWGTWTAYCCYCC"
PRIMER_R_RC="TGRTTYTTTGGTCACCCTGAAGTTTACTCC"
PRIMER_R="GGAGTAAACTTCAGGGTGACCAAARAAYCA"
PRIMER_F_RC="GGRGGRTAWACWGTTCAWCCWGTWCCATGA"

THREADS=8
VSEARCH=$(which vsearch2)

#
echo
echo Checking FASTQ file

$VSEARCH --fastq_chars $(ls -1 *.fastq | head -1)

echo
echo Paired-end merging

for i in *_R1_*.fastq; do
	$VSEARCH --fastq_mergepairs ${i} --reverse ${i/_R1_/_R2_} --fastqout ${i/R1_*/merged.fq}
done

echo
echo Removing primers
for f in *merged.fq; do
	r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
	s=$(cut -d_ -f1 <<< "$f")
	echo ===========================
	echo Processing sample $s
	echo
	cutadapt -g ^${PRIMER_F} -o ${s}_trimmed1.fq $f
	cutadapt -a ${PRIMER_R_RC}$ -o ${s}_trimmed2.fq ${s}_trimmed1.fq
	cutadapt -g ^${PRIMER_R} -o ${s}_trimmed3.fq ${s}_trimmed2.fq
	cutadapt -a ${PRIMER_F_RC}$ -o ${s}_trimmed4.fq ${s}_trimmed3.fq
done

##Quality filtering and sequence lengths, and dereplication in one loop
##note: sizeout option important for subsequent de novo chimera removal step
echo

for f in *_trimmed4.fq; do
	s=$(cut -d_ -f1 <<< "$f")
	echo keeping sequences above min length, with quality scores above 1
	$VSEARCH --fastq_filter $f --fastq_minlen 300 --fastq_maxee 1 --fastaout ${s}.fa
	echo
	echo Dereplication
	$VSEARCH --derep_fulllength ${s}.fa --sizeout --relabel Uniq --relabel $s. --output ${s}_uniques.fa
done


echo Sum of unique sequences in all samples: $(cat *_uniques.fa | grep -c "^>")

# At this point there should be one uniques fasta file for each sample, quality filtered and dereplicated                              
echo
echo
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo Concatenating all samples and processing together
echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##remove previous combinded derep fasta files.
#Concatenate merged fastq's into one file:
cat *_uniques.fa > all.fa

echo
echo Dereplicate combined file and discard singletons.
$VSEARCH --derep_fulllength all.fa \
	--minuniquesize 2 --sizein --sizeout --uc all.derep.uc --output all.derep.fa

echo
echo Total unique non-singleton sequences: $(grep -c "^>" all.derep.fa)


#By not clustering (or pre-clustering) before chimera detection, uchime will take longer.
#If dataset is large, it might be a better idea to edit the script.
#$VSEARCH --threads $THREADS --cluster_size all.derep.fa \
#	--id 0.98 --sizein --sizeout --uc all.preclustered.uc --centroids all.preclustered.fa


echo
echo De novo chimera detection
$VSEARCH --threads $THREADS \
	 --uchime_denovo all.derep.fa --sizein --sizeout --nonchimeras all.derep.nochimeras.fa

echo
echo Unique sequences after de novo chimera detection: $(grep -c "^>" all.derep.nochimeras.fa)


echo
echo Cluster at 98% into OTUs, relabeling with OTU_
$VSEARCH --threads $THREADS \
	--cluster_size all.derep.nochimeras.fa --id 0.98 --iddef 0 --sizein --sizeout --relabel OTU_ --centroids otus.fa 


echo
echo Formatting raw sequences for OTU table

for i in $(ls *_merged.fq | cut -f 1-3 -d "_"); do sed -i "s/^@M/@${i}|/" ${i}_merged.fq; done
#Concatenate merged fastqs into one:
cat *merged.fq > allmerged.fq
#Convert allmerged fastq to fasta:
vsearch --fastq_filter allmerged.fq -fastaout allmerged.fa

#Edit it so that vsearch chops the headers only after the sample identifiers (will be the column names in table)
sed -i 's/|.*//' allmerged.fa
sed -i 's/_.*//' allmerged.fa
sed -i 's/-/_/g' allmerged.fa

echo
echo Creating OTU table
$VSEARCH --usearch_global allmerged.fa -db otus.fa --id 0.97 --otutabout otu_table.txt

echo Done.
date


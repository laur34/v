#!/bin/sh                                                                       

# This is based on Frederic Mahé's pipeline:
# "This is an example of a pipeline using vsearch to process data in the         
# Mothur 16S rRNA MiSeq SOP tutorial dataset to perform initial paired-end      
# read merging, quality filtering, chimera removal and OTU clustering."

THREADS=6
PERL=$(which perl)
#VSEARCH=$(which vsearch2)
VSEARCH=$(which vsearch)

#PRIMER_F="GGWACWGGWTGAACWGTWTAYCCYCC"
#PRIMER_R_RC="TGRTTYTTTGGTCACCCTGAAGTTTACTCC"
#PRIMER_R="TAAACTTCAGGGTGACCAAARAAYCA"
#PRIMER_F_RC="GGRGGRTAWACWGTTCAWCCWGTWCCATGA"

set -e

date

echo
echo "Checking FASTQ format version for one file"

$VSEARCH --threads $THREADS \
    --fastq_chars $(ls -1 *.fastq | head -1)

# Process samples.
# change some arguments in the merging and in the
# quality filtering if they are too stringent.
echo                                                        

for f in *_R1_*.fastq; do

    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
    s=$(cut -d_ -f1-3 <<< "$f")

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================

    $VSEARCH --threads $THREADS \
        --fastq_mergepairs $f \
        --reverse $r \
        --fastq_minovlen 50 \
        --fastq_maxdiffs 20 \
        --fastqout ${s}_merged.fq \
        --fastq_eeout

    # Commands to reverse-complement reads, combine with originals,
    # and then cutadapt will discard all that are in wrong direction.
    echo
    echo Reverse-complementing sample $s
    echo    
    $VSEARCH --threads $THREADS \
        --fastx_revcomp ${s}_merged.fq \
        --fastqout ${s}_mergedrc.fq

    echo
    echo Concatenating
    echo
#    cat $s.merged.fastq $s.merged.rc.fastq > $s.both.fq
done

# Commands to demultiplex and remove tags and primers added here:                       
echo
echo Removing primers
#echo sample $s
#cutadapt --discard-untrimmed -m 100 -g ${PRIMER_F} -o "${s}_trimmed1.fq" "${s}.both.fq"
#cutadapt --discard-untrimmed -m 100 -a ${PRIMER_R_RC} -o "${s}_trimmed2.fq" "${s}_trimmed1.fq"
while read filename PRIMER_F PRIMER_R_RC; do 
	s1=$(cut -d_ -f1-3 <<< "$filename")
	echo "========================================================================"
	echo Processing sample "$s1"
	echo
	cat $filename ${s1}_mergedrc.fq > "${s1}_both.fq"
	cutadapt --discard-untrimmed -g ${PRIMER_F} -o "${s1}_trimmed1.fq" "${s1}_both.fq" 2>&1 | tee -a pdf_trimming_primer_F.txt
	cutadapt --discard-untrimmed -a ${PRIMER_R_RC} -o "${s1}_trimmed2.fq" "${s1}_trimmed1.fq" 2>&1 | tee -a pdf_trimming_primer_R_RC.txt
done < samplesheet.txt

for f in *_R1_*.fastq; do

    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
    s=$(cut -d_ -f1 <<< "$f")
	s1=$(cut -d_ -f1-3 <<< "$f")
                                  
    echo
    echo Calculate quality statistics

    $VSEARCH --threads $THREADS \
        --fastq_eestats ${s1}_trimmed2.fq \
        --output $s1.stats
    echo
    echo Quality filtering

    $VSEARCH --threads $THREADS \
        --fastq_filter ${s1}_trimmed2.fq \
        --fastq_maxee 1 \
        --fastq_minlen 295 \
        --fastq_maxlen 355 \
        --fastq_maxns 0 \
        --fastaout $s1.filtered.fasta \
        --fasta_width 0

    echo
    echo Dereplicate at sample level and relabel with sample_n

    $VSEARCH --threads $THREADS \
        --derep_fulllength $s1.filtered.fasta \
        --strand plus \
        --output $s1.derep.fasta \
        --sizeout \
        --uc $s1.derep.uc \
        --relabel $s1. \
        --fasta_width 0

done

echo Sum of unique sequences in each sample: $(cat *.derep.fasta | grep -c "^>")

# At this point there should be one fasta file for each sample                  
# It should be quality filtered and dereplicated.                               


echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Merge all samples

rm -f all.derep.fasta all.nonchimeras.derep.fasta
cat *.derep.fasta > all.fasta

echo
echo Dereplicate across samples and remove singletons

$VSEARCH --threads $THREADS \
    --derep_fulllength all.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.derep.uc \
    --output all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" all.derep.fasta)

echo
echo Precluster at 98% before chimera detection

$VSEARCH --threads $THREADS \
    --cluster_size all.derep.fasta \
    --id 0.98 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.preclustered.uc \
    --centroids all.preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>" all.preclustered.fasta)

echo
echo De novo chimera detection

$VSEARCH --threads $THREADS \
    --uchime_denovo all.preclustered.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras all.denovo.nonchimeras.fasta \

echo Unique sequences after de novo chimera detection: $(grep -c "^>" all.denovo.nonchimeras.fasta)

echo
#echo Reference chimera detection

#$VSEARCH --threads $THREADS \
#    --uchime_ref all.denovo.nonchimeras.fasta \
#    --db $REF \
#    --sizein \
#    --sizeout \
#    --fasta_width 0 \
#    --nonchimeras all.ref.nonchimeras.fasta

#echo Unique sequences after reference-based chimera detection: $(grep -c "^>" all.ref.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

#$PERL ../map.pl all.derep.fasta all.preclustered.uc all.ref.nonchimeras.fasta > all.nonchimeras.derep.fasta
$PERL ./map.pl all.derep.fasta all.preclustered.uc all.denovo.nonchimeras.fasta > all.nonchimeras.derep.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample

#$PERL ../map.pl all.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta
$PERL ./map.pl all.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta

echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" all.nonchimeras.fasta)

echo
echo Cluster at 97% and relabel with OTU_n, generate OTU table
#First substitute hyphens with underscores:
sed -i 's/-/_/g' all.nonchimeras.fasta

$VSEARCH --threads $THREADS \
    --cluster_size all.nonchimeras.fasta \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel OTU_ \
    --centroids all.otus.fasta \
    --otutabout all.otutab.txt

echo
echo Number of OTUs: $(grep -c "^>" all.otus.fasta)

echo
echo Done

date


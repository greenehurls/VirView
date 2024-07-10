#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh


conda deactivate
conda activate nanite
mkdir $1sequenceresults

for file in $1*
do
  
  gunzip -c $file \
        | NanoFilt -q 8 \
        | gzip \
>> $1sequenceresults/filtered.fastq.gz

done

conda deactivate
conda activate bgminimap2
 minimap2 \
   -ax map-ont \
   --secondary=no \
   --sam-hit-only \
   $2 \
   $1sequenceresults/filtered.fastq.gz \
 > $1sequenceresults/alignment.sam
conda deactivate
conda activate samtools
samtools sort -o $1sequenceresults/alignment.sam $1/sequenceresults/alignment.sam
samtools depth $1sequenceresults/alignment.sam > $1/sequenceresults/coverage.tsv

samtools view $1sequenceresults/alignment.sam | awk '{print length($10)}' | head -1000 | sort -u > $1sequenceresults/readlengths.tsv

samtools view -F 260 $1sequenceresults/alignment.sam | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > $1sequenceresults/counts.txt
gzip -dk $1sequenceresults/filtered.fastq.gz

echo $(cat $1sequenceresults/filtered.fastq|wc -l)/4|bc >> $1/sequenceresults/counts.txt
(awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $2 | tail -n +2) > $1sequenceresults/genomelength.tsv 

Rscript CoverageGraph2.R $1sequenceresults/coverage.tsv $1sequenceresults/genomelength.tsv $1sequenceresults/readlengths.tsv $1sequenceresults/report.pdf $1sequenceresults/coverage.html $3 $4



#!/bin/bash
echo "usage: sh star_align_genome.sh genome_idx_dir genome_fa genome_gff read_fn outputprefix adapter sjdbGTFtagExonParentTranscript"
genome_idx=$1
genome_fa=$2
genome_gff=$3
read_fn=$4
outputprefix=$5
adapter=$6
pname=$7
nproc=15
nmismatch=1
align_params="--clip3pAdapterSeq ${adapter} --seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax ${nmismatch} --outFilterIntronMotifs RemoveNoncanonical"
SAM_params="--outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM"
if  [ ! -d ${genome_idx} ];  then
    echo "building genome index..."
    mkdir -p ${genome_idx}
    STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${genome_idx} --genomeFastaFiles ${genome_fa} --sjdbGTFfile ${genome_gff} --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentTranscript $pname --sjdbOverhang 1
fi
STAR --runThreadN $nproc --genomeDir ${genome_idx} --readFilesIn ${read_fn} --outFileNamePrefix ${outputprefix} ${SAM_params} ${align_params} --quantMode TranscriptomeSAM


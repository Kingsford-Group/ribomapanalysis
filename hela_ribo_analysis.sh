#!/bin/bash
work_dir=/home/hw1/scratch/scratch2/ribomap-playground/hela/
ribomap_dir=/home/hw1/scratch/scratch2/software_testing/ribomap/
fasta_dir=${work_dir}data/fasta/
ref_dir=${work_dir}ref/
#=============================
# step 1: download sra & refs
#=============================
echo "downloading transcriptome reference data..."
gtf_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.annotation.gtf.gz
tfa_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.pc_transcripts.fa.gz
pfa_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_18/gencode.v18.pc_translations.fa.gz
# wget -P ${ref_dir} -N ${gtf_url}
# wget -P ${ref_dir} -N ${tfa_url}
# wget -P ${ref_dir} -N ${pfa_url}
# echo "unzipping data..."
# gunzip -f ${ref_dir}*.gz
# echo "downloading Hela cell reads..."
riboseq_url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546920/suppl/GSM546920_filtered_sequence.txt.gz
rnaseq_url=ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM546nnn/GSM546921/suppl/GSM546921_filtered_sequence.txt.gz
# wget -P ${fasta_dir} -N ${rnaseq_url}
# wget -P ${fasta_dir} -N ${riboseq_url}
#=============================
# step 2: process reference
# build cds range file
#=============================
echo "filtering transcripts...."
gtf=${ref_dir}${gtf_url##*/}
gtf=${gtf%.gz}
tfa=${ref_dir}${tfa_url##*/}
tfa=${tfa%.gz}
pfa=${ref_dir}${pfa_url##*/}
pfa=${pfa%.gz}
# python filter_gencode_transcript.py ${gtf} ${tfa} ${pfa}
#====================================
# step 3: build contaminant sequence
#====================================
# echo "downloading ncRNA from Ensembl..."
# wget -P ${ref_dir} -N ftp://ftp.ensembl.org/pub/release-78/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
# echo "downloading tRNA from gtrnadb..."
# wget -P ${ref_dir} -N http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz
# gunzip -f ${ref_dir}*.gz
# echo "merging rRNA and tRNA..."
# python build_contaminant.py ${ref_dir}
echo "done preparing data for ribomap"
#=============================
# step 4: run ribomap
#=============================
rnaseq_fq=${fasta_dir}${rnaseq_url##*/}
riboseq_fq=${fasta_dir}${riboseq_url##*/}
transcript_fa=${tfa%.*}_filter.fa
contaminant_fa=${ref_dir}human_contaminant.fa
cds_range=${tfa%.*}_cds.txt
star_idx_dir=${work_dir}/StarIndex/gencodev18/
adapter=TCGTATGCCGTCTTCTGCTTG
min_fplen=25
max_fplen=36
offset=offset.txt
${ribomap_dir}/scripts/run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --work_dir ${work_dir} --star_idx_dir ${star_idx_dir} --offset ${offset} --adapter ${adapter} --min_fplen ${min_fplen} --max_fplen ${max_fplen} --force true

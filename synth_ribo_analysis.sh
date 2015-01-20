#!/bin/bash
work_dir=/home/hw1/scratch/scratch2/ribomap-playground/hela/
ribomap_dir=/home/hw1/scratch/scratch2/software_testing/ribomap/
rnaseq_fq=${work_dir}synthetic_data/synth_rnaseq.fq.gz
riboseq_fq_prefix=${work_dir}synthetic_data/synth_riboseq_
transcript_fa=${work_dir}ref/gencode.v18.pc_transcripts_filter.fa
cds_range=${work_dir}ref/gencode.v18.pc_transcripts_cds.txt
star_idx_dir=${work_dir}StarIndex/gencodev18/
min_fplen=25
max_fplen=30
offset=12
for e in 005 01 02; do 
    ${ribomap_dir}scripts/run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq_prefix}$e.fq.gz --transcript_fa ${transcript_fa} --cds_range ${cds_range} --work_dir ${work_dir} --star_idx_dir ${star_idx_dir} --offset ${offset} --min_fplen ${min_fplen} --max_fplen ${max_fplen} --useSecondary true --output_dir ${work_dir}ribomap/ 
    ${ribomap_dir}scripts/run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq_prefix}$e.fq.gz --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --work_dir ${work_dir} --star_idx_dir ${star_idx_dir} --offset ${offset} --adapter ${adapter} --min_fplen ${min_fplen} --max_fplen ${max_fplen} --useSecondary false --output_dir ${work_dir}star_prime/
done

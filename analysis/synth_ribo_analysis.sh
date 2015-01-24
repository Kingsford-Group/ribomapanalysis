#!/bin/bash
hela_dir=/home/hw1/scratch/scratch2/ribomap-playground/hela/
work_dir=${hela_dir}synthetic_data/
ribomap_dir=/home/hw1/scratch/scratch2/software_testing/ribomap/
rnaseq_fq=${work_dir}synth_rnaseq.fq.gz
riboseq_fq_prefix=${work_dir}synth_riboseq_
transcript_fa=${hela_dir}ref/gencode.v18.pc_transcripts_filter.fa
cds_range=${hela_dir}ref/gencode.v18.pc_transcripts_cds.txt
star_idx_dir=${hela_dir}StarIndex/
sm_odir=${work_dir}sm_align_quant
min_fplen=25
max_fplen=30
offset=12
adapter=N # pseudo adapter to avoid star from crashing
ribo_cmd="${ribomap_dir}scripts/run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --transcript_fa ${transcript_fa} --cds_range ${cds_range} --softClipping false --adapter ${adapter} --min_fplen ${min_fplen} --max_fplen ${max_fplen} --offset ${offset} --work_dir ${work_dir} --star_idx_dir ${star_idx_dir} --sailfish_dir ${sm_odir}"

for e in 005 01 02; do 
    ${ribo_cmd} --riboseq_fq ${riboseq_fq_prefix}$e.fq.gz --output_dir ${work_dir}ribomap
    ${ribo_cmd} --riboseq_fq ${riboseq_fq_prefix}$e.fq.gz --output_dir ${work_dir}star_prime --tabd_cutoff -1 --useSecondary false 
done

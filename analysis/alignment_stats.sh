#!/bin/bash
align_dir=/home/hw1/scratch/scratch2/ribomap-playground/hela/alignment/
bam_fn=${align_dir}GSM546920_filtered_sequence_transcript_Aligned.out.bam
./bam_info ${bam_fn} ${align_dir}

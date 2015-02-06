hela_align_dir=~/scratch/scratch2/ribomap-playground/hela/alignment/
hela_genome_bam=${hela_align_dir}GSM546920_filtered_sequence_genome_Aligned.out.bam
hela_transcript_bam=${hela_align_dir}GSM546920_filtered_sequence_transcript_Aligned.out.bam
#./multimapper_info  ${hela_genome_bam} ${hela_transcript_bam}
#./bam_info ${hela_transcript_bam} ${hela_align_dir}

mouse_align_dir=~/scratch/scratch2/ribomap-playground/mouse/alignment/
mouse_genome_bam=${mouse_align_dir}ribo_mesc_yeslif_genome_Aligned.out.bam
mouse_transcript_bam=${mouse_align_dir}ribo_mesc_yeslif_transcript_Aligned.out.bam
#./multimapper_info ${mouse_genome_bam} ${mouse_transcript_bam}
#./bam_info ${mouse_transcript_bam} ${mouse_align_dir}

yeast_align_dir=~/scratch/scratch2/ribomap-playground/yeast/alignment/
yeast_genome_bam=${yeast_align_dir}BY_FP_genome_Aligned.out.bam
yeast_transcript_bam=${yeast_align_dir}BY_FP_transcript_Aligned.out.bam
./multimapper_info ${yeast_genome_bam} ${yeast_transcript_bam}
./bam_info ${yeast_transcript_bam} ${yeast_align_dir}



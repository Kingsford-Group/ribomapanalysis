#!/bin/bash
get_data=false
work_dir=/home/hw1/scratch/scratch2/ribomap-playground/mouse/
ref_dir=${work_dir}ref/
genome_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/GRCm38.p3.genome.fa.gz
annotation_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.annotation.gtf.gz
genome_in=${ref_dir}${genome_url##*/}
genome_in=${genome_in%.gz}
genome_out=${ref_dir}genome.fa
if [ "${get_data}" = true ]; then
    wget -P ${ref_dir} -N ${genome_url}
    wget -P ${ref_dir} -N ${annotation_url}
    gunzip -f ${ref_dir}*.gz
    python filter_gencode_genome.py ${genome_in} ${genome_out}
fi
python filter_gencode_genome.py ${genome_in} ${genome_out}
star_idx=${work_dir}StarIndex/genome/
gtf_fn=${ref_dir}${annotation_url##*/}
gtf_fn=${gtf_fn%.gz}
read_fn=${work_dir}alignment/ribo_mesc_yeslif_rrna_Unmapped.out.mate1
outputprefix=${work_dir}alignment/ribo_mesc_yeslif_genome_
adapter=CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTGAA
pname=transcript_id
./star_align_genome.sh ${star_idx} ${genome_out} ${gtf_fn} ${read_fn} $outputprefix $adapter $pname

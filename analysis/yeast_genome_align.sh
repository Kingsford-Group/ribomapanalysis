#!/bin/bash
get_data=false
work_dir=/home/hw1/scratch/scratch2/ribomap-playground/yeast/
ref_dir=${work_dir}ref/
genome_url=http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
gff_in=${ref_dir}${genome_url##*/}
gff_in=${gff_in%%.*}/saccharomyces_cerevisiae_R64-1-1_20110208.gff
gff_out=${ref_dir}genome_annotations.gff
fa_out=${ref_dir}genome.fa
if [ "${get_data}" = true ]; then
    wget -P ${ref_dir} -N ${genome_url}
    tar zxvf ${gff_in} -C ${ref_dir}
    head_line_cnt=$(grep -n "##FASTA" ${gff_in})
    head_line_cnt=${head_line_cnt%:*}
    echo "generating gff..."
    head -n ${head_line_cnt} ${gff_in} > ${gff_out}
    tot_line_cnt=$(wc -l ${gff_in})
    tot_line_cnt=${tot_line_cnt% *}
    tail_line_cnt=$(expr ${tot_line_cnt} - ${head_line_cnt})
    echo "generating reference fasta..."
    tail -n ${tail_line_cnt} ${gff_in} > ${fa_out}
fi
star_idx=${work_dir}StarIndex/genome/
read_fn=${work_dir}alignment/BY_FP_rrna_Unmapped.out.mate1
outputprefix=${work_dir}alignment/BY_FP_genome_
adapter=CTGTAGGCACCATCAAT
pname=Name
./star_align_genome.sh ${star_idx} ${fa_out} ${gff_out} ${read_fn} $outputprefix $adapter $pname

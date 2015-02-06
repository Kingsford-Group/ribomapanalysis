#!/bin/bash
work_dir=/home/hw1/scratch/scratch2/ribomap-playground/yeast/
ribomap_dir=/home/hw1/scratch/scratch2/software_testing/ribomap/
data_dir=${work_dir}data/
sra_dir=${data_dir}sra/
fasta_dir=${data_dir}fasta/
ref_dir=${work_dir}ref/
get_data=true
if [ "${get_data}" = true ]; then
    #=============================
    # make directories
    #=============================
    mkdir -p ${sra_dir}
    mkdir -p ${fasta_dir}
    mkdir -p ${ref_dir}
    #=============================
    # step 1: download sra
    #=============================
    echo "downloading data..."
    #wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476344/SRR1177156/SRR1177156.sra
    #wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476345/SRR1177157/SRR1177157.sra
    #==============================
    # step 2: convert sra to fastq
    #==============================
    echo "converting sra to fastq.gz..."
    for file in ${sra_dir}*.sra; do
	fastq-dump $file -O ${fasta_dir} --gzip
    done
    #==============================
    # step 3: rename files 
    #==============================
    echo "renaming fastq files to be more informative..."
    mv ${fasta_dir}SRR1177156.fastq.gz ${fasta_dir}BY_mRNA.fastq.gz
    mv ${fasta_dir}SRR1177157.fastq.gz ${fasta_dir}BY_FP.fastq.gz
    exit
    #==============================
    # step 4: build references
    #==============================
    echo "downloading references..."
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_1000.fasta.gz
    wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz
    nc_url=ftp://ftp.ensembl.org/pub/release-78/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz
    trna_url=http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz
    wget -P ${ref_dir} -N ${nc_url}
    wget -P ${ref_dir} -N ${trna_url}
    gunzip -f ${ref_dir}*.gz
    python filter_yeast_transcript.py ${ref_dir} 33
    echo "merging rrna and trna..."
    rrna_fa=${ref_dir}${nc_url##*/}
    rrna_fa=${rrna_fa%.gz}
    trna_fa=${ref_dir}${trna_url##*/}
    trna_fa=${trna_fa%.gz}
    python build_contaminant.py ${rrna_fa} ${trna_fa} Saccharomyces_cerevisiae ${ref_dir}yeast_contaminant.fa
fi
#==============================
# step 4: run ribomap
#==============================
rnaseq_fq=${fasta_dir}BY_mRNA.fastq.gz
riboseq_fq=${fasta_dir}BY_FP.fastq.gz
transcript_fa=${ref_dir}protein_coding_33_filtered.fasta
contaminant_fa=${ref_dir}yeast_contaminant.fa
cds_range=${ref_dir}cds_range.txt
offset=offset.txt
ribo_cmd="${ribomap_dir}scripts/run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --offset ${offset} --work_dir ${work_dir}"
# ribomap
${ribo_cmd} --output_dir ${work_dir}ribomap
# star prime
${ribo_cmd} --output_dir ${work_dir}star_prime --tabd_cutoff -1 --useSecondary false

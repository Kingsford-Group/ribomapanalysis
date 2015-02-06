#!/bin/bash
work_dir=/home/hw1/scratch/scratch2/ribomap-playground/mouse/
ribomap_dir=/home/hw1/scratch/scratch2/software_testing/ribomap/
data_dir=${work_dir}data/
sra_dir=${data_dir}sra/
fasta_dir=${data_dir}fasta/
ref_dir=${work_dir}ref/
get_data=false
#=============================
# urls
#=============================
gtf_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.annotation.gtf.gz
tfa_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.pc_transcripts.fa.gz
pfa_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.pc_translations.fa.gz
gtf=${ref_dir}${gtf_url##*/}
gtf=${gtf%.gz}
tfa=${ref_dir}${tfa_url##*/}
tfa=${tfa%.gz}
pfa=${ref_dir}${pfa_url##*/}
pfa=${pfa%.gz}
if [ "${get_data}" = true ]; then
    #=============================
    # make directories
    #=============================
    mkdir -p ${sra_dir}
    mkdir -p ${fasta_dir}
    mkdir -p ${ref_dir}
    #=============================
    # step 1: download sra & refs
    #=============================
    echo "downloading sra..."
    # ES cell feeder-free, w/ LIF 60 s CYH (100 ug/ml) mrna_mesc_yeslif Illumina GAII
    wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084812/SRR315595/SRR315595.sra
    wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084812/SRR315596/SRR315596.sra
    # # ES cell feeder-free, w/ LIF 60 s CYH (100 ug/ml) ribo_mesc_chx Illumina GAII
    # wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084815/SRR315601/SRR315601.sra
    # wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084815/SRR315602/SRR315602.sra
    # ES cell feeder-free, w/ LIF 60 s CYH (100 ug/ml) ribo_mesc_yeslif Illumina GAII
    wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084824/SRR315624/SRR315624.sra
    wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084824/SRR315625/SRR315625.sra
    wget -P ${sra_dir} -N ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX084%2FSRX084824/SRR315626/SRR315626.sra
    echo "downloading references..."
    wget -P ${ref_dir} -N ${gtf_url}
    wget -P ${ref_dir} -N ${tfa_url}
    wget -P ${ref_dir} -N ${pfa_url}
    echo "unzipping data..."
    gunzip -f ${ref_dir}*.gz
    #==============================
    # step 2: convert sra to fastq
    #==============================
    echo "converting sra to fastq..."
    for file in ${sra_dir}*.sra; do
	fastq-dump $file -O ${fasta_dir}
    done
    #==============================
    # step 3: zip fastq
    #==============================
    echo "zipping fastq files..."
    gzip -f -c ${fasta_dir}SRR315595.fastq ${fasta_dir}SRR315596.fastq > ${fasta_dir}mrna_mesc_yeslif.fq.gz
    gzip -f -c ${fasta_dir}SRR315624.fastq ${fasta_dir}SRR315625.fastq ${fasta_dir}SRR315626.fastq > ${fasta_dir}ribo_mesc_yeslif.fq.gz
    #=============================
    # step 4: process reference
    # build cds range file
    #=============================
    echo "filtering transcripts...."
    python filter_gencode_transcript.py ${gtf} ${tfa} ${pfa}
    #====================================
    # step 5: build contaminant sequence
    #====================================
    nc_url=ftp://ftp.ensembl.org/pub/release-78/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz
    trna_url=http://gtrnadb.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz
    echo "downloading ncRNA from Ensembl..."
    wget -P ${ref_dir} -N ${nc_url}
    echo "downloading tRNA from gtrnadb..."
    wget -P ${ref_dir} -N ${trna_url}
    gunzip -f ${ref_dir}*.gz
    echo "merging rRNA and tRNA..."
    rrna_fa=${ref_dir}${nc_url##*/}
    rrna_fa=${rrna_fa%.gz}
    trna_fa=${ref_dir}${trna_url##*/}
    trna_fa=${trna_fa%.gz}
    python build_contaminant.py ${rrna_fa} ${trna_fa} Mus_musculus ${ref_dir}mouse_contaminant.fa
    echo "done preparing data for ribomap"
fi
#=============================
# step 6: run ribomap
#=============================
rnaseq_fq=${fasta_dir}mrna_mesc_yeslif.fq.gz
riboseq_fq=${fasta_dir}ribo_mesc_yeslif.fq.gz
transcript_fa=${tfa%.*}_filter.fa
contaminant_fa=${ref_dir}mouse_contaminant.fa
cds_range=${tfa%.*}_cds.txt
adapter=CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTGAA
min_fplen=25
max_fplen=35
offset=offset.txt
ribo_cmd="${ribomap_dir}scripts/run_ribomap.sh --rnaseq_fq ${rnaseq_fq} --riboseq_fq ${riboseq_fq} --transcript_fa ${transcript_fa} --contaminant_fa ${contaminant_fa} --cds_range ${cds_range} --adapter ${adapter} --min_fplen ${min_fplen} --max_fplen ${max_fplen} --offset ${offset} --work_dir ${work_dir}"
# ribomap
${ribo_cmd} --output_dir ${work_dir}ribomap
# star prime
${ribo_cmd} --output_dir ${work_dir}star_prime --tabd_cutoff -1 --useSecondary false

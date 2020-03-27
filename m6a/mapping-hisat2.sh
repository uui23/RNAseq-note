#!/bin/sh
#PBS -q middleq
#PBS -l mem=10gb,walltime=500:00:00,nodes=1:ppn=5
#PBS -o ./script_stout/gzip_stdout.txt
#PBS -e ./script_stout/gzip_stderr.txt
#HSCHED -s m6A_miCLIP+gzip+humo

cd /share_bio/unisvx1/yangyg_group/yangxin/ZhouJY/xiada/CLIP_SEQ_data/clip_seq_data/ctrl/1.rawdata/CTRL_TKR180900812


cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o CTRL_TKR_1-c.fq.gz -p CTRL_TKR_2-c.fq.gz CTRL_TKR180900812_1.fq.gz CTRL_TKR180900812_2.fq.gz

java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 CTRL_TKR_1-c.fq.gz CTRL_TKR_2-c.fq.gz CTRL_TKR_1-c1.fq.gz ./unpaired_1.fastq.gz CTRL_TKR_2-c1.fq.gz ./unpaired_2.fastq.gz LEADING:10 TRAILING:10 SLIDINGWINDOW:4:10 MINLEN:20

/software/biosoft/software/hisat2/hisat2-2.0.1-beta/hisat2 -p 5 -N 1 --dta -x /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Human/ENSEMBLE_72/index/genome.fa -1 CTRL_TKR_1-c1.fq.gz -2 CTRL_TKR_2-c1.fq.gz -S CTRL_TKR.sam


cat CTRL_TKR.sam|grep -E "^@|NH:i:1$|NH:i:1[^0-9]"|samtools view -S - -b -o CTRL_TKR.bam

bedtools bamtobed -i CTRL_TKR.bam -split  >CTRL_TKR.bed

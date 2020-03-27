#!/bin/sh
#PBS -q tmp
#PBS -l mem=40gb,walltime=500:00:00,nodes=1:ppn=8
#PBS -o ./log/ko4-3.log
#PBS -e ./log/ko4-3.err
#HSCHED -s empyro+fastqc+zebrafish

################################################################################################


###/pnas/public_data/Reference/temp/Reference/Mouse/ucsc/mm10/Sequence/WholeGenomeFasta/mm10.fa

#/pnas/yangyg_group/yangxin/reference/mouse/ENSEMBLE_68/index/genome.fa
#/pnas/yangyg_group/yangxin/reference/mouse/ENSEMBLE_68/GTF/Mus_musculus.GRCm38.68.gtf


#/pnas/yangyg_group/sunbf/PuJM/IOZ-LiuFeng-zyf/1.rawdata/cKO-4_FRAS190322633-1a/cKO-4_FRAS190322633-1a_1.fq.gz

ref_genome=/pnas/yangyg_group/yangxin//reference/mouse/ENSEMBLE_68/index/genome.fa
#ref_45sRNA=/pnas/yangyg_group/zhangting/Reference/human/ENSEMBLE_72/45S/45S_Homo.fasta
refGene=/pnas/yangyg_group/yangxin/reference/mouse/ENSEMBLE_68/GTF/Mus_musculus.GRCm38.68.chr.gtf

R1=/pnas/yangyg_group/sunbf/PuJM/IOZ-LiuFeng-zyf/1.rawdata/cKO-4_*/*_1.fq.gz
R2=/pnas/yangyg_group/sunbf/PuJM/IOZ-LiuFeng-zyf/1.rawdata/cKO-4_*/*_2.fq.gz
 

cd /pnas/yangyg_group/pujm/mHPSC_RNAseq/ko4

<<!
mkdir clean_data
mkdir ./clean_data/fastqc
mkdir log
cd clean_data

#fastqc $R1  -t 8 -o ./fastqc
#fastqc $R2  -t 8 -o ./fastqc

cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o trim1_1.fq.gz -p trim1_2.fq.gz $R1 $R2

fastqc trim1_1.fq.gz  -t 8 -o ./fastqc
fastqc trim1_2.fq.gz  -t 8 -o ./fastqc

# ###trim卡质量
java -Xmx4g -jar /pnas/yangyg_group/yangxin/software/new/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 ./trim1_1.fq.gz ./trim1_2.fq.gz ./trim2_1.fq.gz ./unpaired_1.fq.gz ./trim2_2.fq.gz ./unpaired_2.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35

fastqc trim2_1.fq.gz  -t 8 -o ./fastqc
fastqc trim2_2.fq.gz  -t 8 -o ./fastqc

rm trim1_1.fq.gz trim1_2.fq.gz

#cd ../
#mkdir tophat_45s
#cd tophat_45s
#hisat2 -p 8 --dta -x $ref_45sRNA -1 ../clean_data/trim2_1.fq.gz -2 ../clean_data/trim2_2.fq.gz -S ./hisat2.sam


cd ../
mkdir hisat2
cd hisat2

hisat2 -p 8 --dta -x $ref_genome -1 ../clean_data/trim2_1.fq.gz -2 ../clean_data/trim2_2.fq.gz -S ./hisat2.sam
samtools view hisat2.sam -h -b > accepted_hits.bam
rm hisat2.sam
samtools view accepted_hits.bam -h -q 20 -b | samtools sort -@ 8 -m 3G -l 9 -o sort.bam -
samtools index sort.bam
bedtools bamtobed -i sort.bam -split > uniqmap.bed

samtools flagstat sort.bam > readsCount

echo "" >> readsCount
echo "" >> readsCount
echo "" >> readsCount

intersectBed -a uniqmap.bed -b /pnas/yangyg_group/yangxin/reference/mouse/ENSEMBLE_68/GeneList/single-transcript-chr -wa -wb -f 0.5 > tmp.unique
echo 'unique mapping reads:' >> readsCount
awk '!x[$4]++' uniqmap.bed|wc >> readsCount
echo 'unique nc+mRNA reads:' >> readsCount
awk '!x[$4]++' tmp.unique|wc >> readsCount
echo 'unique ncRNA reads:' >> readsCount
awk '$13!="mRNA"' tmp.unique|awk '!x[$4]++'|wc >> readsCount
echo 'unique mRNA region:' >> readsCount
cat tmp.unique|awk '$13=="mRNA"'|awk '!x[$4]++'|cut -f14|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount
rm tmp.unique
!

#HT_seq
cd /pnas/yangyg_group/pujm/mHPSC_RNAseq/ko4/hisat2
mkdir ./htSeq
/software/biosoft/software/python/python2.7_2018_12/bin/htseq-count -m union -f bam -s no sort.bam $refGene > htSeq/no_union.out
/software/biosoft/software/python/python2.7_2018_12/bin/htseq-count -m union -f bam -s yes sort.bam $refGene > htSeq/yes_union.out
/software/biosoft/software/python/python2.7_2018_12/bin/htseq-count -m union -f bam -s reverse sort.bam $refGene > htSeq/rev_union.out

<<!
mkdir ./bedgraph
genomeCoverageBed -bga -ibam ./sort.bam -g /pnas/yangyg_group/yangxin/reference/mouse/ENSEMBLE_68/GeneList/mm-size -split|sort - -k1,1 -k2,2n > ./bedgraph/sort.bedgraph
bedGraphToBigWig ./bedgraph/sort.bedgraph /pnas/yangyg_group/yangxin/reference/mouse/ENSEMBLE_68/GeneList/mm-size ./bedgraph/norm.bw

awk 'NR==FNR{sum+=($3-$2)*$4;len+=$3-$2}NR>FNR{print $1,$2,$3,$4*len/sum}' ./bedgraph/sort.bedgraph ./bedgraph/sort.bedgraph > ./bedgraph/norm.bedgraph
bedGraphToBigWig ./bedgraph/norm.bedgraph /pnas/yangyg_group/yangxin/reference/mouse/ENSEMBLE_68/GeneList/mm-size ./bedgraph/norm.bw
!



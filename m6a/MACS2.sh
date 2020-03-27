#!/bin/sh
#PBS -q middleq
#PBS -l mem=10gb,walltime=500:00:00,nodes=1:ppn=5
#PBS -o ./script_stout/gzip_stdout.txt
#PBS -e ./script_stout/gzip_stderr.txt
#HSCHED -s m6A_miCLIP+gzip+humo


##peak calling
##主要调整的参数：
##1、--gsize=1.87e9	根据物种参考基因组大小来调整（也有参考转录组的）
##2、--keep-dep		是否保留重复序列（reads位置完全一样，就保留一条。这个酌情考虑，peak数量上会有很大差别）
##3、-q 0.05/-f 0.001	fdr小于0.05／p小于0.0001（同样酌情考虑，得根据call到的peak数目自己衡量。建议多用fdr）

cd /share_bio/unisvx1/yangyg_group/yangxin/ZhouJY/ck/m6A_OV_0

/software/biosoft/software/python2.7.2/bin/macs2 callpeak -t uniqmap-m6a.bam -c uniqmap-input.bam  -f BAM --nomodel -g hs --keep-dup all -n m6a -q 0.05


/software/biosoft/software/python2.7.2/bin/macs2 callpeak -t uniqmap-m6a.bed -c uniqmap-input.bed -f BED --nomodel -g hs -n m6a1 -q 0.05

cat m6a_peaks.bed|awk -v OFS="\t" '{print $1,$2,$3,$4,$5}'|intersectBed -a - -b /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Mouse/ENSEMBLE_68/GeneList/single-transcript-chr -wa -wb -f 0.5|awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$9}'|sort|uniq |cut -f1,2,3,6 >ok.bed

perl distribution_plot_human.pl

shuffleBed -i ok.bed -incl /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Human/ENSEMBLE_72/GeneList/single-transcript-chr -g /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Human/ENSEMBLE_72/GeneList/chr.size  > shuffle.bed

/share_bio/unisvx1/yangyg_group/yangxin/software/homer/bin/findMotifsGenome.pl ok.bed /share_bio/unisvx1/yangyg_group/yangxin/REFERENCE/Human/ENSEMBLE_72/index/genome.fa ./TOTAL-5nt/ -size given -p 1 -len 5,6,7,8 -rna -chopify -norevopp -cache 1000 -bg shuffle.bed
























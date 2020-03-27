# step2 pickDEG.R

#input RPKM.selected edger deseq2

# input	deseq2/edger
# 		rpkm >1 expSet
# output 	diffgene_cutoff/up/down(list&data.frame_with rpkm) (enseble)(UP&DOWM) 
# 		diffgene_ranked

# save(expset.mean,
#	 RPKM.selected,
#	 file='res1_getDEG.R')

# # Usage
# pickDEG(deseq2)
# pickDEG(edger)
# source(step2_pickDEG.R)

rm(list = ls()) 
options(stringsAsFactors = F)

load("res1_getDEG.R")

deseq2 <- read.delim("./src/deseq2.txt")  #deseq2 get from server

FC_cutoff = 1.5

expSet.selected = RPKM.selected

# pickDEG <- function(FC_cutoff,pro='test'){

# 	if(pro =="deseq2")
# 	return(pickDEG_deseq2(deseq2,RPKM.selected,FC_cutoff))

# 	# if(pro =="edger")
# 	# return(pickDEG_edger(edger,RPKM.selected,FC_cutoff))

# 	# save(deseq2,
# 	# deseq2.ranked,
# 	# deseq2.cutoff,
# 	# edger,
# 	# edger.ranked,
# 	# edger.cutoff,file="res2_pickDEG.R")
# }



pickDEG_deseq2 <- function(deseq2,
							expSet.selected,
							FC_cutoff){

	deseq2$change <- as.factor(
		ifelse(deseq2$padj < 0.05 & abs(deseq2$log2FoldChange) > log2(FC_cutoff), 
			ifelse(deseq2$log2FoldChange > log2(FC_cutoff),'UP','DOWN'),'NOT'))

	deseq2.rpkm.seleced<-deseq2[rownames(expSet.selected),];                     # > 10k
	deseq2.cutoff <- subset(deseq2.rpkm.seleced, deseq2.rpkm.seleced$change %in% c("UP","DOWN")); #selected by rpkm #2k # for valcona plot
	deseq2.ranked <- as.data.frame(row.name = rownames(deseq2), deseq2[, c("log2FoldChange")]) #for GSEA
	save(deseq2,deseq2.ranked,deseq2.cutoff,file="res2_pickDEG_deseq2.R")
}








pickDEG_edger <- function(edger,
							expSet.selected,
							FC_cutoff){

	edger$change = as.factor(ifelse(edger$PValue < 0.05 & abs(edger$logFC) > log2(FC_cutoff),
								ifelse(edger$logFC > log2(FC_cutoff) ,'UP','DOWN'),'NOT'))

	edger.ranked <- as.data.frame(row.name = rownames(edger), edger[, c("logFC")]) #for GSEA
	edger.rpkm.seleced<-edger[rownames(RPKM.selected),]                     # > 10k
	edger.cutoff <- subset(edger.rpkm.seleced, edger.rpkm.seleced$change %in% c('UP', 'DOWN')) #selected by rpkm #2k # for valcona plot
	save(edger,
		edger.ranked,
		edger.cutoff,file="step2_pickDEG.edger.R")
}


# paste0(pro,'_kk.up.csv')


# ###3.1的
# edger.final <- merge(edger.covert.name.r, edger.deg.fc.mean[, c("logFC" ,"mean_ctrl" , "mean_ko")], by=0)




# ??????????????????????????????
# 	edger.deg.fc.mean = df1
# 	df1 <- data.frame(df1[,-1], row.names = df1[,1])
# 	df1 = edger.deg.fc.mean


# # Usage
# pickDEG(deseq2)
# pickDEG(edger)
# source(step2_pickDEG.R)

pickDEG_deseq2(deseq2,
	expSet.selected,
	FC_cutoff
	)
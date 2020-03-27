#DEG_function




pickDEG_deseq2 <- function(deseq2,
							expSet.selected,
							FC_cutoff,
							n='DEseq2'){

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





#input
# ensembl_gene=rownames(deseq2.cutoff)
# rankedgene = deseq2.ranked

ensemblToentrezid <- function(ensembl_gene,    ####passed
								rankedgene,
								pro='test'){ # pro= edger or deseq2 
	require("biomaRt")
	mart <- useMart("ENSEMBL_MART_ENSEMBL")
	mart <- useDataset("mmusculus_gene_ensembl", mart)
	colnames(rankedgene) = "logFC"                #colnames(X)[2] <- "superduper"


	ensLook <- gsub("\\.[0-9]*$", "", c(ensembl_gene))
	
	annotLook <- getBM(mart=mart,
							attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name", "entrezgene_id"), #annotLook.down[,4] is entrezid
							filter="ensembl_gene_id",
							values=ensLook,
							uniqueRows=TRUE)                         #num data.frame
			df1 <- annotLook
			df1 <- data.frame(df1[,-1], row.names = df1[,1])
			annotLook.r <- df 

	entrezid <- annotLook[,4]  #  pure integret 

	fc.mean <- merge(rankedgene, expset.mean[, c("mean_ctrl", "mean_ko")], by=0)   #fc.mean rowname is number
			df1 <- fc.mean
			df1 <- data.frame(df1[,-1], row.names = df1[,1])
			fc.mean.r <- df1


	final <- merge(annotLook, fc.mean.r[, c("logFC" ,"mean_ctrl" , "mean_ko")], by=0) ##edger   

	save(final,entrezid,file =paste0('res3_',pro,'_final_DEG.Rdata'))  #rename
}










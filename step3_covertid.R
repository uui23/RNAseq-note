# step3	covertid
# input 	diffgene_cutoff/up/ down (list&data.frame)(enseble) 
# output	diffgene_/up/down (list&data.frame)(with name)(entriz)

# 	step3.1	makeSummaryDEG
# 	input 	means expSet
# 			diffgene_data.frame(entriz)(cutoff)
# 			diffgene_ranked
# 	output 	diffgene_data_summary (name,FC,mean_RC)



	# save(edger,
	# 	edger.ranked,
	# 	edger.cutoff,file="step2_pickDEG.edger.R")


# rm(list = ls()) 
# options(stringsAsFactors = F)

# load('res2_pickDEG.R')
# load('res1_getDEG.R')

# input <- function(pro='test'){ #pro should be deseq2 or edger
# 	if(pro =="deseq2")
# 	ensembl_gene=rownames(deseq2)
# 	rankedgene = edger.ranked
# 	colnames(rankedgene)[2] = "logFC"                #colnames(X)[2] <- "superduper"
# 	if(pro =="edger")
# 	ensembl_gene=rownames(edger)
# 	rankedgene = edger.ranked
# }





ensemblToentrezid <- function(pro='test'){ # pro= edger or deseq2 
	require("biomaRt")
	mart <- useMart("ENSEMBL_MART_ENSEMBL")
	mart <- useDataset("mmusculus_gene_ensembl", mart)


	ensembl_gene=rownames(deseq2.cutoff)
	rankedgene = deseq2.ranked
	colnames(rankedgene) = "logFC"                #colnames(X)[2] <- "superduper"


	# if(pro =="edger")
	# ensembl_gene=rownames(edger)
	# rankedgene = edger.ranked


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

	save(final,entrezid,file =paste0(pro, 'res3_final_DEG.Rdata'))  #rename
}




# write.csv(final, file= "edger.summary.csv")

# "./src/edger.covert.name.csv"

# annotLook <- read.csv("./src/edger.covert.name.csv", row.names = 1)
# rankedgene <-read.csv("./src/edger.rank.deg.csv", row.names = 1)


# #usage

# ensemblToentrezid(deseq2)
# ensemblToentrezid(edger)
#########################################################################################################################


# ensemblToentrezid <- function(pro='test'){ # pro= edger or deseq2 
	
	
# 	ensembl_gene=rownames(deseq2)

# 	rankedgene = deseq2.ranked

# 	colnames(rankedgene)[2] = "logFC"                #colnames(X)[2] <- "superduper"


# 	# if(pro =="edger")
# 	# ensembl_gene=rownames(edger)
# 	# rankedgene = edger.ranked


# 	ensLook <- gsub("\\.[0-9]*$", "", c(ensembl_gene))
	
# 	annotLook <- getBM(mart=mart,
# 							attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name", "entrezgene_id"), #annotLook.down[,4] is entrezid
# 							filter="ensembl_gene_id",
# 							values=ensLook,
# 							uniqueRows=TRUE)                         #num data.frame
# 			df1 <- annotLook
# 			df1 <- data.frame(df1[,-1], row.names = df1[,1])
# 			annotLook.r <- df 

# 	entrezid <- annotLook[,4]  #  pure integret 

# 	fc.mean <- merge(rankedgene, expset.mean[, c("mean_ctrl", "mean_ko")], by=0)   #fc.mean rowname is number
# 			df1 <- fc.mean
# 			df1 <- data.frame(df1[,-1], row.names = df1[,1])
# 			fc.mean.r <- df1


# 	final <- merge(annotLook, fc.mean.r[, c("logFC" ,"mean_ctrl" , "mean_ko")], by=0) ##edger   

# 	save(final,entrezid,file =paste0(pro, 'res3_final_DEG.Rdata'))  #rename
# }


#find TF 

mmu.tf <- read.delim("./src/Mus_musculus_TF.txt")
mmu.tf.co <- read.delim("./src/Mus_musculus_TF_cofactors.txt")

df1 <- mmu.tf
df1 <- data.frame(df1[,-1], row.names = df1[,3])
mmu.tf.r <- df1

df1 <- final
df1 <- data.frame(df1[,-1], row.names = df1[,1])
final.r <- df1

rownames(final.r)
rownames(mmu.tf.r)

edger.tf <- subset(mmu.tf.r, rownames(final.r) %in% rownames(mmu.tf.r))   ### 1k  

edger.tf <- subset(mmu.tf.r,rownames(mmu.tf.r) %in% rownames(final.r))   #### diff here



final.tf <- merge(final.r, edger.tf[, c("Symbol", "Family")], by=0)
#########################

df1 <- mmu.tf.co
df1 <- data.frame(df1[,-1], row.names = df1[,3])
mmu.tf.co.r <- df1

df1 <- final
df1 <- data.frame(df1[,-1], row.names = df1[,1])
final.r <- df1

rownames(final.r)
rownames(mmu.tf.co.r)

edger.tf.co <- subset(mmu.tf.co.r, rownames(final.r) %in% rownames(mmu.tf.co.r))   ### only 8

edger.tf.co <- subset(mmu.tf.co.r,rownames(mmu.tf.co.r) %in% rownames(final.r))   #### 114
#both 114 but 

final.tf.co <- merge(final.r, edger.tf.co[, c("Symbol", "Family")], by=0)

write.csv(final.tf.co, file="./output/edger.tf.co.csv")
write.csv(final.tf, file="./output/edger.tf.csv")


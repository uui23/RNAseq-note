#step1  NE

rm(list = ls()) 
options(stringsAsFactors = F)

genelength <- read.table("./src/mm10.length", stringsAsFactors = F)
exprSet <- read.delim("./src/expSet_htSeq.txt", row.names=1, stringsAsFactors = F)

#manual
exprSet.sub<-exprSet[1:37315,]
exprSet.sub[,5] <- 0

exprSet.sub[genelength[,1],5] <- genelength[,2]

depth<-colSums(exprSet.sub[,-5])
RPKM<-t(apply(exprSet.sub, 1, function(x){x[1:4]/(x[5]*depth)*10^9}))

RPKM.selected <- RPKM[apply(RPKM, 1, function(x){sum(x>1)})==4,]   #matrix
exprSet.rpkm <- RPKM
a <- as.data.frame(RPKM.selected) 
a$mean_ko= (a$ko4 + a$ko5) / 2
a$mean_ctrl = (a$ch + a$c7) / 2
expset.mean <- a

save(expset.mean,exprSet.rpkm,RPKM.selected,file='res1_getDEG.R')


#step2  WO
rm(list = ls()) 
options(stringsAsFactors = F)

source("DEG_function.R")
load("res1_getDEG.R")


deseq2 <- read.delim("./src/deseq2.txt")  #deseq2 get from server

FC_cutoff = 1.5
expSet.selected = RPKM.selected

pickDEG_deseq2(deseq2,
	expSet.selected,
	FC_cutoff
	)


save(deseq2,deseq2.ranked,deseq2.cutoff,file="res2_pickDEG_deseq2.R"
# edger

#step3  WO
rm(list = ls()) 
options(stringsAsFactors = F)

source("DEG_function.R")

load("res1_getDEG.R")
load("res2_pickDEG_deseq2.R")


ensembl_gene=rownames(deseq2.cutoff)
rankedgene = deseq2.ranked

ensemblToentrezid(ensembl_gene,rankedgene,pro='deseq2')

save(final,entrezid,file =paste0('res3_',pro,'_final_DEG.Rdata'))




# step1	getDEG

# input	expSet 
# output	deseq2/edger (data.frame)      
# 		rpkm>1 expSet
# 		means_expSet

# step1


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




load(file='res1_getDEG.R')






check


load(file='res1_getDEG.R')

exprSet=exprSet.rpkm
group_list=colnames(exprSet)
table(group_list)






if(T){
  colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet))
  group_list
  tmp=data.frame(g=group_list)
  rownames(tmp)=colnames(exprSet)
  # 组内的样本的相似性理论上应该是要高于组间的
  # 但是如果使用全部的基因的表达矩阵来计算样本之间的相关性
  # 是不能看到组内样本很好的聚集在一起。
  pheatmap::pheatmap(cor(exprSet),annotation_col = tmp)
  dim(exprSet)
  # 所以我这里初步过滤低表达量基因。

  RPKM.selected <- RPKM[apply(RPKM, 1, function(x){sum(x>1)})==4,] 

  exprSet.test=exprSet[apply(exprSet,1, function(x) sum(x>1) > 4),]

    exprSet.test <- exprSet[apply(exprSet, 1, function(x){sum(x>1)})==4,] 




  dim(exprSet)
  
  exprSet=log(edgeR::cpm(exprSet)+1)
  dim(exprSet)
  # 再挑选top500的MAD值基因
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  M=cor(log2(exprSet+1))
  tmp=data.frame(g=group_list)
  rownames(tmp)=colnames(M)
  pheatmap::pheatmap(M,annotation_col = tmp)
  # 现在就可以看到，组内样本很好的聚集在一起
  # 组内的样本的相似性是要高于组间
  pheatmap::pheatmap(M,annotation_col = tmp,filename = 'cor.png')
  
  
  
  library(pheatmap)
  pheatmap(scale(cor(log2(exprSet+1))))
  
}
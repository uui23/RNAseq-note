pick gene accoding cutoff from deseq2 output


dim(pdata[!is.na(pdata[,2]) & !is.na(pdata[,6]) & pdata[,6]<0.05 & pdata[,2]>log2(1.5),]); #上调基因数目
up.regulat<-rownames(pdata[!is.na(pdata[,2]) & !is.na(pdata[,6]) & pdata[,6]<0.05 & pdata[,2]>log2(1.5),]); #上调基因list



dim(pdata[!is.na(pdata[,2]) & !is.na(pdata[,6]) & pdata[,6]<0.05 & pdata[,2]< -log2(1.5),]); #下调基因数目
down.regulat<-rownames(pdata[!is.na(pdata[,2]) & !is.na(pdata[,6]) & pdata[,6]<0.05 & pdata[,2]< -log2(1.5),]); #下调基因list



write.table(up.regulat, "rpkm_up_regulat.txt", sep="\t", quote=F, row.names=F, col.name=F);
write.table(down.regulat, "rpkm_down_regulat.txt", sep="\t", quote=F, row.names=F, col.name=F);



write.table(annotLookup[,6], "EntrezID_up_regulat.txt", sep="\t", quote=F, row.names=F, col.name=F);

write.table(rownames(RPKM.selected.df), "background", sep="\t", quote=F, row.names=F, col.name=F);

write.csv(rownames(RPKM.selected.df), file = "background"

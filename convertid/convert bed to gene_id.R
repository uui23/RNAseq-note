#convert bed to gene_id


library(rtracklayer)

file=file.choose()

bed=import(file)

library("TxDb.Mmusculus.UCSC.mm10.knownGene")
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

bed.gene = CompGO::annotateBedFromDb( gRanges = bed, db = txdb,window = 5000)
	
bed.gene.list = unique(as.data.frame(bed.gene$gene_id)$value)

write.table(bed.gene.list, "bed.gene.list.txt", sep="\t", quote=F, row.names=F, col.name=F)






BP<-read.table("H:/up-go.txt",head=T,stringsAsFactors = F,sep="\t",quote="");

BP_pvalue <- (-log10(as.numeric(BP[,8])))
names(BP_pvalue)<-BP[,3]

BP_pvalue<-BP_pvalue[1:10];     #####从data_pvalue前50个功能中挑选、合并选出10个function data_pvalue.txt

write.table(BP_pvalue,"H:/dme_BP-10.txt",sep="\t",quote = FALSE,col.names = F);


data<-read.table("H:/dme_BP-10.txt",sep="\t",head=F,stringsAsFactors=F,quote="");



data_pvalue<-data[,2];  #####转化成序列

names(data_pvalue)<-data[,1];

all_pv<-sort(unique(data_pvalue));
col_data_pvalue <- NULL;

for(i in 1:length(data_pvalue))
{
	col_data_pvalue<-c(col_data_pvalue,which(all_pv==data_pvalue[i]));
}


par(mar=c(5,20,4,2),lwd=3);
barplot(rev(data_pvalue),horiz=T,las=1,cex.names=1,xlim=c(0,8),space=0.5,col=2,xlab="-log10(pvalue)",xaxt="n");
axis(1,at=c(0,4,8),lwd=3,tcl=-1);


################################################################


data_pvalue<-gprofile.up[,5];  #up

names(data_pvalue)<-gprofile.up[,2];

all_pv<-sort(unique(data_pvalue));
col_data_pvalue <- NULL;

for(i in 1:length(data_pvalue))
{
	col_data_pvalue<-c(col_data_pvalue,which(all_pv==data_pvalue[i]));
}


par(mar=c(5,20,4,2),lwd=3);
barplot(rev(data_pvalue),border = NA,horiz=T,las=1,cex.names=1,xlim=c(0,20),space=0.5,col="firebrick2",xlab="-log10(pvalue)",xaxt="n");
axis(1,at=c(0,5,10,15,20),lwd=3,tcl=-1);################# firebrick2
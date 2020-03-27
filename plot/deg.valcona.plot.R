# valcona plot
#deg.seleced : RPKM <1 gene list  for deseq2

edger.deg.rpkm.seleced  # edger.deg.rpkm.seleced[,1]


png("test.png")

plot(deg.seleced[,2], -log10(deg.seleced[,6]), pch=16, col="grey80", xlim=c(-4,4), ylim=c(0,10), xaxt="n", yaxt="n", xlab="", ylab="")

points(deg.seleced[deg.seleced[7]=="UP",2], -log10(deg.seleced[deg.seleced[7]=="UP",6]), pch=16, col="firebrick2")
points(deg.seleced[deg.seleced[7]=="DOWN",2], -log10(deg.seleced[deg.seleced[7]=="DOWN",6]), pch=16, col="royalblue")

abline(v=log2(1.5), lwd=3, lty=2, col="grey40");
abline(v=-log2(1.5), lwd=3, lty=2, col="grey40");

abline(h=-log10(0.05), lwd=3, lty=2, col="grey40");

axis(1, at=c(-100, -4, -2, 0, 2, 4, 100), labels=rep("", 7), lwd=3, tcl=-1.2);
axis(2, at=c(-100, 0, 5, 10, 100), labels=rep("", 5), lwd=3, tcl=-1.2);
axis(3, at=c(-100, 100), lwd=3);
axis(4, at=c(-100, 100), lwd=3);

dev.off()


#edger.deg.rpkm.seleced # edger.deg.rpkm.seleced[,1]   [5]

#input  edger.deg.rpkm.seleced  # data.frame   rpkm>1 ## more than 10k 

png("edger.test.png")

plot(edger.deg.rpkm.seleced[,1], -log10(edger.deg.rpkm.seleced[,3]), pch=16, col="grey80", xlim=c(-1,0), ylim=c(1,2), xaxt="n", yaxt="n", xlab="", ylab="")

points(edger.deg.rpkm.seleced[edger.deg.rpkm.seleced[5]=="UP",1], -log10(edger.deg.rpkm.seleced[edger.deg.rpkm.seleced[5]=="UP",3]), pch=16, col="firebrick2")
points(edger.deg.rpkm.seleced[edger.deg.rpkm.seleced[5]=="DOWN",1], -log10(edger.deg.rpkm.seleced[edger.deg.rpkm.seleced[5]=="DOWN",3]), pch=16, col="royalblue")

abline(v=log2(1.5), lwd=3, lty=2, col="grey40");
abline(v=-log2(1.5), lwd=3, lty=2, col="grey40");
abline(h=-log10(0.05), lwd=3, lty=2, col="grey40");

axis(1, at=c(-100, -4, -2, 0, 2, 4, 100), labels=rep("", 7), lwd=3, tcl=-1.2);
axis(2, at=c(-100, 0, 5, 10, 100), labels=rep("", 5), lwd=3, tcl=-1.2);
axis(3, at=c(-100, 100), lwd=3);
axis(4, at=c(-100, 100), lwd=3);

dev.off()

##### zengjingming  
## edger


need_DEG = edger.deg.rpkm.seleced
logFC_cutoff <- log(1.5)
n = "edger"

this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',]))

library(ggplot2)
  g = ggplot(data=need_DEG, 
             aes(x=logFC, y=-log10(PValue), 
                 color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle(this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
  print(g)
  ggsave(g,filename = paste0(n,'_volcano.png'))





choose_matrix.omit    
pheatmap(choose_matrix.omit,
        filename = 'rpkm3_diff_heatmap.png') 
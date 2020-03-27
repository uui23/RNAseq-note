# valcona plot
#deg.seleced : RPKM <1 gene list


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
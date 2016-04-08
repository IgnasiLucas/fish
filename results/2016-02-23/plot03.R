library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(9,"Greens"))(34)
pe.gaps <- as.matrix(read.table("pe/gaps.txt", header=TRUE))
se.gaps <- as.matrix(read.table("se/gaps.txt", header=TRUE))

png(filename="gaps.png", width=2000, height=1000)
par(mfrow=c(1,2), mex=3)
plot(pe.gaps[,1], pe.gaps[,2], type="l", col=colors[11],
   main="Non-merged reads", xlab="Number of gaps", ylab="Frequency",
   lwd=3, cex=3, cex.main=3, cex.lab=3, cex.axis=3)
for (i in 3:14) {
   lines(pe.gaps[,1], pe.gaps[,i], col=colors[9 + i], lwd=3)
}

plot(se.gaps[,1], se.gaps[,2], type="l", col=colors[11],
   main="Merged reads", xlab="Number of gaps", ylab="Frequency",
   lwd=3, cex=3, cex.main=3, cex.lab=3, cex.axis=3)
for (i in 3:14) {
   lines(se.gaps[,1], se.gaps[,i], col=colors[9 + i], lwd=3)
}

dev.off()

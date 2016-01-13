len <- as.matrix(read.table("lengths.txt"))

png(filename="lengths.png", width=1000, height=1000)
par(mex=2, cex=2, cex.axis=2, cex.lab=2, cex.main=2)
hist(len, breaks=20, xlab="Assembly lengths (bp)")
dev.off()

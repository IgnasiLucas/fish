len <- read.table('lengths.txt', header=TRUE)
png(filename="sizes.png", width=3000, height=2000)
par(mfrow=c(4,6), mex=2)
for (i in 2:25) {
   plot(len[250:543,1], len[250:543,i], type="l", lwd=2, xlab="Length (bp)", ylab="# Reads", main=names(len)[i],
      cex.main=3, cex.lab=3, cex.axis=3)
}
dev.off()

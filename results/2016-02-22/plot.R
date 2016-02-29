pe <- read.table("pe/summary.txt", header=TRUE)
se <- read.table("se/summary.txt", header=TRUE)

pe.prop <- pe[,2:14]
se.prop <- se[,2:14]

for (i in 1:24) {
   pe.prop[i,] <- pe.prop[i,]/pe[i,15]
   se.prop[i,] <- se.prop[i,]/se[i,15]
}

png(filename="maxN.png", width= 2000, height=1000)
par(mfrow=c(1,2), mex=2)
plot(seq(from=10, to=70, by=5), pe.prop[1,], type="l", xlab="Maximum Ns", ylab="Num. of reads",
   main="Non-merged reads", cex.main=2, cex.lab=2, cex.axis=2, lwd=2)
for (i in c(2:14,16:24)) {
   lines(seq(from=10, to=70, by=5), pe.prop[i,])
}

plot(seq(from=10, to=70, by=5), se.prop[1,], type="l", xlab="Maximum Ns", ylab="Num. of reads",
   main="Merged reads", cex.main=2, cex.lab=2, cex.axis=2, lwd=2)
for (i in c(2:14, 16:24)) {
   lines(seq(from=10, to=70, by=5), se.prop[i,])
}
dev.off()


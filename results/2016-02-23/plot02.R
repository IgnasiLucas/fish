names <- c("StCa0001", "StCa0003", "StCa0006", "StCa0015", "StCa0016", "StCa0019", "StCf0037", "StCf0039",
           "StCf0043", "StCf0044", "StCf0049", "StCf0050", "BlCa0065", "BlCa0076", "BlCa0080", "BlCa0083",
           "BlCa0104", "BlCa0108", "BlCl0091", "BlCl0093", "BlCl0094", "BlCl0095", "BlCl0098", "BlCl0116")
pe.files <- paste(rep("pe/", 24), names, rep(".derep.Ns", 24), sep="")
se.files <- paste(rep("se/", 24), names, rep(".derep.Ns", 24), sep="")

png(filename="Ns.png", width=6000, height=4000)
par(mfrow=c(4,6), mex=4)
for (i in 1:24) {
   pe <- as.matrix(read.table(pe.files[i]))
   se <- as.matrix(read.table(se.files[i]))
   plot(ecdf(pe[,1]), do.points=FALSE, col="red", xlim=c(0, 0.4), main=names[i], xlab="Prop. Ns", cex.axis=4, cex.lab=4, cex.main=4, lwd=4)
   plot(ecdf(pe[,2]), do.points=FALSE, col="orange", add=TRUE, lwd=4)
   plot(ecdf(se), do.points=FALSE, col="blue", add=TRUE, lwd=4)
   legend(0.2, 0.32, c("P.E. first", "P.E. both", "Merged"), lty=c(1,1,1), lwd=c(2,2,2), col=c("red","orange","blue"), cex=4)
}
dev.off()

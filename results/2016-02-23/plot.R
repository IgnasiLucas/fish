pe <- as.matrix(read.table("pe/summary.txt", header=TRUE))
se <- as.matrix(read.table("se/summary.txt", header=TRUE))

x <- c(74,76,78,80,82,84,86,88,90,92,94,96,98,99)
x.pe <- x[1:(dim(pe)[2]/4)]
x.se <- x[1:(dim(se)[2]/4)]

filenames   <- c("clustNum.png", "meanDepth.png", "clustNumD3.png", "meanDepthD3.png")
main.first  <- c("Non-merged reads", "Non-merged reads", "Non-merged. Depth > 2", "Non-merged. Depth > 2")
main.second <- c("Merged reads", "Merged reads", "Merged reads. Depth > 2", "Merged reads. Depth > 2")
ylabs  <- c("Number of clusters", "Mean depth", "Number of clusters", "Mean Depth")

for (j in 1:4) {
   png(filename=filenames[j], width= 2000, height=1000)
   par(mfrow=c(1,2), mex=2)
   plot(c(73,100), c(min(pe[pe[,j] > 1, seq(from=j, by=4, to=dim(pe)[2])]), max(pe[,seq(from=j, by=4, to=dim(pe)[2])])), 
      type="n", xlab="Clustering threshold", ylab=ylabs[j],
      main=main.first[j], cex.main=2, cex.lab=2, cex.axis=2, lwd=2)
   for (i in 1:24) {
      lines(x.pe, pe[i,seq(from=j, by=4, to=dim(pe)[2])])
   }

   plot(c(73,100), c(min(se[se[,j] > 1, seq(from=j, by=4, to=dim(se)[2])]), max(se[,seq(from=j, by=4, to=dim(se)[2])])), 
      type="n", xlab="Clustering threshold", ylab=ylabs[j],
      main=main.second[j], cex.main=2, cex.lab=2, cex.axis=2, lwd=2)
   for (i in c(1:24)) {
      lines(x.se, se[i,seq(from=j, by=4, to=dim(se)[2])])
   }
   dev.off()
}

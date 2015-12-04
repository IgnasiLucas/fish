# First, let's make a couple of barplots, for the proportions of read pairs
# that were properly merged, not merged, or discarded.

reads <- read.table('numbers.txt', col.names=c("sample", "run", "merged", "unmerged", "discarded", "total"))
attach(reads)
png(filename="merging.png", width=1000, height=1500)
par(mfrow=c(2,1), mex=2, cex.axis=2, cex.lab=2)

barplot(t(as.matrix(reads[run==1,3:4])), names.arg=sample[run==1],
   cex.main=2, xlab="Samples", main="Run 1", ylab="Number of reads")

barplot(t(as.matrix(reads[run==2,3:4])), names.arg=sample[run==2],
   cex.main=2, xlab="Samples", main="Run 2", ylab="Number of reads",
   legend.text = c("Merged", "Not merged"), args.legend=c(x=29,y=10.8e+05,cex=2))
dev.off()

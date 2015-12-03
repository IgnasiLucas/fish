# The purpose is just to make sure that runs do not affect much the distribution
# of reads among samples, which should be straight forward. I do expect different
# runs to have different mean coverage, but the idea is that all samples should
# experience a similar change in coverage between runs. There is no need to fit
# models. I can express this by quantifying the variation in the fold change of
# the number of reads per sample between runs. 
#
# Actually, this is similar to look at the correlation between runs.
#
# The three columns are the sample identifier, the run (1 or 2) and the
# coverage (total number of reads identified from such sample in such run).
# Note that one of the samples is "Undetermined".

Cov <- read.table("coverage.txt", colClasses=c("factor","factor","integer"),
       col.names=c("smpl","run","cov"))

attach(Cov)

# Exploring the data, I notice that the "Undetermined" sample does not fit well
# in the same pattern as the rest: it is the only category with less coverage in
# run 2 than in run 1. Since changes in error rate can make it more difficult in
# one run to properly assign reads to samples, it makes sense to treat those reads
# separately.

f1 <- run == 1 & smpl != "Undetermined"
f2 <- run == 2 & smpl != "Undetermined"

m0 <- lm(cov[f2] ~ cov[f1])
summary(m0)

m1 <- lm(cov[f2] ~ 0 + cov[f1])
summary(m1)

png(filename="coverage.png", width=1000, height=1000)
par(mex=2)
plot(cov[f1], cov[f2], cex=2, cex.axis=2, cex.lab=2, cex.main=2,
   xlab="Coverage in run 1", ylab="Coverage in run 2")
abline(m1, col="red", lwd=2)
dev.off()

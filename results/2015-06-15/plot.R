library("SimRAD")

cs_5p1 <- "T"
cs_3p1 <- "TAA"
cs_5p2 <- "C"
cs_3p2 <- "CGG"
FragLengths <- numeric(length=0)
SalmoLengths <- numeric(length=0)

for (rep in 1:50) {
	salmo <- ref.DNAseq("../../data/2015-06-15/Salmo_salar.fna", 
		subselect.contigs = TRUE,
		prop.contigs = 0.1)
	SalmoLengths <- append(SalmoLengths, length(salmo), after=length(SalmoLengths))
	salmo.dig <- insilico.digest(salmo, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose = TRUE)
	salmo.adapted <- adapt.select(salmo.dig, type = "AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
	salmo.sized <- size.select(salmo.adapted, min.size = 400, max.size = 700, graph = FALSE, verbose = TRUE)
	FragLengths <- append(FragLengths, width(salmo.sized), after=length(FragLengths))
	rm(salmo, salmo.dig, salmo.adapted, salmo.sized)
}

FragFreq <- tabulate(FragLengths)
FragFreq <- FragFreq/5
png(filename = "FragFreq.png")
plot(c(400,700), c(0, max(FragFreq)), type="n", xlab="Fragment size (bp)", ylab="Frequency")
lines(FragFreq)
dev.off()

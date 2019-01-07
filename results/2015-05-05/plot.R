data <- read.table("../../data/2015-05-05/data", sep="\t", header=TRUE, row.names=1)
norm <- data
norm[,3:10] <- scale(data[,3:10]) # This normalizes the numeric columns
attach(norm)
# The names of columns are: "Sp", "Lake", "Mass", "Num.Digest", "Vol.digest", "Vol.total", "Elution", "Ave.Size", "Conc", "Molarity".

case <- norm[is.na(Ave.Size), 1:7]

mod <- lm(cbind(Ave.Size, Conc, Molarity) ~ Sp + Lake + Mass + Num.Digest + Vol.digest + Vol.total + Elution)

predictions <- predict(mod, case)

# These are the values imputed to the samples St0044-2 and Bl0093-2:
# Average Size (bp):
predictions["St0044-2", "Ave.Size"] * sd(data$Ave.Size, na.rm=TRUE) + mean(data$Ave.Size, na.rm=TRUE)
predictions["Bl0093-2", "Ave.Size"] * sd(data$Ave.Size, na.rm=TRUE) + mean(data$Ave.Size, na.rm=TRUE)

# Concentration (note the samples were diluted to 1/3; in pg/µl):
predictions["St0044-2", "Conc"] * sd(data$Conc, na.rm=TRUE) + mean(data$Conc, na.rm=TRUE)
predictions["Bl0093-2", "Conc"] * sd(data$Conc, na.rm=TRUE) + mean(data$Conc, na.rm=TRUE)

# Molarity (pmol/l):
predictions["St0044-2", "Molarity"] * sd(data$Molarity, na.rm=TRUE) + mean(data$Molarity, na.rm=TRUE)
predictions["Bl0093-2", "Molarity"] * sd(data$Molarity, na.rm=TRUE) + mean(data$Molarity, na.rm=TRUE)

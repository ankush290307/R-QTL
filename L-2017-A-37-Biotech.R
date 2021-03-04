
#Importing the library and file on the Rstudio
library(qtl)
sug <- read.cross("csv",dir = "D:/LINUX" , "sug1.csv",genotypes=c("A", "B", "X"), alleles=c("A","B"))
jittermap(sug, amount = 1e-6)
nind(sug)
nchr(sug)
totmar(sug)
nmar(sug)
nphe(sug)
plot(sug)
plotMissing(sug)

plotMap(sug)
plotPheno(sug, pheno.col=1)
plotPheno(sug, pheno.col=2)
plotPheno(sug, pheno.col=5)

#Single QTL analysis
#step1 QTL genotype probabilities
sug <- calc.genoprob(sug, step=1)
#standard interval mapping
out.em <- scanone(sug)
plot(out.em)
summary(out.em)
#chr with LOD > 3
summary(out.em, threshold = 3)
#chr with LOD score more than 2
summary(out.em, threshold = 2)
#genome scan via Haley-Knott regression
out.hk <- scanone(sug, method = "hk")
plot(out.em, out.hk, col = c("red","green"))

#mUltiple Qtl
sug <- calc.genoprob(sug, step=1)
qtl <- makeqtl(sug, chr=c(7,15), pos=c(47.7, 12), what="prob")
out.fq <- fitqtl(sug, qtl=qtl, method="hk") 
rqtl <- refineqtl(sug, qtl=qtl, method= "hk")
plotLodProfile(rqtl)

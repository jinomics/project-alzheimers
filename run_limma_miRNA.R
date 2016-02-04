setwd("/home/sun-gou/Dropbox/Insight/Sage_Biosciences_Consulting_Project/miRNA/")

library(limma)
library(edgeR)

data <- read.csv("mirna_sorted2.csv",header=T,row.names=1)
ref <- read.csv("mirna_clin2.csv",header=T)
rownames(ref) <- ref$projid

# option 1

ref$az <- NA
ref$az[ref$cogdx == 4] <- 1
ref$az[ref$cogdx != 4] <- 0

ref$ok <- NA
ref$ok[ref$cogdx == 1] <- 1
ref$ok[ref$cogdx != 1] <- 0

ref$mci <- NA
ref$mci[ref$cogdx == 2] <- 1
ref$mci[ref$cogdx != 2] <- 0


ref <- subset(ref, select = -c(projid,cogdx))

fit <- lmFit(data,ref)

# To make all pair-wise comparisons between the three groups the appropriate contrast matrix can becreated by
contrast.matrix <- makeContrasts(ok-az, ok-mci, mci-az, levels=ref)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=T)

res <- topTable(fit2, coef=1, adjust="fdr", n = Inf, sort = "p", p = 0.05)

# The outcome of each hypothesis test can be assigned using
results <- decideTests(fit2,adjust.method="fdr",p.value=0.05)

# A Venn diagram showing numbers of genes significant in each comparison can be obtained from
vennDiagram(results)


pattern <- c(0,0,0)
results[results[,1]!=pattern[1] & results[,2]!=pattern[2] & results[,3]!=pattern[3],]

de <- rownames(results[results[,1]!=pattern[1] | results[,2]!=pattern[2] | results[,3]!=pattern[3],])
write(de,"mirna_deg2_fdr0.05.txt")


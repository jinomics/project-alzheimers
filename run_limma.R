setwd("/home/sun-gou/Dropbox/Insight/Sage_Biosciences_Consulting_Project/RNA-seq/")

library(limma)
library(edgeR)

data <- read.csv("fpkm_mrna.csv",header=T,row.names=1)
ref <- read.csv("mrna_clin.csv",header=T)
rna <- read.table("AMP_AD_ROSMAP_Broad-Rush_RNA_Seq_RIN.v3.txt",header=T)
ref <- merge(ref,rna,by="projid")
rownames(ref) <- ref$projid

# option 1

ref$az <- NA
ref$az[ref$cogdx == 4] <- 1
ref$az[ref$cogdx != 4] <- 0

ref$con <- NA
ref$con[ref$cogdx == 1] <- 1
ref$con[ref$cogdx != 1] <- 0

ref$mci<- NA
ref$mci[ref$cogdx == 2] <- 1
ref$mci[ref$cogdx != 2] <- 0


ref <- subset(ref, select = -c(projid,ID,Sampleid,cogdx))

fit <- lmFit(log2(t(data+0.25)),ref)

# To make all pair-wise comparisons between the three groups the appropriate contrast matrix can becreated by
contrast.matrix <- makeContrasts(con-az, con-mci, mci-az, levels=ref)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=T)

res <- topTable(fit2, coef=1, adjust="fdr", n = Inf, sort = "p", p = 0.05)

# The outcome of each hypothesis test can be assigned using
results <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# A Venn diagram showing numbers of genes significant in each comparison can be obtained from
vennDiagram(results)


pattern <- c(0,0,0)
results[results[,1]!=pattern[1] & results[,2]!=pattern[2] & results[,3]!=pattern[3],]

de <- rownames(results[results[,1]!=pattern[1] | results[,2]!=pattern[2] | results[,3]!=pattern[3],])
write(de,"fpkm_log2_mrna_deg_fdr0.05.txt")


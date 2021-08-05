rm(list = ls())
setwd("/Users/reubenthomas/Dropbox (Gladstone)/Projects/BC-SH-Jan18/")
# load("CP-beating-with-Brg1flfl_HopachData.RData")
PhenoData <- read.csv("FullPhenoData.csv", header = TRUE)
Comparison <- "Brm_specific_analyses_W_D4"
# Comparison <- "Brg1flfl"

##read-in the count matrix
PeakCounts <- read.table(paste(Comparison,"_BedEveryMergeCounts_Updated.txt", sep = ""), header = TRUE)

TempIndices <- match(colnames(PeakCounts)[-c(1,2)], PhenoData$FileName)

SubPhenoData <- PhenoData[TempIndices, ]

row.names(PeakCounts) <- as.character(PeakCounts$Geneid)
PeakLengths <- PeakCounts$Length
PeakCounts <- PeakCounts[,-c(1,2)]


##filter regions with low counts
filter <- apply(PeakCounts, 1, function(x) length(x[x>5])>=2)
FilterPeakCounts <- PeakCounts[filter,]



#formal edgeR analyses
require('edgeR')
Condition <- SubPhenoData$Genotype
Condition <- droplevels(Condition)
print(levels(Condition))
Condition <- factor(Condition,levels(Condition)[c(2,1)])
print(levels(Condition))
Stage <- SubPhenoData$Stage
print(levels(Stage))
Stage <- droplevels(Stage)
print(levels(Stage))

Stage <- factor(Stage,levels(Stage)[c(3,2,1)])
print(levels(Stage))
Bmp <- SubPhenoData$BmpLevel
Bmp <- droplevels(Bmp)
print(levels(Bmp))
Bmp <- factor(Bmp,levels(Bmp)[c(2,1)])
print(levels(Bmp))
y <- DGEList(counts=FilterPeakCounts, group=NULL)
y <- calcNormFactors(y)
group <- Condition
design <- model.matrix(~Condition + Stage + Bmp + Condition:Stage + Stage:Bmp +Condition:Bmp)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
Coef <- fit$coefficients
lrt <- glmLRT(fit, coef=2:10)
cpm <- cpm(y)
top <- topTags(lrt, n=dim(FilterPeakCounts)[1])$table

DiffPeaks <- row.names(top[top$FDR<0.05,])

rldMat <- log2(cpm + 1)
rldMat1 <- log2(cpm + 0.01)


write.csv(rldMat, paste(Comparison,"_LogNormalized_Data.csv",sep = ""))
rv <- apply(rldMat1,1,var)
select <- order(rv, decreasing = TRUE)[seq_len(min(10000, length(rv)))]

pca <- prcomp(t(rldMat[select,]), scale=T, center = T)


stddev <- pca$sdev
pc1_var <- round(100*stddev[1]^2/sum(stddev^2))
pc2_var <- round(100*stddev[2]^2/sum(stddev^2))
pc3_var <- round(100*stddev[3]^2/sum(stddev^2))
require('ggplot2')
PlotData <- data.frame(cbind(PC1 = pca$x[,1], PC2 = pca$x[,2]))

GetSequencingBatch <- function(x) {
  return(strsplit(x,"_")[[1]][1])
}
Batch <- sapply(row.names(PlotData), GetSequencingBatch)


PlotData <- cbind(PlotData, Genotype=SubPhenoData$Genotype, Bmp=SubPhenoData$BmpLevel, Stage=SubPhenoData$Stage,  Batch=Batch)
pdf(paste(Comparison,"Genotype_Stage_Updated_PCA_plot.pdf",sep = ""))
sp <- ggplot(PlotData, aes(x=PC1, y=PC2, color=Genotype, shape=Stage)) + geom_point(size=4.5)
print(sp +  xlab(paste("PC1:", pc1_var, "% variance")) + ylab(paste("PC2:", pc2_var, "% variance")) + guides(size=FALSE) + theme(legend.text = element_text(size = rel(1.5)))+ theme(title=element_text(size=rel(1.5))) + theme(axis.text=element_text(size=rel(1.5))))
dev.off()

pdf(paste(Comparison,"Genotype_Bmp_Updated_PCA_plot.pdf",sep = ""))
sp <- ggplot(PlotData, aes(x=PC1, y=PC2, color=Genotype, shape=Bmp)) + geom_point(size=4.5)
print(sp +  xlab(paste("PC1:", pc1_var, "% variance")) + ylab(paste("PC2:", pc2_var, "% variance")) + guides(size=FALSE) + theme(legend.text = element_text(size = rel(1.5)))+ theme(title=element_text(size=rel(1.5))) + theme(axis.text=element_text(size=rel(1.5))))
dev.off()


pdf(paste(Comparison,"Genotype_SequencingBatch_Updated_PCA_plot.pdf",sep = ""))
sp <- ggplot(PlotData, aes(x=PC1, y=PC2, color=Batch, shape=Genotype)) + geom_point(size=4.5)
print(sp +  xlab(paste("PC1:", pc1_var, "% variance")) + ylab(paste("PC2:", pc2_var, "% variance")) + guides(size=FALSE) + theme(legend.text = element_text(size = rel(1.5)))+ theme(title=element_text(size=rel(1.5))) + theme(axis.text=element_text(size=rel(1.5))))
dev.off()
write.csv(PlotData, "Updated_PCA_Plot_Data.csv")
##Output all results
TempIndices <- match(row.names(top), row.names(rldMat1))
FullResults <- data.frame(cbind(PeakId=row.names(top), log2(cpm[TempIndices, order(SubPhenoData$Genotype)]+0.01), top))
write.csv(FullResults, "Updated_cpm_Brm_specific_analyses_w_D4_DiffPeakAnalyses_logCPM.csv", row.names = F)

# Load Packages #######
# Packages ####
# bioconductor
library(oligo)
library(arrayQualityMetrics)
library(pd.hugene.2.0.st)
library(RSQLite)
library(DBI)

library(limma)
library(oligoData)
library(pd.hg.u95av2)
library(sva)
library(gcrma)

# CRAN
library(Rcpp)
library(RSQLite)
library(DBI)
library(VennDiagram)
library(ggplot2)
library(scales)
library(reshape)

sampleInfo <- read.csv("MARBLES_Microarray/MARBLES_micro_covariates_shorten.csv", header=T)

sampleInfo_bTrim1 <- sampleInfo[ which(sampleInfo$fw_timepoint=='bTrim1'),]
table(duplicated(sampleInfo_bTrim1))
sampleInfo_cTrim2 <- sampleInfo[ which(sampleInfo$fw_timepoint=='cTrim2'),]
table(duplicated(sampleInfo_cTrim2))
sampleInfo_dTrim3 <- sampleInfo[ which(sampleInfo$fw_timepoint=='dTrim3'),]
table(duplicated(sampleInfo_dTrim3))
# merge three trimester together
sampleInfo_MB = rbind(sampleInfo_bTrim1, sampleInfo_cTrim2, sampleInfo_dTrim3)
dim(sampleInfo_MB)
# Just take about maternal blood
sampleInfo = sampleInfo_MB
dim(sampleInfo)
head(sampleInfo)

rawData <- read.celfiles(Samples_all)
saveRDS(object = rawData, file = "MARBLES_Microarray/R_objects/rawData_B123.rds")

################################ all samples ################################
# Normalize with RMA
rawData <- readRDS("MARBLES_Microarray/R_objects/rawData_B123.rds")
sumData <- oligo::rma(rawData, target = "core")
saveRDS(object = sumData, file = "MARBLES_Microarray/R_objects/sumData_B123.rds")

# Add Feature Data ####
# Add features
featureData(sumData) <- getNetAffx(sumData, "transcript")

# Remove control probes
tmp <- featureData(sumData)@data$category
keep <- grep("main", tmp)
sumData2 <- sumData[keep,]

# Take out genes without annotation
tmp <- featureData(sumData2)@data$geneassignment
drop <- which(is.na(tmp))
sumData2 <- sumData2[-drop,]

# Output Normalized Expression Data with Annotations
gene <- sumData2
colnames(gene) <- substr(colnames(gene), 1, 17)
genes <- matrix(unlist(lapply(strsplit(featureData(sumData2)@data$geneassignment, split = "//", fixed = TRUE), function(x) paste(x[1:2]))), ncol = 2, byrow = TRUE)
outdat <- cbind(rownames(gene), genes, exprs(gene))
colnames(outdat)[c(1:3)] <- c("Probe", "GeneID", "GeneName")
write.table(outdat, "MARBLES_Microarray/tables/Normalized_MARBLES_Expression_Data_all_123.txt", quote = FALSE, sep = "\t", row.names = FALSE)
rm(gene, genes, keep, tmp, outdat, sumData2, drop)


# Array Quality Metrics ####
# RMA Normalized data
pData(sumData)$Batch <- c(rep("Batch_1", length(Samples_B1)), rep("Batch_2", length(Samples_B2)), rep("Batch_3", length(Samples_B3)))
arrayQualityMetrics(expressionset=sumData, outdir="MARBLES_Microarray/All RMA-Normalized with Batch B123",force=TRUE, do.logtransform = FALSE, intgroup = "Batch")
# outlier: 109069, 100508, 111545, 115554, 105573, 103201, 103653, 113709


# PCA Plots ####
# Setup pData
pData <- pData(sumData)
pData$CEL_FileName <- rownames(pData)
pData <- merge(pData, sampleInfo, by="CEL_FileName")
head(pData)
#pData <- pData[,c("index", "IBC", "CEL_FileName", "Batch.x", "Dx_alg", "COI_GENDER")]
#colnames(pData) <- c("index", "IBC", "CEL_FileName", "Batch", "Dx_alg", "COI_GENDER")
pData(sumData) <- pData

# Make plots
RMAeset <- exprs(sumData)
RMAeset.data<-as.data.frame(RMAeset)
order.RMA<-RMAeset.data[,as.character(pData$CEL_FileName)]
prin <- princomp(order.RMA)
pervar<-prin$sdev^2 / sum(prin$sdev^2)
pervar[1:5]
# 1      2        3        4       5
# 0.976  0.0026   0.0017   0.001   0.0001

loadingsPC1 <- prin$loadings[,"Comp.1"]
thresholdPC1 <- quantile(loadingsPC1, 0.95)
outliersPC1 <- loadingsPC1[loadingsPC1 > thresholdPC1]
names(outliersPC1)

pdf("MARBLES_Microarray/figures/All_RMA Normalized PCA Plots B123.pdf")
screeplot(prin, col="dodgerblue", xlab="Principal Components of RMA Normalized Beta Values All", main="", cex.lab=1.3)

plot.new()
myColors<-c("seagreen3","dodgerblue", "darkorchid", "firebrick1", "darkorange", "khaki1", "azure4", "lightcyan1", "brown", "black")
palette(myColors)
par(mar=c(0, 0, 0, 0))
legend("bottom", legend=levels(as.factor(pData$Batch)), fill=myColors, title="Principal Components by Batch")
pairs(prin$loadings[,1:6], col=as.factor(pData$Batch), labels=c("PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch=1, cex=0.6)

plot.new()
legend("bottom", legend=levels(as.factor(pData$Dx_alg)), fill=myColors, title="Principal Components by Diagnosis")
pairs(prin$loadings[,1:6], col=as.factor(pData$Dx_alg), labels=c("PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch=1, cex=0.6)

plot.new()
legend("bottom", legend=levels(as.factor(pData$COI_GENDER)), fill=myColors, title="Principal Components by Sex")
pairs(prin$loadings[,1:6], col=as.factor(pData$COI_GENDER), labels=c("PC1", "PC2", "PC3", "PC4","PC5", "PC6"), pch=1, cex=0.6)

dev.off()

# KS Test ####
preparedData <- prepdata(expressionset = sumData, do.logtransform = FALSE, intgroup = c("Batch"))
boxplot <- aqm.boxplot(preparedData, subsample = 20000, outlierMethod = "KS")
box.outliers <- boxplot@outliers
box.list <- box.outliers@which
box.stat <- box.outliers@statistic
threshold <- box.outliers@threshold # 0.02282063 # IHer-P1A1-B1_M110835_533201-HuGene2.0_(HuGene-2_0-st).CEL  
write.csv(box.list,"MARBLES_Microarray/tables/All KS boxplot outlier RMA B123.csv")
write.csv(box.stat,"MARBLES_Microarray/tables/All KS boxplot KS stat RMA B123.csv")

# MA diagnostics, d statistic ####
ma <- aqm.maplot(preparedData, subsample = 20000, Dthresh = 0.15, maxNumArrays = 8, nrColumns = 4)
ma.outliers <- ma@outliers
ma.list <- ma.outliers@which #No outliers
ma.stat <- ma.outliers@statistic
ma.threshold <- ma.outliers@threshold #0.15
write.csv(ma.stat,"MARBLES_Microarray/tables/All MA stat RMA B123.csv")

# Upper Quartile Outliers ####
out.upper <- outliers(order.RMA, method = "upperquartile")
stat.upper <- out.upper@statistic
out.threshold <- out.upper@threshold # 6.336067
upper.list <- out.upper@which 

write.csv(stat.upper, "MARBLES_Microarray/tables/All quartile stat RMA B123.csv")

# Prep Data ####
saveRDS(pData, file = "MARBLES_Microarray/R_objects/child_pData_B123.rds")

# NUSE Values and Boxplots (Laptop) ####
# New R session
setwd()
library(oligo)
library(VennDiagram)
memory.limit(50000)
## Not able to do this part
rawData
## fit <- fitProbeLevelModel(rawData)
rm(rawData)
NUSEvalues <- NUSE(fit, type = "values")
write.csv(NUSEvalues, "Tables/All NUSE Values B123.csv")
rm(fit)
NUSEvalues <- as.data.frame(NUSEvalues)
medians <- sapply(NUSEvalues, median, na.rm = TRUE)
names(medians[medians > 1.05])
# No outliers
summary(medians)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9842  0.9915  1.0005  1.0028  1.0105  1.0464 
rm(medians, NUSEvalues)

# There is no samples that failed outlier analysis for all comparisons, so remain all samples in the analysis

# Seperate sumData2 into three trimesters
dim(sumData2)
head(sumData2)
#sumData2 <- readRDS("MARBLES_Microarray/R_objects/sumData2_all.rds")

sumData2_T1 = sumData2[, sumData2$fw_timepoint=="bTrim1"]
dim(sumData2_T1)
table(sumData2_T1$Dx_alg)
sampleInfo_bTrim1 <- sampleInfo[ which(sampleInfo$fw_timepoint=='bTrim1'),]
pData <- pData(sumData2_T1)
pData$CEL_FileName <- row.names(pData)
pData <- merge(pData, sampleInfo_bTrim1, by="CEL_FileName", sort=FALSE, all.x=TRUE, all.y=FALSE)
pData = pData[c(1:16)]
pData <- pData[,c("index", "IBC.x", "CEL_FileName", "Batch.x", "Dx_alg.x", "COI_GENDER.x", "fw_timepoint.x")]
colnames(pData) <- c("index", "IBC", "CEL_FileName", "Batch", "Dx_alg", "COI_GENDER", "fw_timepoint")
colnames(pData)
pData(sumData2_T1) <- pData
pheno <- pData(sumData2_T1)
saveRDS(sumData2_T1, "MARBLES_Microarray/R_objects/sumData2_T1.rds")

sumData2_T2 = sumData2[, sumData2$fw_timepoint=="cTrim2"]
dim(sumData2_T2)
table(sumData2_T2$Dx_alg)
sampleInfo_bTrim1 <- sampleInfo[ which(sampleInfo$fw_timepoint=='cTrim2'),]
pData <- pData(sumData2_T2)
pData$CEL_FileName <- row.names(pData)
pData <- merge(pData, sampleInfo_bTrim1, by="CEL_FileName", sort=FALSE, all.x=TRUE, all.y=FALSE)
pData = pData[c(1:16)]
pData <- pData[,c("index", "IBC.x", "CEL_FileName", "Batch.x", "Dx_alg.x", "COI_GENDER.x", "fw_timepoint.x")]
colnames(pData) <- c("index", "IBC", "CEL_FileName", "Batch", "Dx_alg", "COI_GENDER", "fw_timepoint")
colnames(pData)
pData(sumData2_T2) <- pData
pheno <- pData(sumData2_T2)
saveRDS(sumData2_T2, "MARBLES_Microarray/R_objects/sumData2_T2.rds")

sumData2_T3 = sumData2[, sumData2$fw_timepoint=="dTrim3"]
dim(sumData2_T3)
table(sumData2_T3$Dx_alg)
sampleInfo_bTrim1 <- sampleInfo[ which(sampleInfo$fw_timepoint=='dTrim3'),]
pData <- pData(sumData2_T3)
pData$CEL_FileName <- row.names(pData)
pData <- merge(pData, sampleInfo_bTrim1, by="CEL_FileName", sort=FALSE, all.x=TRUE, all.y=FALSE)
pData = pData[c(1:16)]
pData <- pData[,c("index", "IBC.x", "CEL_FileName", "Batch.x", "Dx_alg.x", "COI_GENDER.x", "fw_timepoint.x")]
colnames(pData) <- c("index", "IBC", "CEL_FileName", "Batch", "Dx_alg", "COI_GENDER", "fw_timepoint")
colnames(pData)
pData(sumData2_T3) <- pData
pheno <- pData(sumData2_T3)
saveRDS(sumData2_T3, "MARBLES_Microarray/R_objects/sumData2_T3.rds")

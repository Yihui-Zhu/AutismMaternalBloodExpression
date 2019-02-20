############# Differential Gene Expression Analysis ############
library(sva)
library(limma)
library(oligo)
library(ggplot2)
library(scales)
library(reshape2)
library(VennDiagram)

### All samples without trimesters difference ####
sumData2_all <- sumData2
pheno_all <- pData(sumData2_all)
# pheno$Dx_alg <- factor(pheno$Dx_alg, levels=c(0, 1, 2), labels = c("TD", "ASD", "NonTD"), ordered=FALSE)
write.table(pheno_all, "MARBLES_Microarray/tables/All_Pheno.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
edata_all <- exprs(sumData2_all)
write.table(edata_all, "MARBLES_Microarray/tables/All_edata.txt", sep="\t", quote=FALSE)

# Differential Expression Analysis with One Model for Diagnosis ####
# Setup model matrices
pheno_all$Dx_alg <- relevel(as.factor(pheno_all$Dx_alg), ref = "TD")
mod <- model.matrix(~Dx_alg, data = pheno_all)
colnames(mod) # "(Intercept)" "Dx_algASD"   "Dx_algNonTD"
colnames(mod) <- c("Intercept", "Dx_alg_ASD", "Dx_alg_NonTD")
mod0 <- model.matrix(~1, data = pheno_all)
colnames(mod0) <- "Intercept"
# edata = as.matrix(edata)

# Run SVA
n.sv <- num.sv(edata_all, mod, method = "be") 
# 28
svobj <- sva(edata_all, mod, mod0, n.sv = n.sv)
sv <- svobj$sv

# Differential Expression Analysis with limma
modSv <- cbind(mod, svobj$sv)
mod0Sv <- cbind(mod0, svobj$sv)
fit <- lmFit(edata_all, modSv)
colnames(modSv)[4:length(colnames(modSv))] <- paste("Sv", 1:28, sep = "")
colnames(mod0Sv)[2:length(colnames(mod0Sv))] <- paste("Sv", 1:28, sep = "")
contrast.matrix <- makeContrasts(Dx_alg_ASD, Dx_alg_NonTD, levels = modSv)
fitContrasts <- contrasts.fit(fit,contrast.matrix)
eb <- eBayes(fitContrasts)

# Make results table
ftable <- topTable(eb, adjust="BH", confint=TRUE, sort.by="none", number=Inf)
ftable$Variance <- apply(edata_all, 1, var)
asd <- topTable(eb, coef="Dx_alg_ASD", adjust="BH", confint=TRUE, sort.by="none", number=Inf)
NonTD <- topTable(eb, coef="Dx_alg_NonTD", adjust="BH", confint=TRUE, sort.by="none", number=Inf)
outTable <- cbind(ftable, asd[,c("CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B")], NonTD[,c("CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B")])
colnames(outTable)[8:13] <- paste("ASD", colnames(outTable)[8:13], sep="_")
colnames(outTable)[14:19] <- paste("NonTD", colnames(outTable)[14:19], sep="_")
write.table(outTable, "MARBLES_Microarray/tables/All_Dx_alg_limma_results_adjSVAonly_B123.txt", sep="\t", quote=FALSE)

# Add GeneNames
genes <- as.data.frame(cbind(as.character(row.names(sumData2_all)), 
                             matrix(unlist(lapply(strsplit(featureData(sumData2_all)@data$geneassignment, split = "//", fixed = TRUE), function(x) paste(x[1:2]))), 
                                    ncol = 2, byrow = TRUE)))
colnames(genes) <- c("Probe", "GeneID", "GeneName")
stats <- outTable
stats$Probe <- row.names(stats)
stats <- merge(stats, genes, by="Probe", all.x=TRUE, all.y=FALSE, sort=FALSE)
stats_All = stats
write.table(stats, "MARBLES_Microarray/tables/All_MARBLES All Genes Stats Dx_alg SVAonly B123.txt", sep = "\t", quote = FALSE)
write.csv(stats_All,"MARBLES_Microarray/tables/stats_All.csv")

############################## bTrim1 #########################################
pheno_T1 <- pData(sumData2_T1)
# pheno$Dx_alg <- factor(pheno$Dx_alg, levels=c(0, 1, 2), labels = c("TD", "ASD", "NonTD"), ordered=FALSE)
write.table(pheno_T1, "MARBLES_Microarray/tables/T1_Pheno.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
edata_T1 <- exprs(sumData2_T1)
write.table(edata_T1, "MARBLES_Microarray/tables/T1_edata_T1.txt", sep="\t", quote=FALSE)

# Differential Expression Analysis with One Model for Diagnosis ####
# Setup model matrices
pheno_T1$Dx_alg <- relevel(as.factor(pheno_T1$Dx_alg), ref = "TD")
mod <- model.matrix(~Dx_alg, data = pheno_T1)
colnames(mod) # "(Intercept)" "Dx_algASD"   "Dx_algNonTD"
colnames(mod) <- c("Intercept", "Dx_alg_ASD", "Dx_alg_NonTD")
mod0 <- model.matrix(~1, data = pheno_T1)
colnames(mod0) <- "Intercept"
# edata = as.matrix(edata)

# Run SVA
n.sv <- num.sv(edata_T1, mod, method = "be") 
# 12
svobj <- sva(edata_T1, mod, mod0, n.sv = n.sv)
sv <- svobj$sv

# Differential Expression Analysis with limma
modSv <- cbind(mod, svobj$sv)
mod0Sv <- cbind(mod0, svobj$sv)
fit <- lmFit(edata_T1, modSv)
colnames(modSv)[4:length(colnames(modSv))] <- paste("Sv", 1:12, sep = "")
colnames(mod0Sv)[2:length(colnames(mod0Sv))] <- paste("Sv", 1:12, sep = "")
contrast.matrix <- makeContrasts(Dx_alg_ASD, Dx_alg_NonTD, levels = modSv)
fitContrasts <- contrasts.fit(fit,contrast.matrix)
eb <- eBayes(fitContrasts)

# Make results table
ftable <- topTable(eb, adjust="BH", confint=TRUE, sort.by="none", number=Inf)
ftable$Variance <- apply(edata_T1, 1, var)
asd <- topTable(eb, coef="Dx_alg_ASD", adjust="BH", confint=TRUE, sort.by="none", number=Inf)
NonTD <- topTable(eb, coef="Dx_alg_NonTD", adjust="BH", confint=TRUE, sort.by="none", number=Inf)
outTable <- cbind(ftable, asd[,c("CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B")], NonTD[,c("CI.L", "CI.R", "t", "P.Value", "adj.P.Val", "B")])
colnames(outTable)[8:13] <- paste("ASD", colnames(outTable)[8:13], sep="_")
colnames(outTable)[14:19] <- paste("NonTD", colnames(outTable)[14:19], sep="_")
write.table(outTable, "MARBLES_Microarray/tables/T1_Dx_alg_limma_results_adjSVAonly_B123.txt", sep="\t", quote=FALSE)

# Add GeneNames
genes <- as.data.frame(cbind(as.character(row.names(sumData2_T1)), 
                             matrix(unlist(lapply(strsplit(featureData(sumData2_T1)@data$geneassignment, split = "//", fixed = TRUE), function(x) paste(x[1:2]))), 
                                    ncol = 2, byrow = TRUE)))
colnames(genes) <- c("Probe", "GeneID", "GeneName")
stats <- outTable
stats$Probe <- row.names(stats)
stats <- merge(stats, genes, by="Probe", all.x=TRUE, all.y=FALSE, sort=FALSE)
stats_T1 = stats
write.table(stats, "MARBLES_Microarray/tables/T1_MARBLES All Genes Stats Dx_alg SVAonly B123.txt", sep = "\t", quote = FALSE)
write.csv(stats_T1,"MARBLES_Microarray/tables/stats_T1.csv")

# Differential Probes ####
# ASD Differential Probes
ASD_diff <- as.character(stats$Probe[stats$ASD_P.Value < 0.01 & abs(stats$Dx_alg_ASD) > 0.1])
ASD_diff <- sort(unique(ASD_diff))
ASD_diff_T1 = ASD_diff
length(ASD_diff) # 376
write.table(x=ASD_diff, file="MARBLES_Microarray/T1_Differentially Expressed Probes Dx_alg B123/MARBLES_ASD_diff_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

ASD_up <- as.character(stats$Probe[stats$ASD_P.Value < 0.01 & stats$Dx_alg_ASD > 0.1])
ASD_up <- sort(unique(ASD_up))
ASD_up_T1 = ASD_up
length(ASD_up) # 196
write.table(x=ASD_up, file="MARBLES_Microarray/T1_Differentially Expressed Probes Dx_alg B123/MARBLES_ASD_up_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

ASD_down <- as.character(stats$Probe[stats$ASD_P.Value < 0.01 & stats$Dx_alg_ASD < -0.1])
ASD_down <- sort(unique(ASD_down))
ASD_down_T1 = ASD_down
length(ASD_down) # 180
write.table(x=ASD_down, file="MARBLES_Microarray/T1_Differentially Expressed Probes Dx_alg B123/MARBLES_ASD_down_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

# NonTD Differential Probes
NonTD_diff <- as.character(stats$Probe[stats$NonTD_P.Value < 0.01 & abs(stats$Dx_alg_NonTD) > 0.1])
NonTD_diff <- sort(unique(NonTD_diff))
NonTD_diff_T1 = NonTD_diff
length(NonTD_diff) # 364
write.table(x=NonTD_diff, file="MARBLES_Microarray/T1_Differentially Expressed Probes Dx_alg B123/MARBLES_NonTD_diff_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

NonTD_up <- as.character(stats$Probe[stats$NonTD_P.Value < 0.01 & stats$Dx_alg_NonTD > 0.1])
NonTD_up <- sort(unique(NonTD_up))
NonTD_up_T1 = NonTD_up
length(NonTD_up) # 204
write.table(x=NonTD_up, file="MARBLES_Microarray/T1_Differentially Expressed Probes Dx_alg B123/MARBLES_NonTD_up_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

NonTD_down <- as.character(stats$Probe[stats$NonTD_P.Value < 0.01 & stats$Dx_alg_NonTD < -0.1])
NonTD_down <- sort(unique(NonTD_down))
length(NonTD_down_T1) # 160
NonTD_down_T1 = NonTD_down
write.table(x=NonTD_down, file="MARBLES_Microarray/T1_Differentially Expressed Probes Dx_alg B123/MARBLES_NonTD_down_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

# All probes
all <- sort(unique(as.character(stats$Probe)))
length(all)
all_T1 = all
write.table(x=all, file="MARBLES_Microarray/T1_Differentially Expressed Probes Dx_alg B123/MARBLES_all_probes_Dx_alg_B123.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)

## Process this analysis for each trimesters


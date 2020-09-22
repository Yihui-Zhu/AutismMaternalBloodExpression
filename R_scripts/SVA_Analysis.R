# Packages ####
sapply(c("tidyverse", "scales", "variancePartition", "sm", "biomaRt", "reshape2", "wesanderson", "WGCNA", "matrixStats"), 
       require, character.only = TRUE)

# Differential Expression Analysis with One Model for Diagnosis ####
# Setup model matrices
pheno_all$Dx_alg <- relevel(as.factor(pheno_all$Dx_alg), ref = "TD")
mod <- model.matrix(~Dx_alg, data = pheno_all)
colnames(mod) <- c("Intercept", "Dx_alg_ASD", "Dx_alg_NonTD")
mod0 <- model.matrix(~1, data = pheno_all)
colnames(mod0) <- "Intercept"

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
write.table(outTable, "SVAadjusted.txt", sep="\t", quote=FALSE)

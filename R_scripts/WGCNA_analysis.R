##########  WGCNA analysis  ########
# This script is adapted from WGCNA package 
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Langfelder, Peter, and Steve Horvath. "WGCNA: an R package for weighted correlation network analysis." BMC bioinformatics 9.1 (2008): 559.
# read in R package
library(MASS)
library(class)
library(cluster)
library(impute)
library(WGCNA)
options(stringsAsFactors = F)


# All_combat_edata don't have batch effect
dim(All_combat_edata)
typeof(All_combat_edata)
All_combat_edata_B = All_combat_edata
All_combat_edata_B = as.data.frame(All_combat_edata_B)

All_combat_edata_B$rn <- rownames(All_combat_edata_B)

genes$rn <- genes$Probe

All_combat_edata_B_all <- join_all(list(genes,All_combat_edata_B), by = 'rn', type = 'full')
dim(All_combat_edata_B_all)

row.names(All_combat_edata_B_all) = All_combat_edata_B_all$rn
dim(All_combat_edata_B_all)
head(All_combat_edata_B_all)
rownames(All_combat_edata_B_all)
All_combat_edata_B_all$rn <- NULL
colnames(All_combat_edata_B_all)

head(All_combat_edata_B_all)
dim(genes)

write.csv(All_combat_edata_B_all, "MARBLES_Microarray/tables/All_combat_edata_B_all.csv")

# this contains information on the genes
# the general information about the data frame, gene symbol and location, color
dat0 = All_combat_edata_B_all
datSummary = dat0[,1:3]
dim(dat0)

# the following data farame contains
# the gene expression data: columns are genes, rows are arrays(samples)
datExpr = t(dat0[,4:303])
head(datExpr)
colnames(datExpr)
rownames(datExpr)

no.samples = dim(datExpr)[[1]] # number of samples
dim(datExpr)

library(WGCNA)
options(stringsAsFactors = FALSE)
dim(sampleInfo)
All_combat_edata_B_all
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
# Take a quick look at what is in the data set:
dim(All_combat_edata_B_all);
names(All_combat_edata_B_all);

gene_annotation = All_combat_edata_B_all[,c(1:3)]
datExpr0 = as.data.frame(t(All_combat_edata_B_all[, -c(1:3)]));
names(datExpr0) = All_combat_edata_B_all$Probe;
rownames(datExpr0) = names(All_combat_edata_B_all)[-c(1:3)];
colnames(datExpr0)

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
dim(datExpr)

# write.csv(sampleInfo,"MARBLES_Microarray/sampleInfo.csv")
sampleInfo = read.csv("MARBLES_Microarray/sampleInfo.csv", header=TRUE)
# M: 1, F: 2
# TD: 1, ASD: 2, NonTD: 3
traitData = sampleInfo
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -c(2,3,4,7,8,9,11,13,14)];
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

MicroSamples = rownames(datExpr);
traitRows = match(MicroSamples, allTraits$CEL_FileName);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "MARBLES_Microarray/R_objects/Micro_WGCNA_dataInput.RData")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "MicroAdjustTOM", 
                       verbose = 3)
table(net$colors)
net$MEs

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "MARBLES_Microarray/ME_Adjusted_networkConstruction-auto.RData")

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable weight containing the weight column of datTrait
# weight here is diagnosis
weight = as.data.frame(datTraits$Dx_alg);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

table(moduleColors)
module = "orange"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Diagnosis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(datExpr)

modNames
names(datExpr)[moduleColors=="brown"]

annot = gene_annotation
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$Probe)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$GeneName[probes2annot],
                       LocusLinkID = annot$GeneID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor);
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "MARBLES_Microarray/geneInfo.csv")

hubs    =
  chooseTopHubInEachModule(
    datExpr, 
    moduleColors, 
    omitColors = "grey", 
    power = 2, 
    type = "signed")
hubs
length(hubs)
write.csv(hubs, "MARBLES_Microarray/hubs.csv")
genes

hubs_genes = as.data.frame(hubs)
hubs_genes$module = rownames(hubs_genes)
typeof(hubs_genes$hubs)
genes$Probe
hubs_genes_2 = merge(hubs_genes, genes, by.x=c("hubs"), by.y =c("Probe"))
write.csv(hubs_genes_2,"MARBLES_Microarray/hubs_genes_2.csv")

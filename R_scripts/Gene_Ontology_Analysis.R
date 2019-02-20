## Gene Ontology Analysis for each trimester separately 
# Packages ####
library(sva)
library(limma)
library(oligo)
library(ggplot2)
library(scales)
library(reshape2)
library(VennDiagram)

## Trimester 1
all_probes <- as.character(T1_probes$Probe)
all_genes <- sort(unique(as.character(T1_probes$GeneName)))

ASD_T1_diff <- subset(T1_probes, abs(T1_probes$Dx_alg_ASD) > 0.1 & ASD_P.Value < 0.01)
dim(ASD_T1_diff)
ASD_T1_up <- subset(ASD_T1_diff, Dx_alg_ASD > 0)
ASD_T1_down <- subset(ASD_T1_diff, Dx_alg_ASD < 0)

ASD_T1_diff_probes <- as.character(ASD_T1_diff$Probe)
ASD_T1_up_probes <- as.character(ASD_T1_up$Probe)
ASD_T1_down_probes <- as.character(ASD_T1_down$Probe)

ASD_T1_diff_genes <- as.character(ASD_T1_diff$GeneName)
ASD_T1_up_genes <- as.character(ASD_T1_up$GeneName)
ASD_T1_down_genes <- as.character(ASD_T1_down$GeneName)

# T1 NonTD
NonTD_T1_diff <- subset(T1_probes, abs(T1_probes$Dx_alg_NonTD) > 0.1 & NonTD_P.Value < 0.01)
dim(NonTD_T1_diff)
NonTD_T1_up <- subset(NonTD_T1_diff, Dx_alg_NonTD > 0)
NonTD_T1_down <- subset(NonTD_T1_diff, Dx_alg_NonTD < 0)

NonTD_T1_diff_probes <- as.character(NonTD_T1_diff$Probe)
NonTD_T1_up_probes <- as.character(NonTD_T1_up$Probe)
NonTD_T1_down_probes <- as.character(NonTD_T1_down$Probe)

NonTD_T1_diff_genes <- as.character(NonTD_T1_diff$GeneName)
NonTD_T1_up_genes <- as.character(NonTD_T1_up$GeneName)
NonTD_T1_down_genes <- as.character(NonTD_T1_down$GeneName)

genes_array <-list("ASD_T1_diff_genes"=ASD_T1_diff_genes, "ASD_T1_up_genes"=ASD_T1_up_genes, "ASD_T1_down_genes"=ASD_T1_down_genes, 
                   "NonTD_T1_diff_genes"=NonTD_T1_diff_genes, "NonTD_T1_up_genes"=NonTD_T1_up_genes, 
                   "NonTD_T1_down_genes"=NonTD_T1_down_genes)
genes_array <- lapply(genes_array, function(x){gsub(" ", "", x, fixed = TRUE)})

# Get Gene IDs from Ensembl ####
probes <- list("ASD_T1_diff_probes"=ASD_T1_diff_probes, "ASD_T1_up_probes"=ASD_T1_up_probes, "ASD_T1_down_probes"=ASD_T1_down_probes, 
               "NonTD_T1_diff_probes"=NonTD_T1_diff_probes, "NonTD_T1_up_probes"=NonTD_T1_up_probes, 
               "NonTD_T1_down_probes"=NonTD_T1_down_probes, "all_probes"=all_probes)

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="apr2018.archive.ensembl.org", dataset="hsapiens_gene_ensembl", ensemblRedirect = FALSE, verbose=TRUE, archive=FALSE)
genes <- sapply(probes, function(x){getBM(attributes="hgnc_symbol", filters="affy_hugene_2_0_st_v1", values=x, mart=ensembl, verbose=TRUE)})
length(genes$ASD_T1_diff_probes.hgnc_symbol)
length(probes$ASD_T1_diff_probes)

write.table(sort(unique(genes$ASD_T1_diff_probes.hgnc_symbol)),
            file="T1_ASD_Differential_Genes_ensembl_B123.txt", sep="\t", quote=FALSE, row.names=FALSE, 
            col.names=FALSE)
write.table(sort(unique(genes$NonTD_T1_diff_probes.hgnc_symbol)),
            file="T1_NonTD_Differential_Genes_ensembl_B123.txt", sep="\t", quote=FALSE, row.names=FALSE, 
            col.names=FALSE)
# match probes to gene list
## Repeat for all three trimesters

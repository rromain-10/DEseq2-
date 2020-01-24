#######################################

#__Date:__ December
# #__Author:__ Erin Osborne Nishimura
# #__Script:__ 191210_GomezOrte.R
# #__Project:__ To analyze RNA-seq in a project that compares yeast grown in rich media v. acetic acid media at different time points
# #__Requires:__ 
# # 
# # + DESeq2 (1.20.0)   
# # + corrplot (0.84)
# # + RColorBrewer (1.1.2)
# # + pheatmap (1.0.10)
# 
# ######################################
# 
# ######### FOR FIRST TIME USE ONLY ##############
# ######### After use, comment this section ##############
# 
# # Intall required packages/libraries:

#Install Bioconductor R version 3.5 or less
#source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

#OR

# Install bioconductor R version 3.6 or greater
#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install(version = "3.10")

## When prompted to install (a/s) all or some of the dependent packages, move your cursor down to the console and type "a" for all

# Install DESeq2:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
# #BiocManager::install("DESeq2")

# Install 'apeglm'
#BiocManager::install("apeglm")

# install corrplot:
#install.packages("corrplot")

#Install pretty heatmap - pheatmap: https://cran.r-project.org/web/packages/pheatmap/index.html
#install.packages("pheatmap")

#Install RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
#install.packages("RColorBrewer")
# 
# ######### DONE WITH: FIRST TIME USE ONLY SECTION ##############
# ################################


###########  LOAD PACKAGE  #####################

# Do this every time.
# Load packages:
library(DESeq2)
library(corrplot)
library(RColorBrewer)
library(pheatmap)

#################################################




###########  READ IN THE DATA  #####################

#Must be in your code repository directory, and the counts and metadata file must be in that directory:
getwd()

# import the counts data
countsData <- read.table(file = "20200116_robs_counts.txt", header = TRUE, row.names = 1, skip = 2) 
head(countsData)
dim(countsData)
class(countsData)

# Read in the metadata
metadata <- read.table(file = "RR_ARPE_Deluca_Collab_manifest.txt", header = FALSE) # import the data
metadata
colnames(metadata) <- c("fasta1", "fasta2", "names", "rep","sample")
metadata

# Organize the countsData file.
# Notice that the countsData file doesn't have any column headers:
head(countsData)

# The last names will be names for each sample. We can pull those names from metadata1:
as.vector(metadata$names)
as.vector(metadata$sample)

# Name countsData columns headers:
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata$names))
colnames(countsData)
head(countsData)

################### COUNT MATRIX INPUT ###################
# In this section we will prepare our input data for analysis.
# In the instruction for DESeq2, it states: "We read in a count matrix, which we will name cts, and the sample information table, which we will name coldata."

# generate a table called "cts" out of the countsData table.
# Subset the countsData 
head(countsData)
dim(countsData)
head(countsData[,6:22])
# Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:22])
head(cts)

# Next we need to make an information called coltable. We can make this out of the metadata table.
class(metadata)
# Reorganize the metadata table so the names1 column are now row headers
metadata
rownames(metadata)<- metadata$names
metadata

coldata <- metadata[,c("sample", "rep", "names")]
#coldata = coldata[order(coldata$names),]
coldata
coldata$sample <- as.factor(coldata$sample)
coldata$names <- as.factor(coldata$names)
coldata$rep <- as.factor(coldata$rep)
rownames(coldata) <- metadata$names
coldata$sample
coldata$names
colnames(cts)


# One thing we need to explicitly check. The rownames of coldata need to exactly match the colnames of cts.
#Check that the names match  --> Should be TRUE
all(rownames(coldata) == colnames(cts))

#BiocManager::install("ensembldb")
#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
library(ensembldb)
library(biomaRt)
ensembl = useDataset("hsapiens_gene_ensembl",mart=useMart("ENSEMBL_MART_ENSEMBL"))
ctsnames = rownames(cts)

# figure out filter,values,and attributes at http://www.ensembl.org/biomart/
gene_names = getBM(mart = ensembl, filter='ensembl_gene_id', value=ctsnames, attributes=c('external_gene_name', 'ensembl_gene_id'))
dim(gene_names) # 20002     2
length(ctsnames) # 20003
ctsnames = ctsnames[ ctsnames %in%  gene_names$ensembl_gene_id] # only one gene not found- ENSG00000254462, and it's all 0s
length(ctsnames) # 20002
rownames(gene_names) <- gene_names$ensembl_gene_id
gene_names = gene_names[ctsnames,]
cts=cts[rownames(cts) != 'ENSG00000254462',]# take out the 0 count gene
rownames(cts) <- gene_names$external_gene_name
head(cts)

########## dds ################
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ sample + rep)
                             
################# PRE-FILTERING -- FILTER FOR PRESENT GENES: ################# 
# Not necessary, but helps keep things fast. 
# Exclude all samples that have 0 reads:

# keep <- rowSums(counts(dds1)) >= 1
# dds1 <- dds1[keep,]

keep <- rowSums(counts(dds)) >=1
dds <- dds[keep,]

# How many did we exclude?
dim(dds)
str(dds)
coldata

################### NOTE ON FACTOR LEVELS ###################
# Organize the categories based on what makes sense:
dds$names<- factor(dds$names)
dds$names

# PERFORM DESEQ2 analysis:
####### This is where the magic happens! ########### 
# This will transform dds into a specialized object with many more fields filled in.
# dds1<- DESeq(dds1)
dds<- DESeq(dds)
class(dds)
str(dds)
dim(dds)
plotDispEsts(dds)
head(dds)

# Here is a demonstration of the size Factor scaling that was calculated (sizeFactor):
# dds1$sizeFactor
# head(counts(dds, normalized = TRUE))
# head(counts(dds, normalized = FALSE))
dds$sizeFactor
head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))
normalized_genecounts <-counts(dds, normalized = TRUE)
#save normalized cts file to directory
#write.csv(normalized_genecounts, "/Users/romarioromain/OneDrive - Colostate/RR_ARPE_DELUCA_COLLAB/DEseq/DEseq2\\2020021_normalized_genecounts.csv")

#Take r-stabilized log transformations of all the normalized count data. This will help with the problem that the data is noisy and it will help with the problem that the data is spread across a wide range of values.
rld <- rlog(dds, blind=FALSE)  #Take the r-stabilized log transformed data:

# Calculate the distances between each sample
sampleDists <- dist(t(assay(rld))) # calculate distance matrices:
sampleDistMatrix <- as.matrix(sampleDists) #convert from data.frame -> matrix
rownames(sampleDistMatrix) <- colnames(rld) # Add some labels
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #Pick some pretty colors

# Draw the heatmap
par(mfrow=c(1,1))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, width = 4, height = 4,
         cluster_rows = TRUE) # Plot the heatmap

############## PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION #####################

# An input requirement of the lfcShrink function is a coef term. This is pulled from the resultsNames of dds:
# resultsNames(dds1)
resultsNames(dds)


############## DIFFERENTIAL EXPRESSION ANALYSIS #####################
# calculate the statistically significant differences between resultNames(dds)

######ARPE19 vs Aktmyr#########
res_ARPE19_vs_Aktmyr <- results(dds,
                                lfc = 1.2,
                                contrast=c("sample","ARPE19","Aktmyr"))

summary(res_ARPE19_vs_Aktmyr)

resLFC_ARPE19_vs_Aktmyr <- lfcShrink(dds,
                                     coef="sample_ARPE19_vs_Aktmyr", type='apeglm')

summary(resLFC_ARPE19_vs_Aktmyr)

##########ARPE19 vs Aktmyr plots#########
par(mfrow=c(1,1))
plotMA(res_ARPE19_vs_Aktmyr, main="ARPE19 vs neg Aktmyr\nunshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized ARPE19 / Aktmyr)",
       xlab = "means of normalized counts")

plotMA(resLFC_ARPE19_vs_Aktmyr, main="ARPE19 vs Aktmyr\nshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized ARPE19 / Aktmyr)",
       xlab = "means of normalized counts")

#check known genes (enter gene names where after the ..=="xxx")
plotCounts(dds, gene=which(rownames(resLFC_ARPE19_vs_Aktmyr)=="BUB1"),intgroup = c("sample"))
plotCounts(dds, gene=which(rownames(resLFC_ARPE19_vs_Aktmyr)=="PPP2R5C"),intgroup = c("sample"))

#Identify genes on the plot ARPE19 vs Aktmyr
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_ARPE19_vs_Aktmyr$baseMean,resLFC_ARPE19_vs_Aktmyr$log2FoldChange, labels = gene_names$external_gene_name)
#  Step4 -> click here to see what you got!
rownames(resLFC_ARPE19_vs_Aktmyr)[idx]

#########RasV12 vs Aktmyr#############
res_RasV12_vs_Aktmyr <- results(dds,
                                lfc = 1.2,
                                contrast=c("sample","RasV12","Aktmyr"))

summary(res_RasV12_vs_Aktmyr)

resLFC_RasV12_vs_Aktmyr <- lfcShrink(dds,
                                     coef="sample_RasV12_vs_Aktmyr", type='apeglm')

summary(resLFC_RasV12_vs_Aktmyr)

########RasV12 vs Aktmyr plots##########
par(mfrow=c(1,1))
plotMA(res_RasV12_vs_Aktmyr, main="RasV12 vs Aktmyr\nunshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized RasV12 / Aktmyr)",
       xlab = "means of normalized counts")

plotMA(resLFC_RasV12_vs_Aktmyr, main="RasV12 vs Aktmyr\nshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized RasV12 / Aktmyr)",
       xlab = "means of normalized counts")

#check known genes (enter gene names where after the ..=="xxx")
plotCounts(dds, gene=which(rownames(resLFC_RasV12_vs_Aktmyr)=="BUB1"),intgroup = c("rep","sample"))
plotCounts(dds, gene=which(rownames(resLFC_RasV12_vs_Aktmyr)=="BUB1"),intgroup = c("sample"))

#Identify genes on the plot ARPE19 vs Aktmyr
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_RasV12_vs_Aktmyr$baseMean,resLFC_RasV12_vs_Aktmyr$log2FoldChange, labels = gene_names$external_gene_name)
#  Step4 -> click here to see what you got!
rownames(resLFC_RasV12_vs_Aktmyr)[idx]

#MekDD vs Aktmyr
res_MekDD_vs_Aktmyr <- results(dds,
                                lfc = 1.2,
                                contrast=c("sample","MekDD","Aktmyr"))

summary(res_MekDD_vs_Aktmyr)

resLFC_MekDD_vs_Aktmyr <- lfcShrink(dds,
                                     coef="sample_MekDD_vs_Aktmyr", type='apeglm')

summary(resLFC_MekDD_vs_Aktmyr)

##########MekDD vs Aktmyr plots##########
par(mfrow=c(1,1))
plotMA(res_MekDD_vs_Aktmyr, main="MekDD vs Aktmyr\nunshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized MekDD / Aktmyr)",
       xlab = "means of normalized counts")

plotMA(resLFC_MekDD_vs_Aktmyr, main="MekDD vs Aktmyr\nshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized MekDD / Aktmyr)",
       xlab = "means of normalized counts")

#Identify genes on the plot ARPE19 vs Aktmyr
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_MekDD_vs_Aktmyr$baseMean,resLFC_MekDD_vs_Aktmyr$log2FoldChange, labels = gene_names$external_gene_name)
#  Step4 -> click here to see what you got!
rownames(resLFC_MekDD_vs_Aktmyr)[idx]


#T53D4 vs Aktmyr
res_T53D4_vs_Aktmyr <- results(dds,
                               lfc = 1.2,
                               contrast=c("sample","T53D4","Aktmyr"))

summary(res_T53D4_vs_Aktmyr)

resLFC_T53D4_vs_Aktmyr <- lfcShrink(dds,
                                    coef="sample_T53D4_vs_Aktmyr", type='apeglm')

summary(resLFC_T53D4_vs_Aktmyr)

########T53D4 vs Aktmyr#########
par(mfrow=c(1,1))
plotMA(res_T53D4_vs_Aktmyr, main="T53D4 vs Aktmyr\nunshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized T53D4 / Aktmyr)",
       xlab = "means of normalized counts")

plotMA(resLFC_T53D4_vs_Aktmyr, main="T53D4 vs Aktmyr\nshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized T53D4 / Aktmyr)",
       xlab = "means of normalized counts")

plotCounts(dds, gene=which(rownames(resLFC_T53D4_vs_Aktmyr)=="BUB1"),intgroup = c("sample"))

#Identify genes on the plot ARPE19 vs Aktmyr
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_T53D4_vs_Aktmyr$baseMean,resLFC_T53D4_vs_Aktmyr$log2FoldChange, labels = gene_names$external_gene_name)
#  Step4 -> click here to see what you got!
rownames(resLFC_T53D4_vs_Aktmyr)[idx]

#negclt(water) vs Aktmyr
res_Water_vs_Aktmyr <- results(dds,
                               lfc = 1.2,
                               contrast=c("sample","negclt","Aktmyr"))

summary(res_Water_vs_Aktmyr)

resLFC_Water_vs_Aktmyr <- lfcShrink(dds,
                                    coef="sample_negclt_vs_Aktmyr", type='apeglm')

summary(resLFC_Water_vs_Aktmyr)

########Water vs Aktmyr#########
par(mfrow=c(1,1))
plotMA(res_Water_vs_Aktmyr, main="Water vs Aktmyr\nunshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized Water / Aktmyr)",
       xlab = "means of normalized counts")

plotMA(resLFC_Water_vs_Aktmyr, main="Water vs Aktmyr\nshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized Water / Aktmyr)",
       xlab = "means of normalized counts")

plotCounts(dds, gene=which(rownames(resLFC_Water_vs_Aktmyr)=="BUB1"),intgroup = c("sample"))

#Identify genes on the plot ARPE19 vs Aktmyr
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_Water_vs_Aktmyr$baseMean,resLFC_Water_vs_Aktmyr$log2FoldChange, labels = gene_names$external_gene_name)
#  Step4 -> click here to see what you got!
rownames(resLFC_Water_vs_Aktmyr)[idx]


##### VOLCANO PLOTS:###################
# Volcano plots are nice ways of displaying the fold change against the p-value.
resultsNames(dds)

#ARPE19 vs Aktmyr
res_volcano_plot_ARPE19vsAktmyr <- lfcShrink(dds,coef="sample_ARPE19_vs_Aktmyr", type='apeglm') # Calculate padj without a lower lfc limit
significantLFC_volcano_plot_ARPE19vsAktmyr <- subset(res_volcano_plot_ARPE19vsAktmyr, padj < 0.05) # Identify significantly changing genes
significant_points_to_plot_ARPE19vsAktmyr <-res_volcano_plot_ARPE19vsAktmyr [which(rownames(res_volcano_plot_ARPE19vsAktmyr ) %in% rownames(significantLFC_volcano_plot_ARPE19vsAktmyr)),] 
# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout_ARPE19vsAktmyr <- subset(res_volcano_plot_ARPE19vsAktmyr , padj < 10e-100)

#Draw volcano plot for ARPE19 vs Aktmyr
par(mfrow=c(1,1)) # one plot only 
# Draw the plot
plot(res_volcano_plot_ARPE19vsAktmyr $log2FoldChange, -log10(res_volcano_plot_ARPE19vsAktmyr $padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-6,6), col = "#00000030")
# Add points
points(maxedout_ARPE19vsAktmyr$log2FoldChange, rep(102, length(maxedout_ARPEvsAktmyr$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot_ARPE19vsAktmyr$log2FoldChange, -log10(significant_points_to_plot_ARPE19vsAktmyr$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")
# Add lines
abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")

#RasV12 vs Aktmyr
res_volcano_plot_RasV12vsAktmyr <- lfcShrink(dds,coef="sample_RasV12_vs_Aktmyr", type='apeglm') # Calculate padj without a lower lfc limit
significantLFC_volcano_plot_RasV12vsAktmyr <- subset(resLFC_RasV12_vs_Aktmyr, padj < 0.05) # Identify significantly changing genes
significant_points_to_plot_RasV12vsAktmyr <-res_volcano_plot_RasV12vsAktmyr [which(rownames(res_volcano_plot_RasV12vsAktmyr ) %in% rownames(significantLFC_volcano_plot_RasV12vsAktmyr)),] 
# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout_RasV12vsAktmyr <- subset(res_volcano_plot_RasV12vsAktmyr , padj < 10e-100)

#Draw volcano plot for ARPE19 vs Aktmyr
par(mfrow=c(1,1)) # one plot only 
plot(res_volcano_plot_RasV12vsAktmyr $log2FoldChange, -log10(res_volcano_plot_RasV12vsAktmyr $padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-6,6), col = "#00000030")
# Add points
points(maxedout_RasV12vsAktmyr$log2FoldChange, rep(102, length(maxedout_RasV12vsAktmyr$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot_RasV12vsAktmyr$log2FoldChange, -log10(significant_points_to_plot_RasV12vsAktmyr$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")
# Add lines
abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")


#MekDD vs Aktmyr
res_volcano_plot_MekDDvsAktmyr <- lfcShrink(dds,coef="sample_MekDD_vs_Aktmyr", type='apeglm') # Calculate padj without a lower lfc limit
significantLFC_volcano_plot_MekDDvsAktmyr <- subset(resLFC_MekDD_vs_Aktmyr, padj < 0.05) # Identify significantly changing genes
significant_points_to_plot_MekDDvsAktmyr <-res_volcano_plot_MekDDvsAktmyr [which(rownames(res_volcano_plot_MekDDvsAktmyr ) %in% rownames(significantLFC_volcano_plot_MekDDvsAktmyr)),] 
# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout_MekDDvsAktmyr <- subset(res_volcano_plot_MekDDvsAktmyr , padj < 10e-100)

#Draw volcano plot for MekDD vs Aktmyr
par(mfrow=c(1,1)) # one plot only 
plot(res_volcano_plot_MekDDvsAktmyr $log2FoldChange, -log10(res_volcano_plot_MekDDvsAktmyr $padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-6,6), col = "#00000030")
# Add points
points(maxedout_MekDDvsAktmyr$log2FoldChange, rep(102, length(maxedout_MekDDvsAktmyr$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot_MekDDvsAktmyr$log2FoldChange, -log10(significant_points_to_plot_MekDDvsAktmyr$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")
# Add lines
abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")

#T53D4 vs Aktmyr
res_volcano_plot_T53D4vsAktmyr <- lfcShrink(dds,coef="sample_T53D4_vs_Aktmyr", type='apeglm') # Calculate padj without a lower lfc limit
significantLFC_volcano_plot_T53D4vsAktmyr <- subset(resLFC_T53D4_vs_Aktmyr, padj < 0.05) # Identify significantly changing genes
significant_points_to_plot_T53D4vsAktmyr <-res_volcano_plot_T53D4vsAktmyr [which(rownames(res_volcano_plot_T53D4vsAktmyr ) %in% rownames(significantLFC_volcano_plot_T53D4vsAktmyr)),] 
# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout_T53D4vsAktmyr <- subset(res_volcano_plot_T53D4vsAktmyr , padj < 10e-100)

#Draw volcano plot for T53D4 vs Aktmyr
par(mfrow=c(1,1)) # one plot only 
plot(res_volcano_plot_T53D4vsAktmyr $log2FoldChange, -log10(res_volcano_plot_T53D4vsAktmyr $padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-6,6), col = "#00000030")
# Add points
points(maxedout_T53D4vsAktmyr$log2FoldChange, rep(102, length(maxedout_T53D4vsAktmyr$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot_T53D4vsAktmyr$log2FoldChange, -log10(significant_points_to_plot_T53D4vsAktmyr$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")
# Add lines
abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")

#Water vs Aktmyr
res_volcano_plot_WatervsAktmyr <- lfcShrink(dds,coef="sample_negclt_vs_Aktmyr", type='apeglm') # Calculate padj without a lower lfc limit
significantLFC_volcano_plot_WatervsAktmyr <- subset(resLFC_Water_vs_Aktmyr, padj < 0.05) # Identify significantly changing genes
significant_points_to_plot_WatervsAktmyr <-res_volcano_plot_WatervsAktmyr [which(rownames(res_volcano_plot_WatervsAktmyr ) %in% rownames(significantLFC_volcano_plot_T53D4vsAktmyr)),] 
# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout_WatervsAktmyr <- subset(res_volcano_plot_WatervsAktmyr , padj < 10e-100)

#Draw volcano plot for Water vs Aktmyr
par(mfrow=c(1,1)) # one plot only 
plot(res_volcano_plot_WatervsAktmyr $log2FoldChange, -log10(res_volcano_plot_WatervsAktmyr $padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-6,6), col = "#00000030")
# Add points
points(maxedout_WatervsAktmyr$log2FoldChange, rep(102, length(maxedout_WatervsAktmyr$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot_WatervsAktmyr$log2FoldChange, -log10(significant_points_to_plot_WatervsAktmyr$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")
# Add lines
abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")




############## GETTING LISTS OF SIGNIFICANTLY CHANGING GENES #####################
# We will use the set of significantly changing genes from our variance shrunken analysis (resLFC):


# Select the significant subset of genes in ARPE19 vs Aktmyr
Up_in_ARPE19vsAktmyr <- subset(resLFC_ARPE19_vs_Aktmyr, padj < 0.05 & log2FoldChange > 0.5)
Up_in_ARPE19vsAktmyr  <- Up_in_ARPE19vsAktmyr [order(Up_in_ARPE19vsAktmyr$padj),] #order them
head(Up_in_ARPE19vsAktmyr) # Check them
dim(Up_in_ARPE19vsAktmyr)
# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_ARPE19vsAktmyr<- subset(resLFC_ARPE19_vs_Aktmyr, padj < 0.05 & log2FoldChange < -0.5)
Down_in_ARPE19vsAktmyr <- Down_in_ARPE19vsAktmyr[order(Down_in_ARPE19vsAktmyr$padj),]
head(Down_in_ARPE19vsAktmyr)
dim(Down_in_ARPE19vsAktmyr)

# Save these lists to output files:
# write(rownames(Up_in_Aktmyr), file = "../03_output/Genes Up in Aktmyr.txt", sep = "\n")
# write(rownames(Down_in_Aktmyr), file = "../03_output/Genes Down Aktmyr.txt", sep = "\n")
# write(rownames(resLFC_ARPE19_vs_Aktmyr), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.

# Select the significant subset of genes MekDD vs Aktmyr
Up_in_MekDDvsAktmyr <- subset(resLFC_MekDD_vs_Aktmyr, padj < 0.05 & log2FoldChange > 0.5)
Up_in_MekDDvsAktmyr  <- Up_in_Aktmyr [order(Up_in_Aktmyr$padj),] #order them
head(Up_in_MekDDvsAktmyr) # Check them
dim(Up_in_MekDDvsAktmyr)
# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_MekDDvsAktmyr<- subset(resLFC_MekDD_vs_Aktmyr, padj < 0.05 & log2FoldChange < -0.5)
Down_in_MekDDvsAktmyr <- Down_in_MekDDvsAktmyr[order(Down_in_MekDDvsAktmyr$padj),]
head(Down_in_MekDDvsAktmyr)
dim(Down_in_MekDDvsAktmyr)

# Save these lists to output files:
# write(rownames(Up_in_Aktmyr), file = "../03_output/Genes Up in Aktmyr.txt", sep = "\n")
# write(rownames(Down_in_Aktmyr), file = "../03_output/Genes Down Aktmyr.txt", sep = "\n")
# write(rownames(resLFC_ARPE19_vs_Aktmyr), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.

# Select the significant subset of genes in RasV12 vs Aktmyr
Up_in_RasV12vsAktmyr <- subset(resLFC_RasV12_vs_Aktmyr, padj < 0.05 & log2FoldChange > 0.5)
Up_in_RasV12vsAktmyr  <- Up_in_RasV12vsAktmyr [order(Up_in_RasV12vsAktmyr$padj),] #order them
head(Up_in_RasV12vsAktmyr) # Check them
dim(Up_in_RasV12vsAktmyr)
# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_RasV12vsAktmyr<- subset(resLFC_RasV12_vs_Aktmyr, padj < 0.05 & log2FoldChange < -0.5)
Down_in_RasV12vsAktmyr <- Down_in_RasV12vsAktmyr[order(Down_in_RasV12vsAktmyr$padj),]
head(Down_in_RasV12vsAktmyr)
dim(Down_in_RasV12vsAktmyr)

# Save these lists to output files:
# write(rownames(Up_in_Aktmyr), file = "../03_output/Genes Up in Aktmyr.txt", sep = "\n")
# write(rownames(Down_in_Aktmyr), file = "../03_output/Genes Down Aktmyr.txt", sep = "\n")
# write(rownames(resLFC_ARPE19_vs_Aktmyr), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.

# Select the significant subset of genes T53D4 vs Aktmyr
Up_in_T53D4vsAktmyr <- subset(resLFC_T53D4_vs_Aktmyr, padj < 0.05 & log2FoldChange > 0.5)
Up_in_T53D4vsAktmyr  <- Up_in_T53D4vsAktmyr [order(Up_in_T53D4vsAktmyr$padj),] #order them
head(Up_in_T53D4sAktmyr) # Check them
dim(Up_in_T53D4vsAktmyr)
# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_T53D4vsAktmyr<- subset(resLFC_T53D4_vs_Aktmyr, padj < 0.05 & log2FoldChange < -0.5)
Down_in_T53D4vsAktmyr <- Down_in_T53D4vsAktmyr[order(Down_in_T53D4vsAktmyr$padj),]
head(Down_in_T53D4vsAktmyr)
dim(Down_in_T53D4vsAktmyr)

# Save these lists to output files:
# write(rownames(Up_in_Aktmyr), file = "../03_output/Genes Up in Aktmyr.txt", sep = "\n")
# write(rownames(Down_in_Aktmyr), file = "../03_output/Genes Down Aktmyr.txt", sep = "\n")
# write(rownames(resLFC_ARPE19_vs_Aktmyr), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.

# Select the significant subset of genes Water vs Aktmyr
Up_in_WatervsAktmyr <- subset(resLFC_Water_vs_Aktmyr, padj < 0.05 & log2FoldChange > 0.5)
Up_in_WatervsAktmyr  <- Up_in_WatervsAktmyr [order(Up_in_WatervsAktmyr$padj),] #order them
head(Up_in_WatervsAktmyr) # Check them
dim(Up_in_WatervsAktmyr)
# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_WatervsAktmyr<- subset(resLFC_Water_vs_Aktmyr, padj < 0.05 & log2FoldChange < -0.5)
Down_in_WatervsAktmyr <- Down_in_WatervsAktmyr[order(Down_in_WatervsAktmyr$padj),]
head(Down_in_WatervsAktmyr)
dim(Down_in_WatervsAktmyr)


# Save these lists to output files:
# write(rownames(Up_in_Aktmyr), file = "../03_output/Genes Up in Aktmyr.txt", sep = "\n")
# write(rownames(Down_in_Aktmyr), file = "../03_output/Genes Down Aktmyr.txt", sep = "\n")
# write(rownames(resLFC_ARPE19_vs_Aktmyr), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.
############## MAKE HIERARCHICALLY CLUSTERED HEATMAPS OF ALL CHANGING GENES #####################

############## MAKE HIERARCHICALLY CLUSTERED HEATMAPS OF ALL CHANGING GENES #####################

# Let's loosen our restrictions on significance to all genes with any log fold change and adjusted p-values less than 0.1 (both are default)

# Get ARPE19 vs RasV12 differentially expressed genes:
res_ARPE19_vs_Aktmyr <- results(dds,
                                lfc = 1.2,
                                contrast=c("sample","ARPE19","Aktmyr"))

#Subset each results table for just the differentially expressed genes:
ARPE19vsAktmyr<- subset(res_ARPE19_vs_Aktmyr , padj < 0.05)
dim(subset(res_ARPE19_vs_Aktmyr , padj < 0.05))

res_MekDD_vs_Aktmyr <- results(dds,
                                lfc = 1.2,
                                contrast=c("sample","MekDD","Aktmyr"))

MekDDvsAktmyr<- subset(res_MekDD_vs_Aktmyr , padj < 0.05)
dim(subset(res_MekDD_vs_Aktmyr , padj < 0.05))

res_RasV12_vs_Aktmyr <- results(dds,
                               lfc = 1.2,
                               contrast=c("sample","RasV12","Aktmyr"))

RasV12vsAktmyr <- subset(res_RasV12_vs_Aktmyr , padj < 0.05)
dim(subset(res_RasV12_vs_Aktmyr , padj < 0.05))

res_T53D4_vs_Aktmyr <- results(dds,
                                lfc = 1.2,
                                contrast=c("sample","T53D4","Aktmyr"))

T53D4vsAktmyr<- subset(res_T53D4_vs_Aktmyr , padj < 0.01)
dim(subset(res_T53D4_vs_Aktmyr , padj < 0.05))


# res_negclt_vs_Aktmyr <- results(dds,
#                                lfc = 1.2,
#                                contrast=c("sample","negclt","Aktmyr"))
# 
# 
# negcltvsAktmyr<- subset(res_negclt_vs_Aktmyr , padj < 0.05)
# dim(subset(res_negclt_vs_Aktmyr , padj < 0.05))

#Determine how many genes were captured and merge them:
changing_genes<- rbind(ARPE19vsAktmyr, MekDDvsAktmyr, RasV12vsAktmyr, T53D4vsAktmyr)  
dim(changing_genes)
length(unique(rownames(changing_genes)))

# Get all r-stabilized log values of all the data:
rdl_all <- assay(rlog(dds, blind=FALSE))

#Subset just the changing genes:
changing_lrt_rdl <- subset(rdl_all, rownames(rdl_all) %in% rownames(changing_genes))
dim(changing_lrt_rdl)
head(changing_lrt_rdl)

# Make sure it is in matrix form:
class(changing_lrt_rdl)

# Draw a heat map
# Scale by row, listed below as scale = "row". This is really important. It sets the mean of every row to 0 and the standard deviation of every row to 1:

#dcols= dist(t(p), method = "minkowski")
p <- pheatmap(changing_lrt_rdl, 
              scale="row", 
              color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(100),
              border_color = TRUE,
              cluster_rows=TRUE, 
              cluster_cols=TRUE, 
              cutree_rows = 3,
              cutree_cols = 5,
              treeheight_row = 100,
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE) #,km= 10 ##kmeans

plot(p$tree_row)
plot(p$tree_col)
plot(p$gtable)

#identify clusters
hc <-p$tree_col
lbl <- cutree(hc, 5) # split gene dendrogram in 5 groups
which(lbl==2) # grab genes of first group
plot(which(lbl==2))

#identify genes in clusters
rownames(changing_lrt_rdl) 
colnames(changing_lrt_rdl) 
rownames(changing_lrt_rdl[p$tree_row[["order"]],])
cutree(p$tree_row,k=3)
sort(cutree(p$tree_row,k=3))
plot(sort(cutree(p$tree_row,k=3)))
cutree(p$tree_col,k=5)
sort(cutree(p$tree_col,k=5))
plot(sort(cutree(p$tree_col,k=5)))



# install.packages("d3heatmap")
library("d3heatmap")
d3heatmap(changing_lrt_rdl,scale = "row", colors= scales::col_quantile("RdYlBu",NULL,8),
          k_row =9,
          k_col = 5, 
          Rowv = T,
          Colv = T,
          symm = T,
          distfun = dist
           )



# p
# help(pheatmap)
# pdf("../03_output/clustered_genes.pdf", width = 10, height = 12)


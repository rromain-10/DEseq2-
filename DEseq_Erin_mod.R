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
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
BiocManager::install("DESeq2")

# Install 'apeglm'
BiocManager::install("apeglm")

# install corrplot:
install.packages("corrplot")

#Install pretty heatmap - pheatmap: https://cran.r-project.org/web/packages/pheatmap/index.html
install.packages("pheatmap")

#Install RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
install.packages("RColorBrewer")
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

# import the counts data
getwd()

#Set this to your working directory:
setwd("/Users/romarioromain/OneDrive - Colostate/RR_ARPE_DELUCA_COLLAB/DEseq/02_scripts/DEseq2")
getwd()
countsData <- read.table(file = "../../01_input/counts.txt", header = FALSE, row.names = 1, skip = 2) 
head(countsData)
dim(countsData)
class(countsData)

# Read in the metadata
metadata1 <- read.table(file = "../../01_input/RR_ARPE_Deluca_Collab_manifest.txt", header = FALSE) # import the data
metadata1
colnames(metadata1) <- c("fasta1", "fasta2", "names1", "rep","sample")
metadata1

# Organize the countsData file.
# Notice that the countsData file doesn't have any column headers:
head(countsData)

# Let's give countsData some columns names. The first names will be... chr', 'start', etc...
# The last names will be names for each sample. We can pull those names from metadata1:
as.vector(metadata1$names1)
as.vector(metadata1$sample)

# Name countsData columns headers:
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata1$names1))
colnames(countsData)
head(countsData)

################### COUNT MATRIX INPUT ###################
# In this section we will prepare our input data for analysis.
# In the instruction for DESeq2, it states: "We read in a count matrix, which we will name cts, and the sample information table, which we will name coldata."

# generate a table called "cts" out of the countsData table.
# Subset the countsData 
head(countsData)
dim(countsData)
head(countsData[,6:16])

# Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:16])
head(cts)

# Next we need to make an information called coltable. We can make this out of the metadata table.
class(metadata1)
# Reorganize the metadata table so the names1 column are now row headers
metadata1
rownames(metadata1)<- metadata1$names1
metadata1

coldata <- metadata1[,c("sample", "rep","names1", "fasta1")]
coldata
coldata$sample <- as.factor(coldata$sample)
coldata$names1 <- as.factor(coldata$names1)
coldata$rep <- as.factor(coldata$rep)
coldata$fasta1 <- as.factor(coldata$fasta1)
rownames(coldata) <- metadata1$names1
coldata$sample
coldata$names1 
colnames(cts)

# One thing we need to explicitly check. The rownames of coldata need to exactly match the colnames of cts.
#Check that the names match  --> Should be TRUE
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ sample + rep)


################# PRE-FILTERING -- FILTER FOR PRESENT GENES: ################# 
# Not necessary, but helps keep things fast. 
# Exclude all samples that have 0 reads:

# keep <- rowSums(counts(dds1)) >= 1
# dds1 <- dds1[keep,]

keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

# How many did we exclude?
dim(dds)
str(dds)
coldata

################### NOTE ON FACTOR LEVELS ###################
# Organize the categories based on what makes sense:
dds2$sample <- factor(dds$sample)
dds2$sample

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
dds2$sizeFactor
head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))

############## PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION #####################

# An input requirement of the lfcShrink function is a coef term. This is pulled from the resultsNames of dds:
# resultsNames(dds1)
resultsNames(dds)


############## DIFFERENTIAL EXPRESSION ANALYSIS #####################

# calculate the statistically significant differences between ARPE19 vs Aktmyr
#ARPE19 vs Aktmyr
res_ARPE19_vs_Aktmyr <- results(dds,
                                lfc = 0.5,
                                contrast=c("sample","ARPE19","Aktmyr"))

summary(res_ARPE19_vs_Aktmyr)

resLFC_ARPE19_vs_Aktmyr <- lfcShrink(dds,
                                     coef="sample_ARPE19_vs_Aktmyr", type='apeglm')

summary(resLFC_ARPE19_vs_Aktmyr)

#RasV12 vs Aktmyr
res_RasV12_vs_Aktmyr <- results(dds,
                                lfc = 0.5,
                                contrast=c("sample","RasV12","Aktmyr"))

summary(res_RasV12_vs_Aktmyr)

resLFC_RasV12_vs_Aktmyr <- lfcShrink(dds,
                                     coef="sample_RasV12_vs_Aktmyr", type='apeglm')

summary(resLFC_RasV12_vs_Aktmyr)



##################  Exploring and exporting results ##################  
#ARPE19 vs Aktmyr
par(mfrow=c(1,1))
plotMA(res_ARPE19_vs_Aktmyr, main="ARPE19 vs neg Aktmyr\nunshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized ARPE19 / Aktmyr)",
       xlab = "means of normalized counts")

plotMA(resLFC_ARPE19_vs_Aktmyr, main="ARPE19 vs Aktmyr\nshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized ARPE19 / Aktmyr)",
       xlab = "means of normalized counts")

#Identify genes on the plot ARPE19 vs Aktmyr
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_ARPE19_vs_Aktmyr$baseMean,resLFC_ARPE19_vs_Aktmyr$log2FoldChange)
#  Step4 -> click here to see what you got!
rownames(resLFC_ARPE19_vs_Aktmyr)[idx]

#RasV12 vs Aktmyr
par(mfrow=c(1,1))
plotMA(res_RasV12_vs_Aktmyr, main="RasV12 vs Aktmyr\nunshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized RasV12 / Aktmyr)",
       xlab = "means of normalized counts")

plotMA(resLFC_RasV12_vs_Aktmyr, main="RasV12 vs Aktmyr\nshrunken", ylim = c(-7,7),
       ylab = "log fold change (ratio of normalized RasV12 / Aktmyr)",
       xlab = "means of normalized counts")

#Identify genes on the plot ARPE19 vs Aktmyr
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_RasV12_vs_Aktmyr$baseMean,resLFC_RasV12_vs_Aktmyr$log2FoldChange)
#  Step4 -> click here to see what you got!
rownames(resLFC_RasV12_vs_Aktmyr)[idx]

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
         col=colors, width = 4, height = 4) # Plot the heatmap

##### VOLCANO PLOTS:###################
# resultsNames(dds1)
resultsNames(dds2)


# Volcano plots are nice ways of displaying the fold change against the p-value.
res_volcano_plot <- lfcShrink(dds2,coef="sample_ARPE19_vs_Aktmyr", type='apeglm') # Calculate padj without a lower lfc limit
significantLFC <- subset(resLFC_ARPE19_vs_Aktmyr, padj < 0.05) # Identify significantly changing genes
significant_points_to_plot <-res_volcano_plot[which(rownames(res_volcano_plot) %in% rownames(significantLFC)),] 
# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout <- subset(res_volcano_plot, padj < 10e-100)

#Draw plot:
par(mfrow=c(1,1)) # one plot only 
# Draw the plot
plot(res_volcano_plot$log2FoldChange, -log10(res_volcano_plot$padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-6,6), col = "#00000030")
# Add points
points(maxedout$log2FoldChange, rep(102, length(maxedout$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot$log2FoldChange, -log10(significant_points_to_plot$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")
# Add lines
abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")

############## GETTING LISTS OF SIGNIFICANTLY CHANGING GENES #####################

# We will use the set of significantly changing genes from our variance shrunken analysis:
#Remember we calculated this above as:
res_ARPE19_vs_RasV12 <- results(dds2,
                                lfc = 0.5,
                                contrast=c("sample","ARPE19","Aktmyr"))

summary(res_ARPE19_vs_RasV12 )
resultsNames(dds2)

resLFC_ARPE19_vs_Aktmyr <- lfcShrink(dds2,
                                     coef="sample_ARPE19_vs_Aktmyr", res =  res_ARPE19_vs_RasV12, type='apeglm' )

summary(resLFC_ARPE19_vs_Aktmyr)

# Check the results table:
head(resLFC_ARPE19_vs_Aktmyr)

# Select the significant subset of genes that are up-regulated in acetic acid
Up_in_Aktmyr <- subset(resLFC_ARPE19_vs_Aktmyr, padj < 0.05 & log2FoldChange > 0.5)
Up_in_Aktmyr  <- Up_in_Aktmyr [order(Up_in_Aktmyr$padj),] #order them
head(Up_in_Aktmyr) # Check them
dim(Up_in_Aktmyr)
# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_Aktmyr<- subset(resLFC_ARPE19_vs_Aktmyr, padj < 0.05 & log2FoldChange < -0.5)
Down_in_Aktmyr <- Down_in_Aktmyr[order(Down_in_Aktmyr$padj),]
head(Down_in_Aktmyr)
dim(Down_in_Aktmyr)

# Save these lists to output files:
# write(rownames(Up_in_Aktmyr), file = "../03_output/Genes Up in Aktmyr.txt", sep = "\n")
# write(rownames(Down_in_Aktmyr), file = "../03_output/Genes Down Aktmyr.txt", sep = "\n")
# write(rownames(resLFC_ARPE19_vs_Aktmyr), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.

############## MAKE HIERARCHICALLY CLUSTERED HEATMAPS OF ALL CHANGING GENES #####################

# Let's loosen our restrictions on significance to all genes with any log fold change and adjusted p-values less than 0.1 (both are default)

# Get ARPE19 vs RasV12 differentially expressed genes:
res_ARPE19_vs_RasV12 <- results(dds2,
                                lfc = 0.5,
                                contrast=c("sample","ARPE19","Aktmyr"))

#Subset each results table for just the differentially expressed genes:
ARPE19vsAktmyr<- subset(res_ARPE19_vs_RasV12 , padj < 0.05)
dim(subset(res_ARPE19_vs_RasV12 , padj < 0.05))

#Determine how many genes were captured and merge them:
changing_genes <- rbind(ARPE19vsAktmyr)
dim(changing_genes)
length(unique(rownames(changing_genes)))

# Get all r-stabilized log values of all the data:
rdl_all <- assay(rlog(dds2, blind=FALSE))

#Subset just the changing genes:
changing_lrt_rdl <- subset(rdl_all, rownames(rdl_all) %in% rownames(changing_genes))
dim(changing_lrt_rdl)
head(changing_lrt_rdl)

# Make sure it is in matrix form:
class(changing_lrt_rdl)

# Draw a heat map
# Scale by row, listed below as scale = "row". This is really important. It sets the mean of every row to 0 and the standard deviation of every row to 1:
p <- pheatmap(changing_lrt_rdl, 
              scale="row", 
              color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=TRUE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE)
# p
# help(pheatmap)
# pdf("../03_output/clustered_genes.pdf", width = 10, height = 12)

###########TOP Go################################
BiocManager::install("topGO")
library(topGO)
data(cts)

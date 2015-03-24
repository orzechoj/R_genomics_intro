## Eaxmples of explorative data analysis

#####
## Download expression data from (Graveley et al. The developmental transcriptome of Drosophila melanogaster, Nature, 2012),
## Supplementary table 9, 2010-09-11375C-Table_S9_revised.xls
library(gdata)
exp.data.raw <- read.xls("2010-09-11375C-Table_S9_revised.xls")
exp.data <- as.matrix(exp.data.raw[,-1*1:4]) ## Expression levels are in column 5 and higher
rownames(exp.data) <- exp.data.raw[,1] ## Gene names are in column 1
colnames(exp.data) <- gsub("X", "", colnames(exp.data)) ## Clean up coumn names
colnames(exp.data) <- gsub("..", ".", colnames(exp.data), fixed=TRUE)


#####
## PCA plot
library(ggplot2)

test.pca <- prcomp(t(log10(exp.data+1)))

# plot how much of the variation the PCs explain
barplot(summary(test.pca)$importance[2,1:10])

annot.col <- c(rep("embryo", 12), rep("larvae", 6), rep("pupae", 6), rep("adult",6))
scores <- data.frame(annot.col, test.pca$x[,1:5])

## Plot PC1 vs PC2
pc <- qplot(x=PC1, y=PC2, data=scores, colour=factor(annot.col), label=rownames(scores))
pc <- pc + geom_text(size=4)
pc


#####
## Clustering
library(gplots)
options("expressions"=15000)

annot.colors <- c(rep("red", 12), rep("green", 6), rep("cyan", 6), rep("purple",6))

heatmap.2(log10(exp.data+1), labRow = NA, trace="none", cexCol=0.7, margins = c(10, 2),	ColSideColors=annot.colors)

# Library
library(lattice)
library(dplyr)
library(reshape2)
library(ggplot2)

Protein.file <- "MS5504_MCap_Protein_Abundances.csv"
Gene.file <- "Protein_matches_in_Mcap.mRNA.fa.salmon.numreads.matrix.csv"
  
# Load matrices
protein.data <- read.table(Protein.file, header=T, check.names = FALSE, sep=',')
gene.data <- read.table(Gene.file, header=T, check.names = FALSE, sep=',')

# Transpose data sets
protein.data.melt <- melt(protein.data)
colnames(protein.data.melt)<-c("Accession","TP_Condition","Protein_count")
gene.data.melt <- melt(gene.data)
colnames(gene.data.melt)<-c("Accession","TP_Condition","Gene_count")

# Merge 2 datasets
data <- merge(protein.data.melt,gene.data.melt,by=c("Accession","TP_Condition"))

# Remove duplicated rows
data.filtered <- data[!duplicated(data$Protein_count), ] 

# Normalize counts
mode(data.filtered$Protein_count)<-'numeric'
mode(data.filtered$Gene_count) <-'numeric'
data.filtered$Norm_gene_count <- (data.filtered$Gene_count*0.0100598304933868)
data.filtered$log.protein.count = log2(data.filtered$Protein_count+1)
data.filtered$log.gene.count = log2(data.filtered$Norm_gene_count+1)

## Plot distribution of 2 datasets based on shared sequences for multiple TPs
a <- xyplot(data.filtered$log.protein.count ~ data.filtered$log.gene.count | data.filtered$TP_Condition , 
       aspect = 1:1,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         fm <- lm(y ~ x - 0)
         panel.lines(x, fitted(fm), col.line = "red", ...)
         panel.text(7.8,1, labels=sprintf("y = %.3f x + %.3f\nR² = %.2f\nnPoints = %.0f", 
                                          coef(fm)[2], coef(fm)[1], summary(fm)$r.squared, length(x)),cex=.7)
       },
       scales = "free",
       par.settings=simpleTheme(col="blue"),
       pch=10, 
       cex=.5,
       xlim = c(0, 10),
       ylim = c(0, 15),
       ylab = list(label=expression("log2 Protein Expression"), fontsize=15),
       xlab = list(label="log2 Gene Expression",fontsize=15),
       layout=c(3,3),
       index.cond=list(c(1,3,5,2,4,6,7))
)
trellis.device(device="pdf", file="Prot_v_Gene_Expression_wField.pdf", paper="a4", height=8.3, width=11.7)
print(a)
dev.off()

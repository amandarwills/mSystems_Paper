# Library
library(devtools)
library(magrittr)
library(lattice)
library(dplyr)
library(reshape2)
library(ggplot2)

Protein.file <- "MS5504_MCap_Protein_Abundances.csv"
Gene.file <- "Mcap.mRNA.fa.salmon.tpm.matrix"

# Load matrices
protein.data <- read.table(Protein.file, header=T, check.names = FALSE, sep=',')
gene.data <- read.table(Gene.file, header=T, check.names = FALSE, sep='\t')

# Get averages of proteins within their treatment and TP
avg.protein <- protein.data %>% 
                        dplyr::mutate(PA = rowMeans(select(protein.data,
                                                           starts_with('T1A')), na.rm = TRUE),
                                      PH = rowMeans(select(protein.data,
                                                           starts_with('T1H')), na.rm = TRUE)
  )

# Select columns to keep in avg.proteins (Accession, 2 avg columns)
avg.protein <- avg.protein[, c(1, 16:17)] 

# Get averages of genes within their treatment and TP
avg.genes <- gene.data %>% 
                      dplyr::mutate(GA = rowMeans(select(gene.data,
                                                         starts_with('TP1_A')), na.rm = TRUE),
                                    GH = rowMeans(select(gene.data,
                                                         starts_with('TP1_H')), na.rm = TRUE)
  )

# Select columns to keep in avg.genes (Accession, 2 avg columns)
avg.genes <- avg.genes[, c(1, 23:24)] 

# Merge proteins with genes
protein.gene.joined <- avg.protein %>% left_join(avg.genes, 
                                             by = c('Accession' = 'Name'))

# Log2 transform protein and gene counts
protein.gene.joined$log.PA = log2(protein.gene.joined$PA+1)
protein.gene.joined$log.PH = log2(protein.gene.joined$PH+1)
protein.gene.joined$log.GA = log2(protein.gene.joined$GA+1)
protein.gene.joined$log.GH = log2(protein.gene.joined$GH+1)

# Remove NAs from each data frame 
na.strings=c(""," ","NA")
protein.gene.joined <- na.omit(protein.gene.joined, c("logGA","logGH"))

# Create separate data frames for thermal stress and control 
PA.GA.data <- protein.gene.joined[,c(1, 6, 8)]
PH.GH.data <- protein.gene.joined[,c(1,7,9)]
rownames(protein.gene.joined)=protein.gene.joined$Accession
protein.gene.joined$Accession <- NULL
protein.gene.joined.matrix <- as.matrix(protein.gene.joined)

# Create matrices
rownames(PA.GA.data)=PA.GA.data$Accession
PA.GA.data$Accession <- NULL
PA.GA.data.matrix <- as.matrix(PA.GA.data)
rownames(PH.GH.data)=PH.GH.data$Accession
PH.GH.data$Accession <- NULL
PH.GH.data.matrix <- as.matrix(PH.GH.data)

# Heatmap 
heatmap(PA.GA.data.matrix,
        labRow = F, labCol = F,
        Colv = NA, Rowv = NA, 
        scale="column",
        margins = c(1, 10),
        main = "T1A Protein-Gene Expression") 



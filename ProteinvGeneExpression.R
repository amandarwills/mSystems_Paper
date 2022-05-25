# Library
library(devtools)
library(magrittr)
library(lattice)
library(dplyr)
library(reshape2)
library(ggplot2)
library(plotrix)
library(ggpubr)

Protein.file <- "MC_DEP_T5.txt"
Gene.file <- "Mcap.mRNA.fa.salmon.tpm.matrix"
Sequence.ID.file <- "MC_FC1.2_DEP_T5_IDs.csv"
output <- "MC_T5_protvgene_correlation_results.csv"
Prot_A <- "T5A"
Prot_H <- "T5H"
Gene_A <- "TP5_A"
Gene_H <- "TP5_H"
Ambient_graph <- "TP5 Ambient"
HiT_graph <- "TP5 Thermal Stress"
fc.1 <- 1.8
fc.2 <- -1.8

# Load matrices
protein.data <- read.table(Protein.file, header=T, check.names = FALSE, sep='\t')
gene.data <- read.table(Gene.file, header=T, check.names = FALSE, sep='\t')
Sequence.IDs <- read.table(Sequence.ID.file, fill=TRUE, header=T, quote="", sep=',', encoding="UTF-8")

### Expression levels of highest diff exp proteins ###

# Filter the data and keep only DEPs with FC
sigProt <- protein.data %>% 
                          arrange(desc(abs(FC)))
top.DEPS <- sigProt %>% 
                      dplyr::filter(FC <= fc.2 | FC >= fc.1)

# Get averages of proteins within their treatment and TP
avg.top.DEPS <- top.DEPS %>% 
                         dplyr::mutate(PA = rowMeans(select(top.DEPS,
                                                           starts_with(Prot_A)), na.rm = TRUE),
                                       PH = rowMeans(select(top.DEPS,
                                                           starts_with(Prot_H)), na.rm = TRUE)
                                       )

# Select columns to keep in avg.top.DEPS (Accession, FC, 2 avg columns)
avg.top.DEPS <- avg.top.DEPS[, c(1, 2, 7, 15:16)] 

# Get averages of genes within their treatment and TP
avg.Genes <- gene.data %>% 
                          dplyr::mutate(GA = rowMeans(select(gene.data,
                                                               starts_with(Gene_A)), na.rm = TRUE),
                                        GH = rowMeans(select(gene.data,
                                                               starts_with(Gene_H)), na.rm = TRUE)
                                        )

# Select columns to keep in avg.genes (Accession, FC, 2 avg columns)
avg.Genes <- avg.Genes[, c(1, 23:24)] 

# Merge top DEPS with genes
top.DEPS.Genes <- avg.top.DEPS %>% left_join(avg.Genes, 
                                                  by = c('Accession' = 'Name'))

# Add sequence IDs 
top.DEPS.Genes <- top.DEPS.Genes %>% left_join(Sequence.IDs)
top.DEPS.Genes$ID[is.na(top.DEPS.Genes$ID)] <- "Uncharacterized"
top.DEPS.Genes <- top.DEPS.Genes[, c(1, 2, 8, 3:7)]

# Log2 transform protein and gene counts
top.DEPS.Genes$log.PA = log2(top.DEPS.Genes$PA+1)
top.DEPS.Genes$log.GA = log2(top.DEPS.Genes$GA+1)
top.DEPS.Genes$log.PH = log2(top.DEPS.Genes$PH+1)
top.DEPS.Genes$log.GH = log2(top.DEPS.Genes$GH+1)

PA.ordered.top.DEPS.Genes <- top.DEPS.Genes %>% 
                                             arrange((abs(log.PA)))
PH.ordered.top.DEPS.Genes <- top.DEPS.Genes %>% 
                                            arrange((abs(log.PH)))

# Write table to csv
write.csv(top.DEPS.Genes, output, row.names = FALSE)

# Plot back to back histograms
par(mar=pyramid.plot(PA.ordered.top.DEPS.Genes$log.PA, PA.ordered.top.DEPS.Genes$log.GA,
                    top.labels=c("Protein counts","Protien coding gene ID","Gene count"),
                    unit = "log2 counts",
                    labelcex = 0.65,
                    space = 0.3,
                    xlim = c(12,12),
                    laxlab = seq(from = 0, to = 12, by = 0.5),
                    raxlab=seq(from = 0, to = 12, by = 0.5),
                    ppmar=c(4,2,4,2),
                    labels = PA.ordered.top.DEPS.Genes$Protein.ID,
                    main=Ambient_graph,
                    lxcol="blue",rxcol= "green",
                    gap=0.8, show.values=F)
    )

par(mar=pyramid.plot(PH.ordered.top.DEPS.Genes$log.PH, PH.ordered.top.DEPS.Genes$log.GH,
                     top.labels=c("Protein counts","Protien coding gene ID","Gene count"),
                     unit = "log2 counts",
                     labelcex = 0.65,
                     space = 0.3,
                     xlim = c(12,12),
                     laxlab = seq(from = 0, to = 12, by = 0.5),
                     raxlab=seq(from = 0, to = 12, by = 0.5),
                     ppmar=c(4,2,4,2),
                     labels = PH.ordered.top.DEPS.Genes$Protein.ID,
                     main=HiT_graph,
                     lxcol="blue",rxcol= "green",
                     gap=0.8, show.values=F)
    )






### All proteins ###

All.Protein.file <- "MS5504_MCap_Protein_Abundances.csv"
All.Gene.file <- "Mcap.mRNA.fa.salmon.tpm.matrix"

# Load matrices
all.protein.data <- read.table(All.Protein.file, header=T, check.names = FALSE, sep=',')
all.gene.data <- read.table(All.Gene.file, header=T, check.names = FALSE, sep='\t')

# Get averages of proteins within their treatment and TP
avg.protein.data <- all.protein.data %>% 
                                    dplyr::mutate(P_TP1A = rowMeans(select(all.protein.data,
                                                                          starts_with("TP1_A")), na.rm = TRUE),
                                                  P_TP1H = rowMeans(select(all.protein.data,
                                                                          starts_with("TP1_H")), na.rm = TRUE),
                                                  P_TP3A = rowMeans(select(all.protein.data,
                                                                           starts_with("TP3_A")), na.rm = TRUE),
                                                  P_TP3H = rowMeans(select(all.protein.data,
                                                                           starts_with("TP3_H")), na.rm = TRUE),
                                                  P_TP5A = rowMeans(select(all.protein.data,
                                                                           starts_with("TP5_A")), na.rm = TRUE),
                                                  P_TP5H = rowMeans(select(all.protein.data,
                                                                           starts_with("TP5_H")), na.rm = TRUE),
                                                  P_Field = rowMeans(select(all.protein.data,
                                                                           starts_with("TP0")), na.rm = TRUE)
                                                  )
avg.protein.data <- avg.protein.data[, c(1, 16:22)]

# Get averages of proteins within their treatment and TP
avg.gene.data <- all.gene.data %>% 
                                    dplyr::mutate(T_TP1A = rowMeans(select(all.gene.data,
                                                                           starts_with("TP1_A")), na.rm = TRUE),
                                                  T_TP1H = rowMeans(select(all.gene.data,
                                                                           starts_with("TP1_H")), na.rm = TRUE),
                                                  T_TP3A = rowMeans(select(all.gene.data,
                                                                           starts_with("TP3_A")), na.rm = TRUE),
                                                  T_TP3H = rowMeans(select(all.gene.data,
                                                                           starts_with("TP3_H")), na.rm = TRUE),
                                                  T_TP5A = rowMeans(select(all.gene.data,
                                                                           starts_with("TP5_A")), na.rm = TRUE),
                                                  T_TP5H = rowMeans(select(all.gene.data,
                                                                           starts_with("TP5_H")), na.rm = TRUE),
                                                  T_Field = rowMeans(select(all.gene.data,
                                                                            starts_with("TP0")), na.rm = TRUE)
                                    )
avg.gene.data <- avg.gene.data[, c(1, 23:29)]

# Merge top DEPS with genes
Prots.Genes <- avg.protein.data %>% left_join(avg.gene.data, 
                                             by = c('Accession' = 'Name'))

# Transform df into matrix
Prots.Genes[is.na(Prots.Genes)] <- 0
rownames(Prots.Genes)=Prots.Genes$Accession
Prots.Genes$Accession <- NULL
Prots.Genes.matrix <- as.matrix(Prots.Genes)

# Log2 transform protein and gene counts and ensure matrix contains numeric data 
#Prots.Genes.matrix <- log2(Prots.Genes.matrix+1) 
#mode(Prots.Genes.matrix)

# Get differences between protein-gene expression levels
Prots.Genes.df <- as.data.frame(Prots.Genes.matrix)

Prots.Genes.df$TP1A.diff <- Prots.Genes.df$P_TP1A - Prots.Genes.df$T_TP1A
Prots.Genes.df$TP1H.diff <- Prots.Genes.df$P_TP1H - Prots.Genes.df$T_TP1H

Prots.Genes.df$TP3A.diff <- Prots.Genes.df$P_TP3A - Prots.Genes.df$T_TP3A
Prots.Genes.df$TP3H.diff <- Prots.Genes.df$P_TP3H - Prots.Genes.df$T_TP3H

Prots.Genes.df$TP5A.diff <- Prots.Genes.df$P_TP5A - Prots.Genes.df$T_TP5A
Prots.Genes.df$TP5H.diff <- Prots.Genes.df$P_TP5H - Prots.Genes.df$T_TP5H

Prots.Genes.df$Field.diff <- Prots.Genes.df$P_Field - Prots.Genes.df$T_Field

# Write table to csv
write.csv(Prots.Genes.df, "MC_full_protvgene_correlation_results.csv", row.names = TRUE, col.names = TRUE)






### FC values of proteins and genes ###

FC.Protein.file <- "MC_DEP_T1.txt"
FC.Gene.file <- "Mcap.mRNA.fa.salmon.numreads.matrix_DiffExprResults.txt"

output <- "MC_T1_FC_protvgene_correlation_results.csv"
DEG.controlSample.name <- "TP1_Amb"
Plot.title <- "TP1 HiT vs Ambient FC correlation between DEPs and DEGs"

# Load matrices
FC.protein.data <- read.table(FC.Protein.file, header=T, check.names = FALSE, sep='\t')
FC.gene.data <- read.table(FC.Gene.file, header=T, check.names = FALSE, sep='\t')

# Pull correct samples from DEG data
DEGs.by.TP <- FC.gene.data %>% 
                              dplyr::filter(controlSample == DEG.controlSample.name)

# Format DEG data to have Accession and FC
DEGs <- DEGs.by.TP[, c(1, 3)]

# Filter DEPs FC
top.DEPS <- FC.protein.data %>% 
                           dplyr::filter(FC <= -0.1 | FC >= 0.1)

# Format DEP data to have Accession and FC
DEPs <- top.DEPS[, c(2, 7)]

# Merge DEGs with DEPs
merged.DEPs.DEGs <- DEPs %>% left_join(DEGs, 
                                           by = c('Accession' = 'seqName'))

# Remove NAs
na.strings=c(""," ","NA")
merged.DEPs.DEGs <- na.omit(merged.DEPs.DEGs, c("log2FoldChange"))

# Order DEPS by degree of FC
merged.DEPs.DEGs <- merged.DEPs.DEGs %>% 
                              arrange(desc(FC))

# Rename columns and melt data frame for plotting
merged.DEPs.DEGs <- merged.DEPs.DEGs %>%
                                      dplyr::rename(Accession = 1,
                                                    Protein = 2,
                                                    Transcript = 3)

melt.merged.DEPs.DEGs <- melt(merged.DEPs.DEGs)

# Format data frame for plotting
final <- melt.merged.DEPs.DEGs %>%
                              dplyr::rename(Accession = 1,
                                            Data.type = 2,
                                            FC = 3)

# Plot stacked bar chart
a <- ggplot(final, aes(fill = Data.type, y = Accession, x = FC)) + 
                          geom_bar(position = 'stack', stat = 'identity') +
                          theme_bw() +
                          labs(x = 'log2 fold change', y = 'Gene') +
                          scale_fill_manual('Position', values = c('red3', 'steelblue'),
                                            labels = c("Protein", "Transcript"))  +
                          scale_x_continuous(breaks = seq(-10, 10, by = 1)) +
                          theme(axis.text.y = element_blank(),
                                axis.ticks.y = element_blank()) +
                          guides(fill = guide_legend(title = NULL)) 

# Combine plots
figure <- ggarrange(a, b, c,
                    labels = c("TP1", "TP3", "TP5"),
                    ncol = 3, nrow = 1,
                    widths = c(1, 1.2, 1.2),
                    font.label = list(size = 12,face = "bold", family = NULL),
                    vjust = 0.6,
                    hjust = 0.1,
                    common.legend = TRUE, legend = "right") +
                    theme(plot.margin = margin(0.3,0.1,0.1,0.5, "cm"))

figure

ggsave(filename="Plots.pdf", figure, width=11, height=8.5)
# End

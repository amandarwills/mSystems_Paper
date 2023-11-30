####################################################################################################
#################################### SCRIPTS FOR MSYSTEMS PAPER ####################################
####################################################################################################

##################### PLOT ALL PROTEINS AND GENES BY ABUNDANCES VIA SCATTERPLOT ####################

# Clear workspace and restart R
rm(list=ls())  
.rs.restartR() 

# Library
library(stringr)
library(data.table)
library(lattice)
library(dplyr)
library(reshape2)
library(ggplot2)

# Set working directory
workingDir <- "/Users/amandawilliams/Desktop/Fix_Proteomics/Input_Files"
setwd(workingDir)

# Set data files names
Protein.file <- "MS5504_MCap_Protein_Abuandance_Report.csv"
Gene.file <- "Mcap.mRNA.fa.salmon.tpm.matrix"

# Load files
protein.data <- read.table(Protein.file, header=T, row.names=1, check.names = FALSE, sep=',')
gene.data <- read.table(Gene.file, header=T, row.names=1, check.names = FALSE, sep='\t')

# Replace NAs with 0s for proteomic data
protein.data[is.na(protein.data)] <- 0
sum(is.na(protein.data))
which(is.na(protein.data))

# Filter out gene data to only include samples/genes present in proteomics data
gene.data.samples.filtered <- gene.data[ ,colnames(protein.data)]
gene.data.filtered <- gene.data.samples.filtered[rownames(protein.data), ]

# Make rownames first column
setDT(gene.data.filtered, keep.rownames = "Accession")
setDT(protein.data, keep.rownames = "Accession")

# Get averages of proteins within their treatment and TP
protein.data <- protein.data %>% 
            dplyr::mutate(T1_Amb = rowMeans(select(protein.data,
                                                   starts_with("MC-289_T1-Amb")), na.rm = TRUE),
                          T1_HiT = rowMeans(select(protein.data,
                                                   starts_with("MC-289_T1-HiT")), na.rm = TRUE),
                          T3_Amb = rowMeans(select(protein.data,
                                                   starts_with("MC-289_T3-Amb")), na.rm = TRUE),
                          T3_HiT = rowMeans(select(protein.data,
                                                   starts_with("MC-289_T3-HiT")), na.rm = TRUE),
                          T5_Amb = rowMeans(select(protein.data,
                                                   starts_with("MC-289_T5-Amb")), na.rm = TRUE),
                          T5_HiT = rowMeans(select(protein.data,
                                                   starts_with("MC-289_T5-HiT")), na.rm = TRUE),
                          Field = rowMeans(select(protein.data,
                                                  starts_with("MC-289_Field")), na.rm = TRUE)
            )

protein.data <- protein.data[, c(1, 16:22)]

# Get averages of proteins within their treatment and TP
gene.data.filtered <- gene.data.filtered %>% 
            dplyr::mutate(T1_Amb = rowMeans(select(gene.data.filtered,
                                                   starts_with("MC-289_T1-Amb")), na.rm = TRUE),
                          T1_HiT = rowMeans(select(gene.data.filtered,
                                                   starts_with("MC-289_T1-HiT")), na.rm = TRUE),
                          T3_Amb = rowMeans(select(gene.data.filtered,
                                                   starts_with("MC-289_T3-Amb")), na.rm = TRUE),
                          T3_HiT = rowMeans(select(gene.data.filtered,
                                                   starts_with("MC-289_T3-HiT")), na.rm = TRUE),
                          T5_Amb = rowMeans(select(gene.data.filtered,
                                                   starts_with("MC-289_T5-Amb")), na.rm = TRUE),
                          T5_HiT = rowMeans(select(gene.data.filtered,
                                                   starts_with("MC-289_T5-HiT")), na.rm = TRUE),
                          Field = rowMeans(select(gene.data.filtered,
                                                  starts_with("MC-289_Field")), na.rm = TRUE)
            )

gene.data.filtered <- gene.data.filtered[, c(1, 16:22)]

# Transpose data sets
protein.data.melt <- melt(protein.data, na.rm = TRUE)
colnames(protein.data.melt) <- c("Accession","TP_Condition","Protein_count")
gene.data.melt <- melt(gene.data.filtered, na.rm = TRUE)
colnames(gene.data.melt) <- c("Accession","TP_Condition","Gene_count")

# Merge 2 datasets
data <- merge(protein.data.melt, gene.data.melt, by = c("Accession", "TP_Condition"))
length(unique(data$Accession))
sum(is.na(data))
which(is.na(data))

# Remove duplicated rows
#data.filtered <- data[!duplicated(data$Accession & data$TP_Condition), ] 

# Simplify names to remove replicates numbers/plugIDs
data$TP_Condition <- sub("MC-289_", "", data$TP_Condition, ignore.case = TRUE)

# data$TP_Condition <- ifelse(str_sub(data$TP_Condition, 1, 1) == "T", 
#                             str_sub(data$TP_Condition, 1, 6), 
#                             data$TP_Condition) # Remove plugIDs
# 
# data$TP_Condition <- ifelse(str_sub(data$TP_Condition, 1, 1) == "F", 
#                             str_sub(data$TP_Condition, 1, 5), 
#                             data$TP_Condition) # Remove field replicate numbers

# Log2 normalize counts
mode(data$Protein_count) <- 'numeric'
mode(data$Gene_count) <- 'numeric'
data$log.protein.count <- log2(data$Protein_count+1)
data$log.gene.count <- log2(data$Gene_count+1)

# Final data check
sum(is.na(data))
which(is.na(data))
length(unique(data$Accession))

# Plot distribution of 2 datasets based on shared sequences for multiple TPs
a <- xyplot(data$log.protein.count ~ data$log.gene.count | data$TP_Condition, 
            aspect = 1:1,
            panel = function(x, y, ...) {
                        panel.xyplot(x, y, ...)
                        fm = lm(y ~ x - 0)
                        panel.lines(x, fitted(fm), col.line = "red", ...)
                        panel.text(11.8,1, labels = sprintf("y = %.3f x + %.3f\nR² = %.2f\nnPoints = %.0f", 
                                                            coef(fm)[2], coef(fm)[1], summary(fm)$r.squared, length(x)), cex=.7)
            },
            scales = "free",
            par.settings = simpleTheme(col = "blue"),
            pch = 10, 
            cex = .5,
            xlim = c(0, 15),
            ylim = c(0, 15),
            ylab = list(label = "Protein Expression [log2(normalized abundance+1)]", fontsize = 15),
            xlab = list(label = "Transcript Expression [log2(TPM+1)]", fontsize = 15),
            layout = c(3, 3),
            index.cond = list(c(2, 4, 6, 3, 5, 7, 1))
)
trellis.device(device="pdf", file="Prot_v_Gene_Expression_wField.pdf", paper="a4", height=8.3, width=11.7)
print(a)
dev.off()
a

write.csv(data, "Prot_Gene_Expression_Corr.csv", row.names = FALSE, col.names = TRUE)




##################### PLOT PROTEINS AND GENES BY FC VIA SCATTERPLOT/BAR CHART ######################

# Clear workspace and restart R
rm(list=ls())  
.rs.restartR() 

# Library
library(devtools)
library(magrittr)
library(lattice)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggpmisc)
library(plotrix)
library(ggpubr)

# Display the current working directory
getwd()

# If necessary, change the path below to the directory where the data files are stored 
workingDir <- "/Users/amandawilliams/Desktop/Fix_Proteomics/Input_Files"
setwd(workingDir)

# Set file names, TP, and file output names
FC.Protein.file <- "MC_DEP_T5_HiTvsAmb.txt"
FC.Gene.file <- "Mcap.mRNA.fa.salmon.numreads.matrix_DiffExprResults.txt"

output <- "MC_T5_FC_protvgene_correlation_results.csv"
DEG.controlSample.name <- "TP5_Amb"
Plot.title <- "T5 HiT vs Ambient FC correlation between DEPs and DEGs"

# Load files
FC.protein.data <- read.table(FC.Protein.file, header=T, check.names = FALSE, sep='\t')
FC.gene.data <- read.table(FC.Gene.file, header=T, check.names = FALSE, sep='\t')

# Pull correct samples from DEG data
DEGs.by.TP <- FC.gene.data %>% 
                        dplyr::filter(controlSample == DEG.controlSample.name)

# Format DEG data to have Accession and FC
DEGs <- DEGs.by.TP[, c(1, 3)]

# Format DEP data to have Accession and FC
DEPs <- FC.protein.data[, c(2, 17)]

# Merge DEGs with DEPs
merged.DEPs.DEGs <- DEPs %>% left_join(DEGs, 
                                       by = c('Accession' = 'seqName'))

# Remove NAs
na.strings <- c(""," ","NA")
merged.DEPs.DEGs <- na.omit(merged.DEPs.DEGs, c("log2FoldChange"))

# Order DEPS by degree of FC
merged.DEPs.DEGs <- merged.DEPs.DEGs %>% 
                                    arrange(desc(FC))

# Rename columns
merged.DEPs.DEGs <- merged.DEPs.DEGs %>%
                                    dplyr::rename(Protein = 2,
                                                  Transcript = 3)


# Designate whether or not each row is DEG/DEP
merged.DEPs.DEGs <- merged.DEPs.DEGs %>% 
            mutate(Regulation = case_when(Transcript >= 1 & Protein >= 1  ~ "DEG and DEP",
                                          Transcript <= -1 & Protein <= -1  ~ "DEG and DEP",
                                          Transcript >= 1 & Protein <= -1   ~ "DEG and DEP",
                                          Transcript <= -1 & Protein >= 1  ~ "DEG and DEP",
                                          Protein >= 1  ~ "DEP only",
                                          Protein <= -1  ~ "DEP only",
                                          Transcript >= 1  ~ "DEG only",
                                          Transcript <= -1  ~ "DEG only",
                                          TRUE ~ "Not differentially expressed"))


# Designate quadrants for ID
x_mid <- 0
y_mid <- 0

final.ID <- merged.DEPs.DEGs %>% mutate(Quadrant = case_when(Transcript > x_mid & Protein > y_mid   ~ "Q1",
                                                             Transcript <= x_mid & Protein > y_mid  ~ "Q2",
                                                             Transcript <= x_mid & Protein <= y_mid ~ "Q3",
                                                             TRUE ~ "Q4"))
final.ID <- final.ID %>% 
                        arrange(Quadrant)

write.csv(final.ID, output, row.names = F, col.names = T)

# Rename columns 
rownames(merged.DEPs.DEGs) <- merged.DEPs.DEGs$Accession
merged.DEPs.DEGs$Accession <- NULL

merged.DEPs.DEGs <- merged.DEPs.DEGs[, c(2, 1, 3)]
merged.DEPs.DEGs <- merged.DEPs.DEGs %>%
                                    dplyr::rename(x = 1,
                                                  y = 2)

# Get column maximums/minimums for plotting
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMax(merged.DEPs.DEGs)

colMin <- function(data) sapply(data, min, na.rm = TRUE)
colMin(merged.DEPs.DEGs)

# Plot scatter with regression line
my.formula = y ~ x

theme_set(theme_bw())

#TP1 <- ggplot(merged.DEPs.DEGs, aes(x = x, y = y, color = Regulation, group = 1)) +
            geom_point(size = 1, alpha = 0.5) +
            guides(colour = guide_legend(title = "")) +
            geom_vline(xintercept = x_mid) + # plot vertical line
            geom_hline(yintercept = y_mid) + # plot horizontal line
            scale_x_continuous(name = "Transcript log2 fold change", 
                               limits = c(-4.2, 7.3), breaks = seq(-4, 7, 1)) +
            scale_y_continuous(name = "Protein log2 fold change", 
                               limits = c(-2, 3), breaks = seq(-2, 3, 1)) +
            stat_poly_eq(formula = my.formula, method = "lm", 
                         geom = "text",  #or 'label'
                         label.x = 4, #set label location
                         label.y = 3, #set label location
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  
                         rr.digits = 3, parse = TRUE, size = 4, col = "black") +
            ggtitle(Plot.title) +
            theme(aspect.ratio = 1)

#TP1 

#TP3 <- ggplot(merged.DEPs.DEGs, aes(x = x, y = y, color = Regulation, group = 1)) +
            geom_point(size = 1, alpha = 0.5) +
            guides(colour = guide_legend(title = "")) +
            geom_vline(xintercept = x_mid) + # plot vertical line
            geom_hline(yintercept = y_mid) + # plot horizontal line
            scale_x_continuous(name = "Transcript log2 fold change", 
                               limits = c(-9.1, 6.5), breaks = seq(-9, 6, 1)) +
            scale_y_continuous(name = "Protein log2 fold change", 
                               limits = c(-3.5, 4), breaks = seq(-3, 4, 1)) +
            stat_poly_eq(formula = my.formula, method = "lm", 
                         geom = "text",  #or 'label'
                         label.x = -5, #set label location
                         label.y = 3, #set label location
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  
                         rr.digits = 3, parse = TRUE, size = 4, col = "black") +
            ggtitle(Plot.title) +
            theme(aspect.ratio = 1)

#TP3

TP5 <- ggplot(merged.DEPs.DEGs, aes(x = x, y = y, color = Regulation, group = 1)) +
            geom_point(size = 1, alpha = 0.5) +
            guides(colour = guide_legend(title = "")) +
            geom_vline(xintercept = x_mid) + # plot vertical line
            geom_hline(yintercept = y_mid) + # plot horizontal line
            scale_x_continuous(name = "Transcript log2 fold change", 
                               limits = c(-4, 8.1), breaks = seq(-4, 8, 1)) +
            scale_y_continuous(name = "Protein log2 fold change", 
                               limits = c(-3.2, 6), breaks = seq(-3, 6, 1)) +
            stat_poly_eq(formula = my.formula, method = "lm", 
                         geom = "text",  #or 'label'
                         label.x = 4, #set label location
                         label.y = 5, #set label location
                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  
                         rr.digits = 3, parse = TRUE, size = 4, col = "black") +
            ggtitle(Plot.title) +
            theme(aspect.ratio = 1)

TP5


# Save plots to one pdf and R data
MyPlots <- list(TP1, TP3, TP5)

pdf("MC_FC_Prot_Gene_Corr_Plots.pdf")
MyPlots
dev.off()


save(TP1, TP3, TP5, file = "MC_FC_ProtvGene_Corr_Plots.RData")

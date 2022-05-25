# Loading library
library(devtools) 
library(magrittr)
library(dplyr)
library(mixOmics)
library(ggplot2)
library(ggbiplot)
library(scales)
library(ggfortify)
library(ggpubr)
library(gridExtra)

# Load files
data <- read.table("Mcap_Pos_dataset_matrix.csv", header=T, check.names=FALSE, row.names = 1, sep=',')
sample.info.df <- read.table("Sample_class.csv", header=F, sep=',')
sample.info.complete <- read.table("Met_sample_info.csv", header=T, row.names = 1, check.names=FALSE, sep=',')

sample.class.number <- 7

# Set sample list as vector
sample.info <- sample.info.df$V1

# Normalize data
data <- sweep(data,2,as.numeric(sample.info.complete$Weight),FUN='/')

# # Transpose the data to have variables as columns
t.data <- t(data)
dim(t.data)

# Check data for correct data mode (numeric) and no NA values
sapply(t.data, mode)
sum(is.na(t.data))

### START PLS-DA ANALYSIS ###

# Check data for correct dimensions
summary(sample.info)
length(sample.info)
dim(t.data)

# Determine the number of components to retain ncomp -> K−1, where K is the number of classes
k <- sample.class.number-1

# Create PLSDA object
data.plsda <- splsda(t.data, sample.info, ncomp = 6) 
set.seed(30)

# Plot initial PLSDA object
plotIndiv(data.plsda, 
          ind.names = FALSE, 
          legend = TRUE,
          ellipse = TRUE, 
          star = TRUE, 
          title = 'Initial untrained sPLS-DA of MC289 polar metabolomic data',
          X.label = 'PLS-DA 1', 
          Y.label = 'PLS-DA 2')

# Use perf to evaluate the performance of PLS-DA for a large number of components, using repeated k-fold cross-validation
perf.plsda <- perf(data.plsda, validation = "Mfold", folds = 3, 
                   progressBar = FALSE, nrepeat = 50) 

perf.plsda$choice.ncomp

plot(perf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

# Use tune.splsda to determine the number of variables (keepX) to select on each component for sparse PLS-DA
list.keepX <- c(5:10,  seq(20, 1000, 10))
list.keepX # to output the grid of values tested
set.seed(30)
tune.splsda <- tune.splsda(t.data, sample.info, ncomp = 6, 
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   

# Extract the classification error rate averaged across all folds and repeats for each tested keepX value
error <- tune.splsda$error.rate
ncomp <- tune.splsda$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX 
plot(tune.splsda, col = color.jet(ncomp))

# Based on tuning results, run final and tuned sPLS-DA model
splsda.final <- splsda(t.data, sample.info, ncomp = ncomp, keepX = select.keepX)

# Plot PLSDA object
a <- plotIndiv(splsda.final, comp = c(1,2),
          ind.names = FALSE, 
          legend = TRUE,
          ellipse = TRUE, 
          star = TRUE, 
          title = 'sPLS-DA comp 1 & 2')

b <- plotIndiv(splsda.final, comp = c(1,3),
               ind.names = FALSE, 
               legend = TRUE,
               ellipse = TRUE, 
               star = TRUE, 
               title = 'sPLS-DA comp 1 & 3')

c <- plotIndiv(splsda.final, comp = c(1,4),
               ind.names = FALSE, 
               legend = TRUE,
               ellipse = TRUE, 
               star = TRUE, 
               title = 'sPLS-DA comp 1 & 4')

d <- plotIndiv(splsda.final, comp = c(1,5),
               ind.names = FALSE, 
               legend = TRUE,
               ellipse = TRUE, 
               star = TRUE, 
               title = 'sPLS-DA comp 1 & 5')

e <- plotIndiv(splsda.final, comp = c(1,6),
               ind.names = FALSE, 
               legend = TRUE,
               ellipse = TRUE, 
               star = TRUE, 
               title = 'sPLS-DA comp 1 & 6')

# Combine plots
ggarrange(a$graph, b$graph, c$graph, d$graph, e$graph,
          ncol = 2, nrow = 3,
          legend = "right",
          common.legend = TRUE)


### START PCA ANALYSIS ###

# Ensure no constant 0 columns
dim(data)
sapply(data, mode)
sum(is.na(data))
data.df <- as.data.frame(data)
filtered.data <- data %>% 
                            replace(is.na(.), 0) %>%
                            dplyr::mutate(sum = rowSums(across(where(is.numeric))))
filtered.data <- filtered.data[!(filtered.data$sum <= 10000),]
filtered.data$sum <- NULL
filtered.t.data <- t(filtered.data)
filtered.t.data.df <- as.data.frame(filtered.t.data)

# Annotate the data with condition group as labels
data_for_PCA <- merge(sample.info.complete, filtered.t.data.df, by = 0)
rownames(data_for_PCA)=data_for_PCA$Row.names
data_for_PCA$Row.names <- NULL

# Create PCA object
data.pca <- prcomp(data_for_PCA[,c(8:10903)],
                   center = TRUE,
                   scale. = TRUE)

# Summary of the prcomp object
summary(data.pca)

# Structure of the PCA object
str(data.pca)

# PCA plot
ggbiplot::ggbiplot(data.pca,
         data = data_for_PCA,
         obs.scale = 1,
         var.scale = 1,
         ellipse=T,
         ellipse.prob = 0.68,
         var.axes=F,
         circle=F,
         groups = data_for_PCA$Replicate) +
         labs(title="PCA of MC289 polar metabolomic data")

# End
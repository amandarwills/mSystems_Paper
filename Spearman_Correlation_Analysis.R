# Library
library(devtools)
library(magrittr)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(stringr)
library(stringi)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

Metabolite.file <- "MC_DAM_T5_2019"
Amplicon.file <- "MC_seqtab.nochim.csv"
Sequence.ID.file <- "MC_taxa.csv"
Metabolite.ID.file <- "MC_DAM_T1_2019_Knowns.csv"

# Load matrices
Metabolite.data <- read.table(Metabolite.file, header=T, check.names = FALSE, sep='\t')
Amplicon.data <- read.table(Amplicon.file, header=T, check.names = FALSE, row.names = 1, sep=',')
Sequence.info <- read.table(Sequence.ID.file, header=T, check.names = FALSE, sep=',')
Metabolite.info <- read.table(Metabolite.ID.file, header=T, check.names = FALSE, sep=',')

### Filter and format metabolite data ###

# Add IDs to metabolite data
Metabolite.data <- Metabolite.data %>% left_join(Metabolite.info)
Metabolite.data$Matched.Compound[is.na(Metabolite.data$Matched.Compound)] <- Metabolite.data$compoundID 
rownames(Metabolite.data)=Metabolite.data$groupId
Metabolite.data$groupId <- NULL
Metabolite.data$Matched.Compound <- NULL

# Filter metabolites for VIP>=1.5 and FC >=1.5
Filtered.metabolite.data <- Metabolite.data[!(Metabolite.data$VIP <= 1.5),]
Filtered.metabolite.data <- Filtered.metabolite.data %>% 
                      dplyr::filter(FC <= -1.5 | FC >= 1.5)

Filtered.metabolite.data <- Filtered.metabolite.data[, c(1, 10:32)]
rownames(Filtered.metabolite.data)=Filtered.metabolite.data$groupId
Filtered.metabolite.data$groupId <- NULL

# Subset genotype for analysis from data frame
genotype.metabolite.data <- Filtered.metabolite.data[ , grepl( "MC-289" , names( Filtered.metabolite.data ) ) ]

# Rename columns of metabolite data to match amplicon data
met.data.4.corr <- genotype.metabolite.data %>%
                                              rename(Amb_1 = 1,
                                                     Amb_2 = 2,
                                                     Amb_3 = 3,
                                                     HT_1 = 4,
                                                     HT_2 = 5,
                                                     HT_3 = 6)

# Create list of filtered variables (metabolites)
metabolites <- rownames_to_column(met.data.4.corr, var = "Metabolite")
metabolites <- data.frame(metabolites$Metabolite)
metabolites <- metabolites %>%
                              rename(Metabolite = 1)

### Filter and format amplicon data ###

# Subset genotype for analysis from data frame
genotype.amplicon.data <- Amplicon.data[ , grepl( "MC-289_T1" , names( Amplicon.data ) ) ]

# Add sequence info to amplicon data
genotype.amplicon.data <- rownames_to_column(genotype.amplicon.data, var = "Sequence")
genotype.amplicon.ID.data <- genotype.amplicon.data %>% left_join(Sequence.info)
rownames(genotype.amplicon.ID.data)=genotype.amplicon.ID.data$Sequence
genotype.amplicon.ID.data$Sequence <- NULL

ordered.genotype.amplicon.ID.data <- genotype.amplicon.ID.data[order(genotype.amplicon.ID.data$Class, genotype.amplicon.ID.data$Family),]

na.strings=c(""," ","NA")
ordered.genotype.amplicon.ID.data <- na.omit(ordered.genotype.amplicon.ID.data, c("Class"))

ordered.genotype.amplicon.ID.data <- ordered.genotype.amplicon.ID.data[, c(1:6, 9)] 

ordered.genotype.amplicon.ID.data <- ordered.genotype.amplicon.ID.data %>%
                                                                        rename(Amb_1 = 1,
                                                                               Amb_2 = 2,
                                                                               Amb_3 = 3,
                                                                               HT_1 = 4,
                                                                               HT_2 = 5,
                                                                               HT_3 = 6)
class.sums <- ordered.genotype.amplicon.ID.data %>% 
                                          group_by(Class) %>% 
                                          summarize(Amb_1 = sum(Amb_1),
                                                    Amb_2 = sum(Amb_2),
                                                    Amb_3 = sum(Amb_3),
                                                    HT_1 = sum(HT_1),
                                                    HT_2 = sum(HT_2),
                                                    HT_3 = sum(HT_3)
                                          )

# Filter amplicon data for Families with >300 reads & remove cyanobacteria and format for correlation
filtered.classes <- class.sums[rowSums(class.sums[2:7])>500,]
filtered.classes <- data.frame(filtered.classes)
rownames(filtered.classes)=filtered.classes$Class
filtered.classes$Class <- NULL
row_names_df_to_remove <- c("Cyanobacteriia")
amp.data.4.corr <- filtered.classes[!(row.names(filtered.classes) %in% row_names_df_to_remove),]

# Create list of filtered variables (amplicon classes)
classes <- rownames_to_column(amp.data.4.corr, var = "Class")
classes <- data.frame(classes$Class)
classes <- classes %>%
                      rename(Class = 1)

# Transpose each data frame 
t.met.data.4.corr <- t(met.data.4.corr)
t.met.data.4.corr.df <- as.data.frame(t.met.data.4.corr)

t.amp.data.4.corr <- t(amp.data.4.corr)
t.amp.data.4.corr.df <- as.data.frame(t.amp.data.4.corr)

# Join data frames
names(t.met.data.4.corr.df)=str_sub(names(t.met.data.4.corr.df),2)
joined.data <- t.amp.data.4.corr.df %>% merge(t.met.data.4.corr.df, by = 0)
rownames(joined.data)=joined.data$Row.names
joined.data$Row.names <- NULL

# Subset joined.data into 2 data frames - one for each condition
Amb.joined.data <- joined.data[1:3,]
HiT.joined.data <- joined.data[4:6,]

# Create data frame of every combination of metabolite and class with leading data frame name
combinations <- tidyr::crossing(metabolites, classes)
combinations.merge <- combinations %>% rename(Metabolite = 1,
                                              Class = 2)
combinations <- t(combinations)
combinations <- data.frame(combinations)
names(combinations)=str_sub(names(combinations),2)


# Run the Spearman correlation test
Amb.cor.results <- lapply(combinations, function(x) {
                                      (cor.test(Amb.joined.data[,x[1]], Amb.joined.data[,x[2]], 
                                                        method = "spearman",
                                                        exact = FALSE))
})

# Run the Spearman correlation test
spearman.output <- lapply(combinations, function(x) {
                                      (cor.test(joined.data[ ,x[1]], joined.data[ ,x[2]], 
                                                method = "spearman",
                                                exact = FALSE))
})

HiT.cor.results <- lapply(combinations, function(x) {
                                      (cor.test(HiT.joined.data[,x[1]], HiT.joined.data[,x[2]], 
                                                method = "spearman",
                                                exact = FALSE))
})

# Format ambient list into data frame
amb.cor.results.df <- data.frame(matrix(unlist(Amb.cor.results), nrow=length(Amb.cor.results), byrow=TRUE))

amb.cor.results.df <- amb.cor.results.df[, c(1:5)] 
amb.cor.results.df <- data.frame(amb.cor.results.df)
amb.cor.results.df <- amb.cor.results.df %>%
                                    rename(S.statistic = 1,
                                           pavalue = 2,
                                           rho = 3,
                                           null.value = 4,
                                           type = 5)

amb.merged.cor.results.df <- amb.cor.results.df %>% merge(combinations.merge, by = 0)

# Format HiT list into data frame
HiT.cor.results.df <- data.frame(matrix(unlist(HiT.cor.results), nrow=length(HiT.cor.results), byrow=TRUE))

HiT.cor.results.df <- HiT.cor.results.df[, c(1:5)] 
HiT.cor.results.df <- data.frame(HiT.cor.results.df)
HiT.cor.results.df <- HiT.cor.results.df %>%
  rename(S.statistic = 1,
         pavalue = 2,
         rho = 3,
         null.value = 4,
         type = 5)

HiT.merged.cor.results.df <- HiT.cor.results.df %>% merge(combinations.merge, by = 0)

# Unmelt ambient table and format for plotting
amb.unmelted.df <- amb.merged.cor.results.df[, c(7, 8, 4)] 
amb.unmelted.df <- amb.unmelted.df[order(amb.unmelted.df$Class, amb.unmelted.df$Metabolite),]
#amb.unmelted.df <- na.omit(amb.unmelted.df, c("rho"))

amb.unmelted.data <- dcast(amb.unmelted.df, Metabolite ~ Class)

amb.final <- amb.unmelted.data
amb.final.metabolite.list <- amb.final[1]

rownames(amb.final)=amb.final$Metabolite
amb.final$Metabolite <- NULL
amb.final <- lapply(amb.final, as.numeric)
amb.final <- amb.final.metabolite.list %>% merge(amb.final, by = 0)
rownames(amb.final)=amb.final$Metabolite
amb.final$Metabolite <- NULL
amb.final$Row.names <- NULL

# Transform df into matrix and replace NA values with 0's
amb.final[is.na(amb.final)] <- 0
amb.final.matrix <- as.matrix(amb.final)
mode(amb.final.matrix)

# Unmelt HiT table and format for plotting
HiT.unmelted.df <- HiT.merged.cor.results.df[, c(7, 8, 4)] 
HiT.unmelted.df <- HiT.unmelted.df[order(HiT.unmelted.df$Class, HiT.unmelted.df$Metabolite),]
#HiT.unmelted.df <- na.omit(HiT.unmelted.df, c("rho"))

HiT.unmelted.data <- dcast(HiT.unmelted.df, Metabolite ~ Class)

HiT.final <- HiT.unmelted.data
HiT.final.metabolite.list <- HiT.final[1]

rownames(HiT.final)=HiT.final$Metabolite
HiT.final$Metabolite <- NULL
HiT.final <- lapply(HiT.final, as.numeric)
HiT.final <- HiT.final.metabolite.list %>% merge(HiT.final, by = 0)
rownames(HiT.final)=HiT.final$Metabolite
HiT.final$Metabolite <- NULL
HiT.final$Row.names <- NULL

# Transform df into matrix and replace NA values with 0's
HiT.final[is.na(HiT.final)] <- 0
HiT.final.matrix <- as.matrix(HiT.final)
mode(HiT.final.matrix)

# Set Heatmap color scheme
f.amb <- colorRamp2(seq(min(amb.final.matrix), max(amb.final.matrix), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
f.HiT <- colorRamp2(seq(min(HiT.final.matrix), max(HiT.final.matrix), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

# Heatmap 
amb.heatmap.plot <- Heatmap(amb.final.matrix, 
                        col = f.amb,
                        heatmap_legend_param = list(title = "Correlation coefficient", title_gp = gpar(fontsize = 8), 
                                                    labels_gp = gpar(fontsize = 8)),
                        column_title_gp = gpar(fontsize = 8),
                        row_names_max_width = unit(4, "cm"),
                        row_names_gp = gpar(fontsize = 6),
                        column_names_max_height = unit(4, "cm"),
                        column_names_gp = gpar(fontsize = 6),
                        rect_gp = gpar(col = "white", lwd = 1),
                        show_row_dend = FALSE,
                        show_column_dend = FALSE,
                        width = unit(8, "cm"), height = unit(18, "cm"),
                        column_title = "Spearman correlation of amplicon and metabolite data Ambient TP5")
draw(amb.heatmap.plot)

HiT.heatmap.plot <- Heatmap(HiT.final.matrix, 
                            col = f.HiT,
                            na_col = "black",
                            show_heatmap_legend = FALSE,
                            column_title_gp = gpar(fontsize = 8),
                            row_names_max_width = unit(4, "cm"),
                            row_names_gp = gpar(fontsize = 6),
                            column_names_max_height = unit(4, "cm"),
                            column_names_gp = gpar(fontsize = 6),
                            rect_gp = gpar(col = "white", lwd = 1),
                            show_row_dend = FALSE,
                            show_column_dend = FALSE,
                            width = unit(8, "cm"), height = unit(18, "cm"),
                            column_title = "Spearman correlation of amplicon and metabolite data HiT TP5")
draw(HiT.heatmap.plot)

draw(HiT.heatmap.plot + amb.heatmap.plot, ht_gap = unit(0.75, "cm"))

write.csv(metabolites, "metabolite_list.csv", col.names = TRUE, row.names = F)

# End
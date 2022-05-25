#!/usr/bin/env Rscript

matrix.file <-"Pos_Mcap_ElMaven_dataset.csv"
out.prefix <- "Mcap_seqtab.nochim.20percent"
percent.to.change <-20

## Load matrix
data.tmp <- read.table(matrix.file, header=T, check.names = FALSE, row.names=1, sep=',')

## Select data range for randomization
data <- data.tmp[,c(6:89)]
data <- as.matrix(data[,-1])
rownames(data) <- data[,1]

## Set data dimensions
dimensions <- dim(data)
print(paste("Loaded matrix with ",dimensions[1]," rows and ",dimensions[2]," columns",sep=''))

## Convert matrix to vector
data.vector <- c(data)
num.elements <- length(data.vector)
num.elements.to.change <- num.elements*(percent.to.change/100)
min.value <- min(data.vector)
max.value <- max(data.vector)

print(paste("Total elements in matrix: ",num.elements,sep=''))
print(paste("Number of elements in matrix to change: ",num.elements.to.change,sep=''))
print(paste("Min value in matrix: ",min.value,sep = ''))
print(paste("Max value in matrix: ",max.value,sep = ''))


## Plot distribution of original data
pdf(paste(out.prefix,".before.pdf", sep = ''))
hist(log2(data.vector))
dev.off()

## Sample without replacement
index.reordered <- sample(1:num.elements,replace=FALSE)

## Iterate X times and replace a different index each time with a new randomly generated number
for (i in 1:num.elements.to.change) {
  index.to.change <- index.reordered[i]
  rand.value <- round(runif(1, min.value, max.value), digits = 3)
  data.vector[index.to.change] <- rand.value
  print(paste("Round: ",i,"; Index to change: ",index.to.change,"; Random value to add: ",rand.value,sep = ''))
}

## Plot distribution of randomized data
pdf(paste(out.prefix,".after.pdf", sep = ''))
hist(log2(data.vector))
dev.off()

## Convert vector into matrix using dimensions of original matrix
data.final <- matrix(data.vector,ncol=dimensions[2])
colnames(data.final) <- colnames(data)
rownames(data.final) <- rownames(data)

## Add extra column to matrix with the row names
to.write <- data.frame("Name"=rownames(data.final),data.final)

## Write data to table
write.table(to.write, file=paste(out.prefix,".matrix", sep = ''),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
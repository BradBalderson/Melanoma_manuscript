source("tiroshData.R")

#just pulling out the malignant melanoma cells
malignantCells = expressionData[2,]==2
malignantExpression = expressionData[,malignantCells]
melanomaCells = malignantExpression[3,]==0
melanomaExpression = malignantExpression[,melanomaCells]
rownames(melanomaExpression) <- rownames(malignantExpression)
tumorInfo <- melanomaExpression[1,]
melanomaExpression = melanomaExpression[4:dim(melanomaExpression)[1],]
melanomaTPM <- melanomaExpression^2

madGenes <- apply(melanomaTPM, 1, mad)
orderedIndices <- sort.list(madGenes, decreasing=TRUE)
madGenes <- madGenes[orderedIndices]

cumulativeMAD <- double(length(madGenes))
cumulativeMAD[1] <- madGenes[1]
for (i in 2:length(madGenes)){
  cumulativeMAD[i] <- cumulativeMAD[i-1]+madGenes[i]
}

nGenes <- 1:length(madGenes)

plot(nGenes, cumulativeMAD)

#Getting the rate of change
dMAD <- double(length(cumulativeMAD)-1)
for (i in 2:length(cumulativeMAD)){
  dMAD[i-1] <- cumulativeMAD[i]-cumulativeMAD[i-1]
}

plot(dMAD)

x<-0
threshold <- 1
for (i in 1:length(dMAD)){
  if (dMAD[i]>threshold && dMAD[i+1]<=threshold){
    print(dMAD[i])
    print(i)
    x<-i
    break
  }
}

abline(h=threshold, v=x)

plot(nGenes, cumulativeMAD)
abline(h=cumulativeMAD[x],v=x)

#Writing out the reduced, top variable genes
melanomaExpression <- melanomaExpression[orderedIndices,]
outputExpression <- rbind(as.data.frame(tumorInfo), 
                          as.data.frame(melanomaExpression[,]))
write.table(outputExpression, file="tiroshMelanomaData_ordered.txt", 
            sep="\t", quote=FALSE)
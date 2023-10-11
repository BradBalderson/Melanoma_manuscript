source("tiroshData.R") #Loading the data
group_info <- read.delim("finalClusters.csv")
group_info <- t(group_info)

matchingCols <- match( group_info["cellCode",], colnames(expressionData))
#reverseMatchCols <- match( group_info[1,], colnames(expressionData))
#unmatched <- group_info[1,is.na(reverseMatchCols)]
#matched <- group_info[1,!is.na(reverseMatchCols)]
#cols <- matchingCols[!is.na(matchingCols)]
malignantExpression <- expressionData[,matchingCols]

colnames(group_info) <- group_info[1,]
group_info <- t(group_info[2,])
x <- tail( rbind( malignantExpression[,], group_info[,] ) )

malignantExpression <- rbind( malignantExpression[,], group_info[,] ) 
rownames(malignantExpression)[length(rownames(malignantExpression))] <- "cluster"
malignantExpression <- malignantExpression[4:(dim((malignantExpression)[1])),]

isCluster0 <- malignantExpression["cluster",]==0
cluster0 <- malignantExpression[1:(dim(malignantExpression)[1]-1),isCluster0]

isCluster1 <- malignantExpression["cluster",]==1
cluster1 <- malignantExpression[1:(dim(malignantExpression)[1]-1),isCluster1]

DEGenes = data.frame(matrix(data=0.0, (nrow=dim(cluster0)[1]), ncol=3))
rownames(DEGenes) <- rownames(cluster0)
colnames(DEGenes) <- c("pvalue","corrected_pvalue","rankDiff(cluster0-cluster1)")#"paired_pvalue","corrected_pvalue","corrected_paired_pvalue")
                     
for (i in 1:nrow(cluster0)) {
  
  threwError <- FALSE
  tryCatch(testResult <- wilcox.test(as.numeric(cluster0[i,]), as.numeric(cluster1[i,]), 
                                     conf.int=TRUE, conf.level=0.9, alternative="two.sided"), 
           warning= function(w) {
             DEGenes[i,1] <<- 100
             DEGenes[i,3] <<- 0.0
             threwError <<- TRUE
           }, 
           error = function(e) {
             DEGenes[i,1] <<- 100
             DEGenes[i,3] <<- 0.0
             threwError <<- TRUE
            }
  )
  
  if (!threwError) {
    DEGenes[i,1] <- testResult$p.value
    DEGenes[i,3] <- testResult$estimate
  }
  #DEGenes[i,2] <- wilcox.test(as.numeric(cluster0[i,]), as.numeric(cluster1[i,]), paired=T)$p.value
}

##Now the corrected pvalues##
DEGenes[,2] <- p.adjust(DEGenes[,1], method=p.adjust.methods[7], n=length(DEGenes[,1]))

#Just writing out the gene info
write.csv(DEGenes, file='TiroshDEGenes_all.csv')

###Now determining which are DE between the groups####
DEGenes[,1] <- as.numeric(DEGenes[,1])
DEGenes[,2] <- as.numeric(DEGenes[,2])
DEGenes[,3] <- as.numeric(DEGenes[,3])

sigGenesBool <- DEGenes[,2]<=0.001
sigGenes <- DEGenes[sigGenesBool,]
nonDEGenes <- DEGenes[!sigGenesBool,]

cluster0GenesBool <- sigGenes[,3]>0
cluster0GenesInfo <- sigGenes[cluster0GenesBool,]

cluster1GenesBool <- sigGenes[,3]<0
cluster1GenesInfo <- sigGenes[cluster1GenesBool,]

cluster0GenesInfo$"cluster" <- matrix(0.0,dim(cluster0GenesInfo)[1],1)
cluster1GenesInfo$"cluster" <- matrix(1,dim(cluster1GenesInfo)[1],1)
nonDEGenes$"cluster" <- matrix('NA',dim(nonDEGenes)[1],1)

orderedResults <- rbind(cluster0GenesInfo, cluster1GenesInfo, nonDEGenes)

write.csv(orderedResults, "TiroshDEGenes_allFormatted.csv", quote=FALSE)

###Getting the DE genes for each state into appropriate format for input to GSEA
cluster0Genes <- as.vector(rownames(cluster0GenesInfo))
cluster1Genes <- as.vector(rownames(cluster1GenesInfo))

write.table(cluster0Genes, "cluster0Genes2.txt", row.names = FALSE, col.names=FALSE, quote=FALSE)
write.table(cluster1Genes, "cluster1Genes2.txt", row.names = FALSE, col.names=FALSE, quote=FALSE)

library("org.Hs.eg.db")

cluster0GenesIDs <- mget(x=cluster0Genes, envir=org.Hs.egALIAS2EG)
#df <- data.frame(matrix(NA, 2, max( sapply( list(cluster0Genes, cluster1Genes), length))))
#df[1,1:length(cluster0Genes)] <- cluster0Genes
#df[2,1:length(cluster1Genes)] <- cluster1Genes
#df

#rownames(df) <- c("0", "1")
#df <- cbind(geneSetDescription=c("DBScan_Cluster_0","DBScan_Cluster_1"), df)
#write.table(df, file = 'TiroshClusterGeneSets.gmt', row.names = TRUE, col.names=FALSE, na = '', sep='\t')

#colnames(group_info) <- group_info[1,]
#group_info <- group_info[2,]
#x <- tail( rbind( malignantExpression[,], group_info[,] ) )

#malignantExpression <- rbind( malignantExpression[,], group_info[,] ) 
#rownames(malignantExpression)[length(rownames(malignantExpression))] <- "cluster"
#malignantExpression <- malignantExpression[4:(dim((malignantExpression)[1])),]

#Factor levels for grouping the data
#meF <- factor(group_info, levels = c("0","1"))

##Ranking the data##
#rankedME <- apply(malignantExpression, 1, rank)

##Summing the ranks according to cluster labels##
#summedRanks <- apply(rankedME, 2, function(x) tapply(x, group_info, sum))

##This is the actual function##

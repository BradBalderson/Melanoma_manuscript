library(amap)

#####reading in the relative expression information#####
relExprs <- read.table('tiroshRelExprs_sigNegCorrs.txt', sep='\t', header=TRUE, row.names=1)

#performing the hiearchical clustering
clustResults <- hcluster(t(relExprs), "euclidean", link="complete", nbproc=20)

#getting labels
members <- cutree(clustResults, k=4)

#saving the new labels
write.table(members, file='relExprsSigNegCorrLabels.txt', quote=F, col.names=F, row.names=F)

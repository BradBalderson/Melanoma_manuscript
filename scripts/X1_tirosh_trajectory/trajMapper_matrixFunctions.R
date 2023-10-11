###Helper functions####
library(foreach)
library(doMC)
library(Matrix)

getGeneOrderings <- function(exprs, genes, nCores = 1){

   nCols <- dim(exprs)[2]
   nGenes <- length(genes)
   
   if (nCores>1){
     registerDoMC(nCores)
   }

   print(nCols)   
   orderings <- foreach (j=1:nCols, .combine='c') %dopar% {
    geneOrderings <- Matrix( rep(0, nGenes*nGenes ), 
                        nGenes, nGenes, sparse = TRUE)
    rownames(geneOrderings) <- genes 
    colnames(geneOrderings) <- genes
    for (i in 2:length(genes)){
      range <- 1:(i-1)
      geneOrderings[range,i] <- exprs[genes[range],j]>exprs[genes[i],j]
    }
    
    c(geneOrderings)
  } 
     
  return(orderings)
}

#Do not re-run, best candidate genes done.
#orderingPrac <- getGeneOrderings(loggedExpression, c("MITF","AXL"), nCores = 42)

###Just need to change up below so that get significant pairs across time points####

#countMatPrac <- Reduce('+', orderingPrac[c(which(timePoint==1),which(timePoint==2))])
#time1Count <- Reduce('+', orderingPrac[which(timePoint==1)])
#time2Count <- Reduce('+', orderingPrac[which(timePoint==2)])
#timeCounts <- list(time1Count, time2Count)
getTimeCounts <- function(geneOrderings, timePoints, nCores=1){
  nPoints <- max(timePoints)
  if (nCores>1){
    registerDoMC(nCores)
  }
  
  timeCounts <- foreach (i=1:nPoints, .combine='c') %dopar% {
    Reduce('+', geneOrderings[which(timePoints==i)])
  }
  
  return(timeCounts)
}

convertToPairCounts <- function(countMatrix, genes) {
  nGenes <- length(genes)
  nPairs <- dim(combn(nGenes,2))[2]
  pairCounts <- integer(nPairs)

  pairNames <- character(nPairs)
  pairN <- 1
  for (i in 2:nGenes){
    for (j in 1:(i-1)){
      pairNames[pairN] <- paste(genes[i], genes[j],sep="_")
      pairCounts[pairN] <- countMatrix[j,i]
    
      pairN <- pairN+1
    }
  }

  names(pairCounts) <- pairNames
  print(pairNames[1])
  return(pairCounts)
}

getOrderingsOfPairs <- function(exprs, pair1, pair2, nCores=1, sep="_"){
	nPairs <- length(pair1)
	nCols <- dim(exprs)[2]

	###getting the pair names##
	pairNames <- character(nPairs)
	for (i in 1:nPairs){
		pairNames[i] <- paste(pair1[i],pair2[i],sep=sep)
	}

	###registering for multithreading##
	if (nCores>1){
		registerDoMC(nCores)
	}

	###getting relative expression###
	orderings <- foreach(n=1:nCols, .combine='rbind') %dopar% {

		###comparing all the genes simulataneously
		as.integer(exprs[pair1,n]>exprs[pair2,n])

	}
	
	###some final formatting
	colnames(orderings) <- pairNames
	rownames(orderings) <- NULL
	orderings <- as(orderings, 'dgCMatrix')
	
	return(orderings)
}

getPairedOrderings <- function(exprs, genes,  nCores=1) {
  nGenes <- length(genes)
  nPairs <- dim(combn(nGenes,2))[2]
  nCols <- dim(exprs)[2]

  ###Getting the pair names###
  pairNames <- character(nPairs)
  pairN <- 1
  for (i in 2:nGenes){
      for (j in 1:(i-1)){
	 pairNames[pairN] <- paste(genes[i], genes[j],sep="_")
	 pairN <- pairN+1
      }
  }

  ###registering for multithreading##
  if (nCores>1){
     registerDoMC(nCores)
  }

  ##getting the ordering##
  orderings <- foreach (n=1:nCols, .combine='rbind') %dopar% {

    pairCounts <- integer(nPairs)#Matrix( rep(0, nPairs ),
                  #      1, nPairs, sparse = TRUE)

    pairN <- 1
    for (i in 2:nGenes){
      for (j in 1:(i-1)){
	
	result <- exprs[genes[i],n]>exprs[genes[j],n] #n is sample

	if (result) { pairCounts[pairN]<-1 } 
	#print(genes[i])
	#print(genes[j])
        pairN <- pairN+1
      }
    }

   #print(pairCounts)
   c(pairCounts)
  }

  #print(orderings)
  #print(pairNames)
  colnames(orderings) <- pairNames
  rownames(orderings) <- NULL
  #print(orderings)
  orderings <- as(orderings, 'dgCMatrix')
  return(orderings)
}

#pairCounts <- convertToPairCounts(countMatPrac, topGenes)
#timePairCounts <- lapply(1:2, function(i) convertToPairCounts(timeCounts[[i]], topGenes))

#overallProbs <- pairCounts/946

z.prop = function(x1,x2,n1,n2){
  numerator = (x1/n1) - (x2/n2)
  p.common = (x1+x2) / (n1+n2)
  denominator = sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
  z.prop.ris = numerator / denominator
  return(z.prop.ris)
}

getSigDiffProps <- function(timePairCounts, timePoints, nCores=1){

  if (nCores>1){
    registerDoMC(nCores)
    print(nCores)
  }

  nPoints <- length(timePairCounts)
  nPairs <- length(timePairCounts[[1]])
  nCells <- length(timePoints)
  
  #print(nPoints)
  #print(nPairs)
  pairWisePVals <- foreach (i=1:(nPoints-1), .combine='c') %dopar% {
    pVals <- Matrix(rep(0, nPairs*(nPoints-i)), nPairs, (nPoints-i))
    rownames(pVals) <- names(timePairCounts[[1]])
    nTime1 <- length(which(timePoints==i))
    
    for (j in (i+1):(nPoints-1)){
      nTime2 <- length(which(timePoints==j))

      for (n in 1:nPairs) {
        sucCounts <- c(timePairCounts[[i]][n], timePairCounts[[j]][n])
        trialCounts <- c(nTime1,nTime2)
        pVals[n,j] <- 2*pnorm(-abs(z.prop(sucCounts[1],sucCounts[2],trialCounts[1],trialCounts[2])))
      }
      
    }
    
    c(pVals)
  }
  
  return(pairWisePVals)
}
  
#pValsCorr <- p.adjust(pVals, method="bonferroni")

#length(which(pValsCorr<0.000001))






library(monocle)
library(Matrix)
source('tiroshData_madOrdered.R')

genes <- rownames(expressionData)!='tumor'
expressionMatrix <- expressionData[genes,]
expressionMatrix <- as((2^as.matrix(expressionMatrix))-1, 'sparseMatrix')

tumorTypes <- as.character(t(expressionData[dim(expressionData)[1],]))
tumors <- data.frame(tumorTypes, stringsAsFactors = TRUE)
tumors <- transform(tumors, tumorTypes = as.character(tumorTypes))
rownames(tumors) <- colnames(expressionMatrix)
phenoData <-  new("AnnotatedDataFrame", data = tumors)
featureData <- new("AnnotatedDataFrame", data =  as.data.frame(rownames(expressionMatrix), nm='gene_short_name'))
rownames(featureData) <- rownames(expressionMatrix)

HSMM <- newCellDataSet(expressionMatrix, phenoData = phenoData, featureData = featureData,
            expressionFamily=tobit(Lower=0.1))

# Next, use it to estimate RNA counts
rpc_matrix <- relative2abs(HSMM, method = "num_genes")

# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = phenoData,
                       featureData = featureData,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM), #genes expressed in more than 10 cells
                      num_cells_expressed >= 100))

###removing low-quality cells###
print(head(pData(HSMM)))

valid_cells <- row.names(subset(pData(HSMM),
                                num_genes_expressed>2000))
HSMM <- HSMM[,valid_cells]
##plotting distributions..##
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(HSMM), color = 'k', geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)
##From the above plot looks like monocle managed to normalise
##Tirosh data!!!!
##removing cells which don't occur within reasonable bounds
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < upper_bound] ##didn't actually remove anything
HSMM <- detectGenes(HSMM, min_expr = 0.1)

##checking data roughly normal
# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")
#YEAAAAH data is roughly normal

###Clustering####
#clustering cells according to type(prolif vs invasive)
exprsMITF <- exprs(HSMM)['MITF',]
mitfOrderedIndices <- sort(exprsMITF, decreasing = TRUE)
exprsMITF <- exprsMITF[mitfOrderedIndices]

table(exprsMITF)

prolif_id <- row.names(subset(fData(HSMM), gene_short_name == "MITF"))
invasive_id <- row.names(subset(fData(HSMM),
                             gene_short_name == 'AXL'))
                               
                            #"MT2A"||'RPS4Y1'||'IFI6'))


cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Prolif", classify_func =
                     function(x) 
                       { x[prolif_id,] > 10 })

cth <- addCellType(cth, "invasive", classify_func = function(x)
{ x[prolif_id,] < 1})

HSMM <- classifyCells(HSMM, cth)
table(pData(HSMM)$CellType)

##Below getting genes that covarry with invasive vs prolif signature
#marker_diff <- markerDiffTable(HSMM[expressed_genes,],
#              cth,
#              residualModelFormulaStr = "~tumorTypes+num_genes_expressed",
#              cores = 3)

#candidate_clustering_genes <-
#  row.names(subset(marker_diff, qval < 0.1))
#marker_spec <-
#  calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
#head(selectTopMarkers(marker_spec, 3), n=9)

##getting genes to use for clustering (basically all)
#semisup_clustering_genes <- 
#  row.names(marker_diff)[order(marker_diff$qval)][1:3000]
semisup_clustering_genes <- as.character(read.table('cellLineDEGenes.txt', header=FALSE)[,1])
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)
#dispersion approach
#disp_table <- dispersionTable(HSMM)
#unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
#HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
#plot_ordering_genes(HSMM)

##looking at the principle components variance
plot_pc_variance_explained(HSMM, return_all = F) # norm_method='log'

#performing clustering based on reduced dimensions
#note HSMMe is below with num_dim=2
HSMM3 <- reduceDimension(HSMM, max_components = 6, num_dim = 3,
                         #norm_method = 'log',
                        reduction_method = "tSNE", 
                        residualModelFormulaStr = 
                          '~tumorTypes',
                        verbose = T)
HSMM3 <- clusterCells(HSMM3, num_clusters = 5)
plot_cell_clusters(HSMM3, 1, 2, color = "Cluster")#, 
clusters <- HSMM3$Cluster
                   #markers=c('MITF','CTNNB1'))

#plot_cell_trajectory(HSMMe, color_by = "CellType")#,
###Learning the tragectory####
##first find DE genes between the clusters to learn tragectory
#diff_test_res <- differentialGeneTest(HSMMe[expressed_genes,],
#                        fullModelFormulaStr = "~Cluster")#,
                       # reducedModelFormulaStr = "~tumorTypes")
#ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

#HSMMe <- setOrderingFilter(HSMMe, ordering_genes)
#plot_ordering_genes(HSMMe)

##Dimension Reduction for tragectory plotting##
semisup_clustering_genes <- as.character(read.table('cellLineDEGenes.txt', header=FALSE)[,1])
HSMM3 <- setOrderingFilter(HSMM3, semisup_clustering_genes)
HSMMDDR <- reduceDimension(HSMM3, max_components = 2, #num_dim=5,
                         reduction_method = 'DDRTree', 
                         residualModelFormulaStr = 
                           '~tumorTypes',
                         verbose = T,
                         auto_param_selection = TRUE)
                         #ncenter=1000)
HSMMOrdered <- orderCells(HSMMDDR)
HSMMOrdered$'Cluster' <- clusters
HSMMOrdered <- orderCells(HSMMOrdered, root_state=7)

plot_cell_trajectory(HSMMOrdered, color_by = "tumorTypes")#,
                    #markers=c('MITF','CTNNB1'), use_color_gradient = TRUE)
source('monocle_plottingFunctions.R')
plot_cell_clusters(HSMMOrdered, 2, 1, color = "State",
                  markers=c('NGFR','AXL'))
plot_complex_cell_trajectory(HSMMOrdered)

###calling DE genes between clusters####
diff_test_res <- differentialGeneTest(HSMMOrdered[expressed_genes,],
                                      fullModelFormulaStr = "~Cluster",
                                      reducedModelFormulaStr="~tumorTypes")

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.05)

###Let's try look at the genes in PCA space####
#ensuring to use only expressed genes...
#matches <- match(semisup_clustering_genes, expressed_genes)
#notNA <- !is.na(matches)
#matches <- matches[notNA]
#importantGenes <- expressed_genes[matches]
#logExprsDiffGenes <- log(exprs(HSMM[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
#scaledDiffGenes <- Matrix::t(scale(Matrix::t(logExprsDiffGenes)))

#melt(Matrix::t(scale(Matrix::t(L))))
#some for matting for plotting
#diffGenesForPlotting <- melt(scaledDiffGenes)

#just checking if distribution makes sense
#qplot(value, geom = "density", data = diffGenesForPlotting) +
#  stat_function(fun = dnorm, size = 0.5, color = 'red') +
#  xlab("Standardized log(FPKM)") +
#  ylab("Density")
####running PCA####



#Plotting on root state
#GM_state <- function(cds){
#  if (length(unique(pData(cds)$State)) > 1){
#    T0_counts <- table(pData(cds)$State, pData(cds)$Cluster)[,"7"]
#    return(as.numeric(names(T0_counts)[which
#                                       (T0_counts == max(T0_counts))]))
#  } else {
#    return (1)
#  }
#}

#HSMMO2 <- orderCells(HSMMOrdered, root_state = 7)
#plot_cell_trajectory(HSMMO2, color_by = "Pseudotime")

#Looking at genes of interest expression between states
HSMMO2 <- HSMMOrdered
blast_genes <- row.names(subset(fData(HSMMO2),
                                gene_short_name %in% c("EZH2", "AXL", "MITF")))
plot_genes_jitter(HSMMO2[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)

#looking at how genes expression changes in pseudo time
HSMM_expressed_genes <-  row.names(subset(fData(HSMMO2),
                                  num_cells_expressed >= 10))
HSMM_filtered <- HSMMO2[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                gene_short_name %in% c("AXL", "MITF", "ERBB3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")

#looking at genes dependent on pseudo time
plot_pseudotime_heatmap(HSMMO2[semisup_clustering_genes,],
                        num_clusters = 5,
                        cores = 1,
                        show_rownames = T)

###Looking at individual branches###
HSMMO3 <- setOrderingFilter(HSMM3, expressed_genes)
branchGenes <- branchTest(HSMMO3, fullModelFormulaStr = "~tumorTypes+sm.ns(Pseudotime, df = 3)*Branch",
           reducedModelFormulaStr = "~tumorTypes+sm.ns(Pseudotime, df = 3)",
           branch_point = 4)
heatMapResults <- plot_genes_branched_heatmap(HSMMO2[row.names(subset(branchGenes,
                                                  qval < 0.05)),],
                            branch_point = 4,
                            num_clusters = 11,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T,
                            return_heatmap = TRUE)

#separating the genes into their clusters
geneClusters <- heatMapResults$annotation_row

getClusters <- function(geneClusters, clusterNumber){
  
  inCluster <- geneClusters[,1]==clusterNumber
  geneCluster <- rownames(geneClusters)[inCluster]
  
  return(geneCluster)
}

cluster1 <- getClusters(geneClusters, 1)
cluster2 <- getClusters(geneClusters, 2)
cluster3 <- getClusters(geneClusters, 3)
cluster4 <- getClusters(geneClusters, 4)
cluster5 <- getClusters(geneClusters, 5)
cluster6 <- getClusters(geneClusters, 6)
cluster7 <- getClusters(geneClusters, 7)
cluster8 <- getClusters(geneClusters, 8)
cluster9 <- getClusters(geneClusters, 9)
cluster10 <- getClusters(geneClusters, 10)
cluster11 <- getClusters(geneClusters, 11)

printCluster <-function(clusterN){
  for (i in 1:length(clusterN)){
    x <- cat(clusterN[i], "\n")
    if (length(x)!=0){
      print(x)
    }
  }
}

#put the diferent cluster stuff into david, from what I can tell, it would appear that
#hypoxia is responsible for inducing the neural-crest state. Meanwhile interactions with
#immune cells induces the undifferentiated state (due to ifn-gamma signalling). 
#But for the melanocyte-like lineage dosn't appear to be any GO terms for signalling
#pathways overrepresented. Therefore, it seems like the base state is the 
#melanocyte, while the stressors induce going from the melanocyte like
#to the undifferentiated forms... 

#I have just reliazed I don't know acutally what these cell states represent yet.
#Going to DE between the cells at the end of the branches between the three states
#to properly classify each. 

orderedExpressionData <- exprs(HSMMOrdered)

#Getting transitory expression
isCluster6 <- HSMMOrdered$State==6
isCluster7 <- HSMMOrdered$State==7
isTransitLike <- logical(length(isCluster6))
for (i in 1:length(isCluster6)){
  if (isCluster6[i] || isCluster7[i]){
    isTransitLike[i]=TRUE
  }
}

transitLikeExpression <- orderedExpressionData[,isTransitLike]

#getting undifferentiated cells
isCluster3 <- HSMMOrdered$State==3
cluster3Pseudotime <- HSMMOrdered$Pseudotime[isCluster3]
cluster3Expression <- orderedExpressionData[,isCluster3]

isNeuralCrestLike <- cluster3Pseudotime>=6.4
neuralCrestLikeExpression <- cluster3Expression[,isNeuralCrestLike]
print(dim(neuralCrestLikeExpression))
#getting the melanocyte like cells
isCluster1 <- HSMMOrdered$State==1
cluster1Pseudotime <- HSMMOrdered$Pseudotime[isCluster1]
cluster1Expression <- orderedExpressionData[,isCluster1]

isMelanocyte <- cluster1Pseudotime>7.8
melanocyteExpression <- cluster1Expression[,isMelanocyte]

#now visualling the cells as I'v placed them to make sure they're representative
#of the end of the branches in the trajectory
neuralCrestCells <- colnames(neuralCrestLikeExpression)
melanocyteCells <- colnames(melanocyteExpression)
transitCells <- colnames(transitLikeExpression)

ncColindex <- match(colnames(orderedExpressionData), neuralCrestCells)
melColindex <- match(colnames(orderedExpressionData), melanocyteCells)
transitColindex <- match(colnames(orderedExpressionData), transitCells)

diffStatus <- character(length(ncColindex))
for (i in 1:length(diffStatus)){
  #if not classified, then must be transitory between one of the states
  if (is.na(ncColindex[i]) && is.na(melColindex[i]) && is.na(transitColindex[i])){
    diffStatus[i]<-"differentiating"
  
  } else if (!is.na(ncColindex[i])){
    diffStatus[i]<-"neural-crest"
    
  } else if (!is.na(melColindex[i])) {
    diffStatus[i]<-"melanocyte"
      
  } else if(!is.na(transitColindex[i])){
    diffStatus[i]<-"transit"
    
  }
}

HSMMOrdered$"diffStatus" = diffStatus
#recolouring the trajectory with the different classifications
plot_cell_trajectory(HSMMOrdered, color_by = "diffStatus") +
  facet_wrap(~diffStatus, nrow = 1)
plot_cell_clusters(HSMMOrdered, x=2, y=1, color_by = "diffStatus")
###Looks like they've been separated nicely.. 

#now going to call DE genes between the states..
##needed to run this on the server since takes ages to run.
##on server in Tirosh/monocle directory is Rscript I ran to get the 
#DE genes between the groups.
#saved now to this directory as .RData so should auto-load
#ncMelSigGenes - neural-crest vs. melanocyte DE genes
#ncUndiffSigGenes - neural-crest vs. undiff DE genes
#melUndiffSigGenes - melanocyte vs. undiff DE genes
#is a mistake, had the melanocyte and undiff group mixed around.
#fixed below
#x <- ncUndiffSigGenes
#ncUndiffSigGenes <- ncMelSigGenes
#ncMelSigGenes <- x
#undiffMelSigGenes <- melUndiffSigGenes
#turns out undiff is nc and nc is a transitory state, fixed below

#let's get the DE genes down to around 100 each representing genes uniquely
#DE in that group
########for below don't forget names of diffStatus mixed around for undiff and melanocyte###
#####transit unique:
#ncUndiff <- subset(ncUndiffSigGenes, `diffStatusneural-crest`>0) #AXL not here
#ncMel <- subset(ncMelSigGenes, `diffStatusundifferentiated`<0) #AXL in this subset
#ncUnique <- intersect(rownames(ncMel), rownames(ncUndiff)) 
#ncUndiffUnique <- setdiff(ncUndiffGenes, ncMelGenes)
#ncMelUnique <- setdiff(ncMelGenes, ncUndiffGenes)

#ncMelQval <- ncMel[ncUnique,]$qval #getting genes with highest average qval
#ncUndiffQval <- ncUndiff[ncUnique,]$qval
#ncMeanQval <- (ncMelQval+ncUndiffQval)/2
#ncUniqueTop <- rownames(ncMel[ncUnique,][order(ncMeanQval),])[1:25]
###DONT RUN AGAIN####
transitNC <- ncUndiff
transitMel <- ncMel
transitUnique <- ncUnique
transitNCUnique <- ncUndiffUnique
transitMelUnique <- ncMelUnique

transitMelQval <- ncMelQval
transitNCQval <- ncUndiffQval
transitMeanQval <- ncMeanQval
transitUniqueTop <- ncUniqueTop

#melanocyte-like unique:
melNC <- subset(ncMelSigGenes, `diffStatusundifferentiated`>0)
melUndiff <- subset(undiffMelSigGenes, `diffStatusundifferentiated`>0)
melUnique <- intersect(rownames(melNC), rownames(melUndiff)) #795, so need qval cutoff

melNCQval <- melNC[melUnique,]$qval #getting genes with highest average qval
melUndiffQval <- melUndiff[melUnique,]$qval
melMeanQval <- (melNCQval+melUndiffQval)/2
melUniqueTop <- rownames(melNC[melUnique,][order(melMeanQval),])[1:25]

#neuralCrest-like unique:
#undiffNC <- subset(ncUndiffSigGenes, `diffStatusneural-crest`<0)
#undiffMel <- subset(undiffMelSigGenes, `diffStatusundifferentiated`<0)
#undiffUnique <- intersect(rownames(undiffNC), rownames(undiffMel)) #795, so need qval cutoff

#undiffNCQval <- undiffNC[undiffUnique,]$qval #getting genes with highest average qval
#undiffMelQval <- undiffMel[undiffUnique,]$qval
#undiffMeanQval <- (undiffNCQval+undiffMelQval)/2
#undiffUniqueTop <- rownames(undiffNC[undiffUnique,][order(undiffMeanQval),])[1:25]
ncTransit <- undiffNC
ncMel <- undiffMel
ncUnique <- undiffUnique

ncTransitQval <- undiffNCQval
ncMelQval <- undiffMelQval
ncMeanQval <- undiffMeanQval
ncUniqueTop <- undiffUniqueTop

#From what I can tell looks more like the undiff category is the neural-crest
#cells due to expression of neural-crest like genes and the neural-crest category is
#something else..potentially some partially differentiated neuronal-like cell'''

#AXL upregulated undiff vs. mel group.
#AXL upregulated in nc(unknown) group compared with mel group.
#MITF upregulated in undiff group compared with nc(unknown") group. 
#MITF upregulated in mel group compared with nc group.

#based off this, looks like:
#undiff = AXL mid MITF mid group
#nc = AXL high MITF low group
#mel = MITF high AXL low group

#mel group is definitely mel group, lots of melanocyte specific genes expressed

#looks like the undiff group is actually the neural-crest like group #cell fate 1
#due to:
#SOX5/SOX10 expression https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302844/
#^^involved in neural-crest migration/specification module
#expression of ANKRD10, a regulator in wnt signalling, signal nc differentiation
#NANOG and POU5F1 expression <- suggests dedifferentiation
#FGF5 and FGFR2 expression <- suggest activity of fibroblast GF pathway initiates nc differentiation
#upregulation of MAPK <- MAP3K13, MAP1Lc3C, MAP3K3, MAPK13, MAPKBP1 <- perhaps suggests BMP activation
#thereby promoting nc formation
#HIF3A expression suggests hypoxia <- which is known to regulate stemness and metastasis
#including regulators such as POU5F1 and TWIST (we know TWIST involved in EMT program in NC)
#POU2AF1 also overexpressed <- suggests expression of OCT1 and OCT2 since it only binds to these
#exprssion of NBPF1 (mainly), NBPF10, and NBPF14 represses beta-catenin for wnt-signalling for 
#proliferation, suggests this might promote nc invasive program with wnt-signalling active
#NBPF1 involved in proliferation/growth suppresion. inihibitor of canonaical wnt-signalling

#####Trying to understand the unknown cluster##### #base state
####differences with melanocyte like group#####
#I think the nc group is something weird due to an interaction effect between
#melanoma dedifferentiation and ifn-gamma signalling from the immune system. 
#This is due to http://www.pnas.org/content/111/6/2301 which shows an effect of
#ifn-gamma on pigmentation dependent on IRF1, which is overexpressed in the 
#undiff and neural-crest group compared with the melanocyte group.
#however, the neural-crest(actually the unknown group) is exposed to
#ifn-gamma according the the DE genes called along that branch due to pseudotime,
#whereas on the undiff(actually neural-crest) this signalling isn't present, instead
#there is hypoxia. Apparently this IRF1 role is independent of MITF expression.
#From the above paper, this ifn-gamma program mediates pigmentation by arresting 
#melanasome in the early maturation process by decreasing flux of proteins trafficked to
#the early melanosomes.
#Looks like there are two hypopigmentation programs in melanoma:
#one mediated by dedifferentiation to neural-crest via hypoxia, 
#the other by ifn-gamm signalling from the immune system resulting in 
#a cell type with stem-like properites (due to being partially dediffereniated) #nope, both ifn-gamma signalling appears upregulated in both
                                                                                #of the clusters that are not melanocytic.. 
#but also slightly high in MITF since is depigmented by different means.. 
#IRF1 upregulated in mel vs unknown
#IRF1","IRF2BP2","IRF3" upregulated in actual nc vs mel
#With visible hypopigmentation under ifn-gamma treatment, showed IRF1 was upregulated.
#with no hypopigmentation IRF1 low. 
#Under ifn-gamma induced hypopigmentation showed HLA-DRB1 extremely highly over expressed.
#in the unknown phenotype in this data, multiple HLAs are overexpressed compared with mel type
#including HLA-DRB1:
#HLA-DMA, HLA-DPA1, HLA-DPB1, HLA-DRA, HLA-DRB1, HLA-DRB5. #these also overexpressed in actual NC cluster
#These genes are MHC class 2 genes

#KLF4 and KLF10 also upregulated in unknown cluster vs. melanocyte 
#TRAF1 upregulated in unknown, is member of TNFR associated factor(TRAF)
#family. Associates with TRAF2 for TNF-alpha meidated activation of MAPK8 and NFkB.
#BIRC2 and BIRC3 act downstream. 
#In the unknown cluster, BIRC3 and BIRC5 are upregulated, suggesting some TNFalpha signalling as well. 
#TNFalpha shown to increase melanoma invasion and oppose attachment to ECM via fibronectin. 

#FOSL1 is upregulated, this binds JUN to form AP-1, which verfaille suggest as a master regulator of melanoma
#metastasis.

####differences with neural-crest group#####
#difficult to interpret.. dosn't make much sense.. 

#I think if i can understand what the intermediate state between the three states mean and
#take the DE between this intermediate and the other 3, should get an indication what
#the third state might be. From what I can tell however, seems to be some kind of abberant
#cell type which arises due to TGF-beta signalling of an intermediate between NC and melanocyte,
#With TGF-beta reducing melanocyte pigmentation independent of MITF and inducing 
#chondrocyte differenitation in NC cells. 



####To help with characterising which state is which, going to use the same 
####differentiation trajectory score determined by Tsoi to see if this corresponds to 
####the states by Tsoi and others somewhat. 

###No instead going to look at an overlap between gene lists with the four groups identified between
###each to see if these somehow map to one another. 
#calculating the DE genes enriched in each group, using wilcoxin rank-sum test.
source('../Tsoi/hiearchialClustering/DEFunctions.R')
notNCColumns <- is.na(match(colnames(orderedExpressionData),neuralCrestCells))
notTransColumns <- is.na(match(colnames(orderedExpressionData),transitCells))
notMelanColumns <- is.na(match(colnames(orderedExpressionData),melanocyteCells))
diffColumns <- logical(length(notMelanColumns))
for (i in 1:length(notMelanColumns)){
  if (notNCColumns[i] && notTransColumns[i] && notMelanColumns[i]){
    diffColumns[i] <- TRUE 
  }
}

notNCExpression <- orderedExpressionData[,notNCColumns]
notTransExpression <- orderedExpressionData[,notTransColumns]
notMelanExpression <- orderedExpressionData[,notMelanColumns]
notDiffExpression <- orderedExpressionData[,!diffColumns]

diffExpression <- orderedExpressionData[,diffColumns]

###Below takes a long time to run, will need to run on server. 
#writing everything to image and writing a separate script to call DE genes called 
#callingDEGenesWilcoxon.R
save(file="expressionData.RData", list = c("neuralCrestLikeExpression", "notNCExpression", 
                                    "transitLikeExpression", "notTransExpression",
                                    "melanocyteExpression", "notMelanExpression",
                                    "diffExpression", "notDiffExpression"))
testType <- 'wilcox'
tiroshNCDE <- getDEGenes(neuralCrestLikeExpression, notNCExpression, cores=3, nGenes=30, testType='wilcox')
tiroshTransDE <- getDEGenes(transitLikeExpression, notTransExpression, cores=3, nGenes=30, testType='wilcox')
tiroshMelanDE <- getDEGenes(melanocyteExpression, notMelanExpression, cores=3, nGenes = 30, testType='wilcox')
tiroshDiffDE <- getDEGenes(diffExpression, notDiffExpression, cores=3, nGenes = 30, testType='wilcox')

cutoff <- 0.005
tiroshNCSigDE <- getSigDEGenes(tiroshNCDE, cutoff, 0)
tiroshTransSigDE <- getSigDEGenes(tiroshTransDE, cutoff, 0)
tiroshMelanSigDE <- getSigDEGenes(tiroshMelanDE, cutoff, 0)
tiroshDiffSigDE <- getSigDEGenes(tiroshDiffDE, cutoff, 0)

###Putting rest of Submap call of mappings into mapTsoitoTirosh.R
# So after exploring the discovered trajectory further, it did not make sense and didnt map properly
# Ive decided its not good. So i then tried another supervised method called ouija, and from this determined genes
# important for switching. However, ouija assumes a linear differentiation with no branching, and so there was 
# some inconsistency whereby it appeared that AXL levels suddenly dropped off toward the AXL side, which corresponded
# to increased SOD3 levels. Further, i did this without mean adjusting for tumors, and then with mean adjusting.
# By mean adjusting, it actually removed the mitf-axl axis, and so I think there is a big effect of the tumor it comes 
# from. Given this, Im retrying with monoly below using the ouija marker genes to order the cells and not correcting
# for tumor effects
##Exploring structure of trajectory further##
semisup_clustering_genes <- as.character(read.table('cellLineDEGenes.txt', header=FALSE)[,1])
HSMM3 <- setOrderingFilter(HSMM3, ouija_markers)
HSMMExplore <- reduceDimension(HSMM3, max_components = 4, #num_dim=5,
                               reduction_method = 'DDRTree', 
                               #residualModelFormulaStr = 
                               #   '~tumorTypes',
                               verbose = T,
                               auto_param_selection = TRUE,
                               #ncenter=1172,
)
#ncenter=1000)
HSMMExplore_O <- orderCells(HSMMExplore)
HSMMExplore_O$'ouijaTime' <- ouijaTime
HSMMExplore_O <- orderCells(HSMMExplore_O, root_state=13)

plot_cell_trajectory(HSMMExplore_O, color_by = "Pseudotime",
                     markers=c('MITF','AXL','APOE','SOD3'), use_color_gradient = TRUE)

plot_cell_trajectory(HSMMExplore_O, color_by = "State")

# Results now show the AXL-MITF axis that is widely accepted in the literature. Whats more, appears to show
# another axis in response to hypoxia whereby instead of becoming dedifferentiated, the cells appear to start expressing
# SOD3, which is an antioxidant. This seems to make the cells not differentiated and pick up expression of MITF. So im
# thinking this shows two different responses to hypoxia, dedifferentiation or over expression of SOD3. Going to call 
# DE genes as a function of pseudo time to see if this hypothesis holds. Furthermore, notably APOE appears to increase in 
# expression as the cells are differentiation to the the AXL high state. Reading about this, APOE has been shown to play
# a role in promoter cancer metastasis in small-cell lung cancer. Aaron thinks the metastatic programs might be similar between
# these cancer types, perhaps this is why?
branchGenesExplore <- branchTest(HSMMExplore_O, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
                          reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
                          branch_point = 1)
sigBranchGenes <- subset(branchGenesExplore, qval < 0.05)
sigBranchGenesOrdered <- sigBranchGenes[order(sigBranchGenes$qval),]
sigGenesExplore <- rownames(sigBranchGenesOrdered)

#getting the genes that tsoi used for ranking..
tsoiRankedGenes <- as.character(read.table('../Tsoi/hiearchialClustering/tsoiRankedGenesByDiffStatus.txt', header=FALSE)[,1]) #tsoi genes ranked by diff status

for (i in 1:length(tsoiRankedGenes)){ #making sure we have no genes 
  if (!tsoiRankedGenes[i]%in%rownames(HSMMExplore_O)){
    print(tsoiRankedGenes[i])
    print(i)
  }
}

#below takes forever to run, going to deploy to the server
tsoiGeneAnnotations <- character()
for (i in 1:length(tsoiRankedGenes)){
  if (i<=match("ACSL5",tsoiRankedGenes)){
    tsoiGeneAnnotations[i] <- "undiff"
  } else if (i<=match("ARL4C", tsoiRankedGenes)){
    tsoiGeneAnnotations[i] <- "undiffNC"
  } else if (i<=match("AIM2", tsoiRankedGenes)){
    tsoiGeneAnnotations[i] <- "NC"
  } else if (i<= match("ALDH1A3",tsoiRankedGenes)){
    tsoiGeneAnnotations[i] <- "NCTrans"
  } else if (i<=match("ALDH1A1", tsoiRankedGenes)){
    tsoiGeneAnnotations[i] <- "Trans"
  } else if (i<=match("HTR2B", tsoiRankedGenes)){
    tsoiGeneAnnotations[i] <- "TransMelan"
  } else if (i<=match("ADAM23", tsoiRankedGenes)){
    tsoiGeneAnnotations[i] <- "Melan"
  }
}
tsoiGeneAnnotations <- as.data.frame(tsoiGeneAnnotations, row.names = tsoiRankedGenes)
colnames(tsoiGeneAnnotations) <- c(" ")

source('/media/WorkingSpace/Brad/Tirosh/monocle_plottingFunctions.R')
heatMapResultsExplore <- plot_genes_branched_heatmapAlt(HSMMExplore_O[sigGenesExplore,],
                                              branch_point = 1,
                                              #num_clusters = 10,
                                              #cores = 30,
                                              use_gene_short_name = T,
                                              show_rownames = F,
                                              add_annotation_row = tsoiGeneAnnotations,
                                              branch_labels = c("SOD3 axis","AXL axis"),
                                              return_heatmap = TRUE
                                              )

###From here on going to characterise these clustered gene states..#### 
trajClusters <- heatMapResultsExplore$annotation_row

SOD3UpGenes <- sigGenesExplore[trajClusters==2]
SOD3UpGenes <- SOD3UpGenes[order(SOD3UpGenes)]

###Going to use the gene signatures provided by the glutamine starvation signature to see
###if they correlate with any of the axis of melanoma differentiation
require(gdata) #provides a reader for xlsx files
geneSignatures <- read.xls('glutamineSSGenes.xlsx', sheet=1, header=TRUE)

GSS <- as.character(geneSignatures$GLUTAMINE.STARVATION.SIGNATURE)
hoekProlif <- as.character(geneSignatures$HOEK_PROLIFERATIVE)
hoekInvasive <- as.character(geneSignatures$HOEK_INVASIVE)
verfProlif <- as.character(geneSignatures$VERFAILLIE_PROLIFERATIVE)
verfInvasive <- as.character(geneSignatures$VERFAILLIE_INVASIVE)
tnf <- as.character(geneSignatures$RIESENBERG_MELANOMA_TNF_RESPONSE)
axl <- as.character(geneSignatures$AXL.PROGRAM)
salS <- as.character(geneSignatures$SALUBRINAL.SIGNATURE)

source('DEFunctions.R')
geChangeObjects <- getGEChangeObjects(HSMMExplore_O, cores=15)

SOD3_exprs <- geChangeObjects$BranchA_exprs
AXL_exprs <- geChangeObjects$BranchB_exprs

plotGSChangesWithPseudotime(AXL_exprs, GSS, xTitle="AXL axis Pseudotime", 
                            yTitle="GSS", title = "AXL pseudotime with starvation signature")
plotGSChangesWithPseudotime(SOD3_exprs, GSS, xTitle="SOD3 axis Pseudotime", 
                            yTitle="GSS", title = "SOD3 pseudotime with starvation signature")
plotGSChangesWithPseudotime(SOD3_exprs, axl, xTitle="SOD3 axis Pseudotime", 
                            yTitle="AXLS", title = "SOD3 pseudotime with tirosh AXL signature")
plotGSChangesWithPseudotime(AXL_exprs, axl, xTitle="AXL axis Pseudotime", 
                            yTitle="AXLS", title = "AXL pseudotime with tirosh AXL signature")
plotGSChangesWithPseudotime(AXL_exprs, tnf, xTitle="AXL axis Pseudotime", 
                            yTitle="TNFS", title = "AXL pseudotime with TNF signature")
plotGSChangesWithPseudotime(SOD3_exprs, tnf, xTitle="SOD3 axis Pseudotime", 
                            yTitle="TNFS", title = "SOD3 pseudotime with TNF signature")
plotGSChangesWithPseudotime(AXL_exprs, hoekProlif, xTitle="AXL axis Pseudotime", 
                            yTitle="hoekProlifS", title = "AXL pseudotime with hoekProlif signature")
plotGSChangesWithPseudotime(SOD3_exprs, hoekProlif, xTitle="SOD3 axis Pseudotime", 
                            yTitle="hoekProlifS", title = "SOD3 pseudotime with hoekProlif signature")
plotGSChangesWithPseudotime(SOD3_exprs, hoekInvasive, xTitle="SOD3 axis Pseudotime", 
                            yTitle="hoekInvasiveS", title = "SOD3 pseudotime with hoekInvasive signature")
plotGSChangesWithPseudotime(AXL_exprs, hoekInvasive, xTitle="AXL axis Pseudotime", 
                            yTitle="hoekInvasiveS", title = "AXL pseudotime with hoekInvasive signature")
plotGSChangesWithPseudotime(AXL_exprs, verfProlif, xTitle="AXL axis Pseudotime", 
                            yTitle="verfProlifS", title = "AXL pseudotime with verfProlif signature")
plotGSChangesWithPseudotime(AXL_exprs, verfInvasive, xTitle="AXL axis Pseudotime", 
                            yTitle="verfInvasiveS", title = "AXL pseudotime with verfInvasive signature")
plotGSChangesWithPseudotime(SOD3_exprs, verfProlif, xTitle="SOD3 axis Pseudotime", 
                            yTitle="verfProlifS", title = "SOD3 pseudotime with verfProlif signature")
plotGSChangesWithPseudotime(SOD3_exprs, verfInvasive, xTitle="SOD3 axis Pseudotime", 
                            yTitle="verfInvasiveS", title = "SOD3 pseudotime with verfInvasive signature")

#really unsure of this SOD3 group. Based on the tnf and axl signatures, seems that cells going along this
#trajectory really belong to the AXL-axis. It makes me wonder if the gating strategy of 
#Tirosh might have resulted in some non-melanoma cells being let through, or something like that.
#Anyway, to resolve this, I'm going to find DE genes in the SOD3 high cells compared to everything 
#else. Then rank the genes by p-value for being upregulated. Then input this gene list into
#GSEA. That should give me a better indication what this group actuall represents. 

#state 1 includes the SOD3 high axis
state1Cells <- which(HSMMExplore_O$State==1)
SOD3Cells <- state1Cells[HSMMExplore_O$Pseudotime[state1Cells]>=3.05]
notSOD3Cells <- setdiff(1:1172, SOD3Cells)

sodStatus <- character(1172)
sodStatus[SOD3Cells] <- "SOD3"
sodStatus[notSOD3Cells] <- "notSOD3"
HSMMExplore_O$'sodStatus' <- sodStatus
plot_cell_trajectory(HSMMExplore_O, color_by = "sodStatus",markers=c('SOD3'), 
                     use_color_gradient = TRUE) +
  facet_wrap(~sodStatus, nrow = 1)

sodExpr <- exprs(HSMMExplore_O)[,SOD3Cells]

otherCells <- which(HSMMExplore_O$State!=1)
notSodExpr <- exprs(HSMMExplore_O)[,otherCells]

sodDE <- getDEGenes(sodExpr, notSodExpr, nGenes=8950, cores=30, testType="wilcox")
sodDEGenes <- getSigDEGenes(sodDE, 1, 0)

write.table(sodDEGenes, file='SOD3SubtypeGenesRanked.rnk',quote=FALSE,row.names = TRUE, sep="\t")

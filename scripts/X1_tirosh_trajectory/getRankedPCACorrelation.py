#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 14:13:56 2018

@author: brad
"""

import os, pandas, numpy, math, scipy, sklearn, matplotlib, sys
from scipy import stats
matplotlib.use('Agg')
from matplotlib import pyplot
from sklearn.cross_validation import train_test_split
from sklearn.cluster import DBSCAN
from sklearn import metrics
#os.chdir("/home/brad/Desktop/Uni_Studies/Honours/Tirosh")

# <codecell>

tiroshDataFile = open("GSE72056_melanoma_single_cell_revised_v2.txt", 'r')

dontNormalise = ['"non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)"',
                 '"malignant(1=no,2=yes,0=unresolved)"',"Cell","tumor"]
#Loading in my data
tiroshData = {}
nLines = 0
for line in tiroshDataFile:
    parts = line.split('\t')
    lineName = parts[0]
    
    parts.remove(parts[0])
    tiroshData[lineName] = parts
    
    if nLines == 15: #Only first 15 lines since is alot of data
        break
    #nLines+=1
    
tiroshDataFile.close()

# <codecell>

#Casting to colData to appropriate types
tiroshDataFrame = pandas.DataFrame.from_dict(tiroshData)
for column in tiroshDataFrame:
    tiroshDataFrame[column] = pandas.to_numeric(tiroshDataFrame[column], errors='ignore')
    
# <codecell>
#Pre-processing the data


geneNames = []

for columnName in tiroshDataFrame:
    if columnName not in dontNormalise:
        geneNames.append(columnName)

#def convertToLogTPM(weirdNormalisedTPM):
#    newTPM = 2**weirdNormalisedTPM
#    newTPM = (newTPM-1)*10
#    newTPM = math.log2(newTPM+1)
#    return newTPM


#for geneName in geneNames:
#    tiroshDataFrame[geneName] = numpy.log2((((2**tiroshDataFrame[geneName])-1)*10)+1)


# <codecell>
#Checking approximately normal
from sklearn import preprocessing

min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0,1))
#data = [tiroshDataFrame[geneNames[0]],tiroshDataFrame[geneNames[1]]]
#minMaxTransformed = min_max_scaler.fit_transform(data)
#ax = pyplot.subplot(221)
#ax2 = pyplot.subplot(221)
#stats.probplot(data[0], plot=ax)
#stats.probplot(minMaxTransformed[0], plot=ax2)
#ax.show()
#ax2.show()
##Just using the DE genes from verfaille
#print(len(tiroshDataFrame))
#for feature in tiroshDataFrame:
#    if feature not in verfDEGenes and feature not in dontNormalise:
#        tiroshDataFrame.drop([feature],axis=1,inplace=True)
#        try:
#            geneNames.remove(feature)
#        except:
#            print(feature)
#print(len(tiroshDataFrame))

#Just pulling out the malignant melanocytes
malignant = tiroshDataFrame[dontNormalise[1]]==2
tiroshMalignantDF = tiroshDataFrame[malignant]

geneNames.append('Cell')
featureValues = tiroshMalignantDF[geneNames]
cellToID = pandas.factorize(featureValues['Cell'])
featureValues['Cell'] = cellToID[0]
tiroshScaled = min_max_scaler.fit_transform(featureValues)

# <codecell>
def placeItemRanked(listLists, p_value, geneName):
    sortList = listLists[0]
    geneList = listLists[1]
    if len(sortList)==0:
        sortList.append(p_value)
        geneList.append(geneName)
        return
    
    for i in range(0,len(sortList)):
        if sortList[i]<p_value and i!=len(sortList)-1:
            continue
        elif sortList[i]>=p_value:
            sortList.insert(i, p_value)
            geneList.insert(i, geneName)
            return
        else:
            sortList.append(p_value)
            geneList.append(geneName)
            
    return
    #print(sortList)
    #print(p_value, geneName)
        
# <codecell>
#loading in the PCA coords
pcaCoords = pandas.read_csv("pcaCoords_allGenes.txt",sep="\t")
pcaCoordsMatrix = pandas.DataFrame.as_matrix(pcaCoords.iloc[:,1:pcaCoords.shape[1]])

pcaCoordsMatrix = numpy.transpose(pcaCoordsMatrix)
tiroshScaled = numpy.transpose(tiroshScaled)

# <codecell>
resultsDict = {}
for PCi in range(0,len(pcaCoordsMatrix)):
    resultsDict["PC"+str(PCi+1)] = [[], []]
    for genei in range(0,len(tiroshScaled)):
        geneName = geneNames[genei]
        if (geneName=="Cell"):
            continue
        
        results = scipy.stats.spearmanr(list(tiroshScaled[genei]), list(pcaCoordsMatrix[PCi]))
        p_value = results[1]
        placeItemRanked(resultsDict["PC"+str(PCi+1)], p_value, geneName)
        
###don't need to correct for multiple hypothesis testing since just want to rank
        
# <codecell>
#writing results to file
results = open("genesRankedWithPC.txt",'w')
for i in range(len(geneNames)-1):
    line = ''
    for PCi in range(len(resultsDict)-1):
        info = resultsDict["PC"+str(PCi+1)]
        line += info[1][i] + "\t" + str(info[0][i]) + "\t"
        
    info = resultsDict["PC"+str(PCi+2)]
    line += info[1][i] + "\t" + str(info[0][i]) + "\n"
    #print(line)
    results.write(line)
        
results.close()
        

    
        










tsoiData = read.table('../Tsoi/data/tsoiLog2TPM.txt')
tiroshData = read.table('tiroshMelanomaData_ordered.txt')

tsoiGenes <- rownames(tsoiData)
tiroshGenes <- rownames(tiroshData)[which(rownames(tiroshData)!='tumor')]

tsoiToTirosh <- match(tiroshGenes, tsoiGenes)
tsoiToTirosh <- tsoiToTirosh[!is.na(tsoiToTirosh)]

tumor <- as.character(tiroshData['tumor',])

overlapGenes <- tsoiGenes[tsoiToTirosh]

tsoiData <- tsoiData[overlapGenes,]

tiroshData <- tiroshData[overlapGenes,]
tiroshData <- rbind(tumor, tiroshData)
rownames(tiroshData)[1] <- 'tumor'

####savingResults####
write.table(tsoiData, file='tsoiTPMwTirosh.txt',quote=FALSE)
write.table(tiroshData, file='tiroshTPMwTsoi.txt',quote=FALSE)

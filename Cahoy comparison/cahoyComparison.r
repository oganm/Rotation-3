#Functions
#################################
dpaste = function (...){
  paste(..., sep='')
}

des = function(sb){
  return(sb[, designR])
}

#return the data of a subset
dat = function(sb){
  return(sb[, dataR])
}

foldChange = function (group1, group2){
  fold = group1 / group2
  return(which(fold <= (1/10) | fold >= 10))
}

#############################
library(affy)
require(gdata)
library(mouse4302.db)
mainDir = 'C:/Users/Ogan/Dropbox/Rotation 3/'
transcriptomeSources = read.xls(dpaste(mainDir,'Data/journal.pone.0016493.s010.xls'),perl='C:/Perl64/bin/perl.exe')
transcriptomeSources$SampleNo=1:64
colnames(transcriptomeSources)=c('officialName','cellType','anatomical','age','method','isolation','amplification','RNAamount','platform','reference','sampleNo')

usedData = (read.csv(dpaste(mainDir, 'Data/usedData'), header = F))
usedData = 2^usedData

geneNames = as.character(read.table(dpaste(mainDir, 'Data/usedGenes'))$V1)
  
cahoyData = read.xls('335_Cahoy_S_Table_S3b_dChip_Data_Table_2007-09-11.xls')
cahoyMatrix = matrix(c(
           as.double(cahoyData$Astrocytes..P7.P8), 
           as.double(cahoyData$Astrocytes.P17),
           as.double(cahoyData$Astrocytes.Gray.P17),
           as.double(cahoyData$Neurons.P7),
           as.double(cahoyData$Neurons.P16),
           as.double(cahoyData$Myelin.Oligos),
           as.double(cahoyData$Oligos),
           as.double(cahoyData$OPCs),
           1:20932),
           nrow = 9, byrow = T
           )
cahoyMatrix[1:8, ] =2^cahoyMatrix[1:8, ] 
cahoyNames = cahoyData$Gene.Symbol

cahoySelect = cahoyMatrix[,which(tolower(cahoyNames) %in% tolower(geneNames)) ]

cahoySelectNames = cahoyNames[ cahoySelect[9, ]]

#some gene symbols did not match because of alternate names. so I am looking at the ranking of matching ones


#now to take Cahoy specific data
usedDataD = cbind(usedData, transcriptomeSources)
dataR=1:1872
designR=1874:1883
usedDataCD = usedData[usedDataD$reference %in% 'Cahoy et al., 2008',]
usedDataCD= rbind(usedDataCD,1:2131)
usedDataCCD = usedDataCD[ ,tolower(geneNames) %in% tolower(cahoySelectNames)]



#repeating previous comparisons
##################################
astro = usedDataCCD[c(2, 3),dataR ]
cahoyAstro = cahoySelect[c(2, 3), ]

mixed = usedDataCCD[c(5),dataR ]
cahoyMixed = cahoySelect[c(5),]

oligo = usedDataCCD[6, dataR]
cahoyOligo = cahoySelect[6,]

#neurons average
mixedAv = mixed
mixedCahoyAv = cahoyMixed

#glial average
glialSet = rbind(astro, oligo)
glialAverage = apply(glialSet, 2, mean)

glialCahoySet = rbind(cahoyAstro, cahoyOligo)
glialCahoyAverage = apply(glialCahoySet, 2, mean)


#proportions

require('VennDiagram')
setHits=geneNames[usedDataCCD[9,foldChange(glialAverage, mixedAv)]]

cahoyHits=unique(cahoyNames[cahoySelect[9, foldChange(glialCahoyAverage, mixedCahoyAv)]])

length(setHits)

length(cahoyHits)

length(intersect(setHits,cahoyHits))
aList=list(setHits=tolower(as.character(setHits)), cahoyHits=unique(tolower(as.character(cahoyHits))))
plot.new()
venn.plot=venn.diagram(aList,filename=NULL,fill=c('red','blue'))
grid.draw(venn.plot)




which(cahoyNames %in% 'Eno2')


cahoyData$Gene.Symbol[which(tolower(cahoyData$Gene.Symbol) %in% tolower(geneNames))] %in% 'Eno2'



#checking if the gene names match the names from the chip... spoiler.. they do. 
affy.data = ReadAffy()
eset.mas5 = mas5(affy.data)
exprSet.nologs = exprs(eset.mas5)

probe.IDs = rownames(exprSet.nologs)
allGenes = unlist(mget(probe.IDs, mouse4302SYMBOL))

geneNames[!tolower(geneNames) %in% tolower(allGenes)]


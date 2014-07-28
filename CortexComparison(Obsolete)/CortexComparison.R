#give names for known genes
giveName = function(numList){
  return(geneNames[which(geneIndex %in% numList)])
}

#function for pasting directories without 'sep' variable
dpaste = function (...){
  paste(..., sep='')
}

#view the design part of a subset of data
des = function(sb){
  return(sb[, designR])
}

#return the data of a subset
dat = function(sb){
  return(sb[, dataR])
}

#look for 10 fold difference
foldChange = function (group1, group2){
  fold = group1 / group2
  return(which(fold <= (1/5) | fold >= 5))
}


require(gdata)
mainDir = 'C:/Users/Ogan/Dropbox/Rotation 3/'
#open all data design file and tweak it a little
#######
transcriptomeSources = read.xls(dpaste(mainDir,'journal.pone.0016493.s010.xls'),perl='C:/Perl64/bin/perl.exe')
transcriptomeSources$SampleNo=1:64
colnames(transcriptomeSources)=c('officialName','cellType','anatomical','age','method','isolation','amplification','RNAamount','platform','reference','sampleNo')

transcriptomeSources$age = gsub('~', '', transcriptomeSources$age)
transcriptomeSources$age = gsub('P', '', transcriptomeSources$age)
transcriptomeSources$age = gsub('7-8', '7.5', transcriptomeSources$age)
transcriptomeSources$age[grepl('(precise age not given)',transcriptomeSources$age)] = 60
transcriptomeSources$age = as.numeric(transcriptomeSources$age)

transcriptomeSources$anatomical = gsub('Corpus Striatum', 'Striatum', transcriptomeSources$anatomical)
transcriptomeSources$anatomical = gsub('Motor Cortex', 'Motor cortex', transcriptomeSources$anatomical)
transcriptomeSources$anatomical = gsub('Layer 5A Cortex', 'Cortex', transcriptomeSources$anatomical)
transcriptomeSources$anatomical = gsub('Layer 5B Cortex', 'Cortex', transcriptomeSources$anatomical)
transcriptomeSources$anatomical = gsub('Layer 6 Cortex', 'Cortex', transcriptomeSources$anatomical)


transcriptomeSources$cellType[grepl('^Neurons', transcriptomeSources$cellType)] = 'Mixed Neurons'
transcriptomeSources$cellType=gsub(', P.*?$','',transcriptomeSources$cellType)
#######

#open gene name files
geneNames = as.character(read.table(dpaste(mainDir, 'usedGenes'))$V1)
geneIndex = as.numeric(read.csv(dpaste(mainDir, 'colsToUse'), header=FALSE))


#data preparation
##############
allDataPre = (read.csv(dpaste(mainDir, 'Data/allData'), header = F))
allData=2^allDataPre
rownames(allData)=dpaste('V',1:64)

allDataD = cbind(allData, transcriptomeSources)
dataR=1:14580
designR=14582:14591
#####################

#acquire relevant data
################################################
#select cortex and above 14 days
cortex = allDataD[allDataD$anatomical == 'Cortex' & allDataD$age >= 14, ]

#cortex oligodendrocytes
cortexOligo = cortex[grepl('Oligodendrocytes', cortex$cellType), ]

#cortex astrocytes
cortexAstro = cortex[cortex$cellType %in% 'Astrocytes', ]

#cortex pyramidal
cortexPyramidal = cortex[grepl('Pyramidal',cortex$cellType), ]

#interneurons
cortexInter = cortex[cortex$cellType %in% 'Interneurons',]

#mix
cortexMix = cortex[cortex$cellType %in% 'Mixed Neurons',]
###############
#cortexOligo
#cortexAstro
#cortexPyramidal
#cortexInter
#cortexMix


#neuron vs glial
##################

#neuron average
neuronSet = rbind(cortexPyramidal, cortexInter, cortexMix)
neuronAverage = apply(dat(neuronSet), 2, mean)

#glial average
glialSet = rbind(cortexAstro, cortexOligo)
glialAverage = apply(dat(glialSet), 2, mean)

neuronVSglial = foldChange(neuronAverage, glialAverage)


###############
#neuronVSglial


#pyramidal vs interneuron
################
interAverage = dat(cortexInter)

pyramidalAverage = apply(dat(cortexPyramidal), 2, mean)

pyramidalVSinter = foldChange(interAverage, pyramidalAverage)
######################
#pyramidalVSinter


#oligo vs astro
#############
oligoAverage = apply(dat(cortexOligo), 2, mean)

astroAverage = apply(dat(cortexAstro), 2, mean)

oligoVSAstro = foldChange(oligoAverage, astroAverage)
##################
#oligoVSAstro


giveName(oligoVSAstro[oligoVSAstro %in% neuronVSglial])

which(pyramidalVSinter %in% neuronVSglial)

giveName(neuronVSglial)
giveName(pyramidalVSinter)
giveName(oligoVSAstro)



giveName = function(numList){
  return(geneNames[which(geneIndex %in% neuronVSglial)])
}


dataCat = data.frame(exp = c(cortexOligo[,301], cortexAstro[,301],cortexPyramidal[,301], cortexInter[,301], cortexMix[, 301]), cellType = c(
                     rep('oligo',nrow(cortexOligo)),rep('astro',nrow(cortexAstro)),rep('pyramidal',nrow(cortexPyramidal)),rep('inter',nrow(cortexInter)), rep('mix',nrow(cortexMix))))
                     
ggplot(dataCAt)

cortexAstro
#cortexPyramidal
#cortexInter
#cortexMix
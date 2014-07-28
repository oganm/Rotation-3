
# initialize with ogbox and matlab ---------
source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')

insist(R.matlab)
insist(ggplot2)

Matlab$startServer()
matlab = Matlab()
isOpen = open(matlab)

predCon = function(samples, types){
  selectGroupData = types
  if (!isOpen)
    throw("MATLAB server is not running")
  
  random = samples
  selectGroupData = types
  random = matrix(random)
  predictions = matrix(0, nrow = ncol(selectGroupData), ncol = ncol(random))
  
  MTakenGeneByType = as.matrix(selectGroupData)
  setVariable(matlab, MTakenGeneByType=MTakenGeneByType)
  
  diagNegTypes = diag(x = -1, nrow = ncol(selectGroupData), ncol = ncol(selectGroupData))
  refZeroTypes = matrix(0, ncol(selectGroupData), ncol = 1)
  setVariable(matlab, refZeroTypes=refZeroTypes)
  setVariable(matlab, diagNegTypes=diagNegTypes)
  
  
  for (i in 1:ncol(random)){
    voxelData = as.matrix(random[, i])
    setVariable(matlab, voxelData=voxelData)
    evaluate(matlab, 'fitVoxelsToTypesBis= lsqlin( MTakenGeneByType, voxelData, diagNegTypes, refZeroTypes );')
    predictions [, i] =  getVariable(matlab, 'fitVoxelsToTypesBis')[[1]]
    
  } 
  
  predictions[predictions <10e-15] = 0
  return(predictions)
}

giveGI = function(daGene){
  return(match(daGene, geneData$Gene.Symbol))  
}

#group replicates together by... mean?
groupOut = function(selectExprData, selectDes){
  selectGroupData = matrix(0, nrow = nrow(selectExprData), ncol = length(unique(selectDes$ourNaming)))
  groupNames = unique(selectDes$ourNaming)
  colnames(selectGroupData) = groupNames
  for (i in 1:length(groupNames)){
    groupIndex = selectDes$ourNaming == groupNames[i]
    groupAverage =  tryCatch({apply(selectExprData[, groupIndex ], 1, mean)},
                             error = function(cond){
                               print('fuu')
                               return(selectExprData[, groupIndex ])
                             })
    selectGroupData[, i] = groupAverage 
  }  
  return(selectGroupData)
}


# As always open the bloody files. ourNaming field added to the design file ----------------
design = read.table('C:/Users/Ogan/Dropbox/Rotation 3/Data/normalizedDesign.csv',header=T,sep='\t')
design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)
allDataPre = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/Data/mostVariableQuantileNormalized', header = T)
geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData = exprData
design = design[match(colnames(exprData),design$sampleName,),]
rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<5)

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]



# Get the list and indexes of final selected genes. Also create the selected matrix --------------------
filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/cheatsyClust/Results/cheatsy2Fix",include.dirs = FALSE)
fileContents = lapply(paste('C:/Users/Ogan/Dropbox/Rotation 3/cheatsyClust/Results/cheatsy2Fix/', filenames, sep = ''), read.table)
geneList = vector(mode = 'character', length = 68000)
step = 1
for (i in fileContents){
  length(i$V1)
  geneList[step:(step + length(i$V1)-1)] = as.character(i$V1)
  step = step+length(i$V1)
}

geneList = trimElement(geneList,'')


length(unique(geneList))
geneList = unique(geneList)
#geneList = c('Mobp', 'Slc1a3', 'Slc1a2', 'Aqp4', 'Mag', 'Mog', 'Tac1', 'Neurod6')
geneIndexes = trimNAs(giveGI(geneList)) #note that the data is already trimmed for low expressed ones so indexes arent from original exprData

# also remove the non desirables
selectGeneData = geneData[geneIndexes, ]
#old criteria
#selectExprData = exprData[geneIndexes, (design$age>=14)&(!design$cellType2=='Astroglia')&(!design$originalIndex %in% c(6,7,8))]
selectExprData = 2^exprData[geneIndexes, !(is.na(design$ourNaming) | design$ourNaming=='Mix')]
selectDes = design[!is.na(design$ourNaming),]


# purity ---------------


purity = matrix(0, nrow = length(unique(selectDes$ourNaming)), ncol = ncol(selectExprData))


#remove mixed neurons
backupSelectExprData = selectExprData


selectExprData = backupSelectExprData
newSelectExprData = selectExprData

n = 1
purityChange = vector(mode = 'list', length = n)
exprDataChange = vector(mode = 'list', length = n)
for (j in 1:n){
  for (i in 1:ncol(selectExprData)){
    samples = selectExprData[, i]
    selectGroupData  = groupOut(selectExprData[, -i], selectDes[-i, ])
    singlePur = predCon(samples, selectGroupData)
    purity [, i]= singlePur
    
    groupDeSample =  which(colnames(selectGroupData)  %in% selectDes$ourNaming[i])
    
    toRemove =apply( 
      t(repCol(singlePur[-groupDeSample],nrow(selectGroupData))) * selectGroupData[, groupDeSample],
      1,
      sum)
    oldSum = sum(selectExprData[, i ])
    
    
    newSelectExprData[,i] = selectExprData[, i ]- toRemove
    newSelectExprData[,i] = newSelectExprData[,i] *oldSum/sum(newSelectExprData[,i])
  }
  purityChange[[j]] = purity
  exprDataChange[[j]] = newSelectExprData
  selectExprData = newSelectExprData
  
}

purityChange[[1]]
rownames(purityChange[[1]]) = colnames(selectGroupData)
colnames(purityChange[[1]]) = make.names(selectDes$cellType2, unique = T)


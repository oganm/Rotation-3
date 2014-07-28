#last execution: age 14<= with undesirables discluded

foldChange = function (group1, group2, f = 10){
  fold = group1 / group2
  chosen =  which(fold <= (1/f) | fold >= f)
  return(
    data.frame(index = chosen, foldChange = fold[chosen])
  )
}

foldChange2 = function (group1, group2, f = 10){
  
  group1 = unlist(group1)
  group2 = unlist(group2)
  groupAverage1 = tryCatch({apply(exprData[group1,], 2, median)},
                           error = function(cond){
                             print('fuu')
                             return(exprData[i,])
                           })
  
  
  
  groupAverage2 = tryCatch({apply(exprData[group2,], 2, median)},
                           error = function(cond){
                             print('fuu')
                             return(exprData[i,])
                           })
  
  g19 = groupAverage1 < 9.5 & groupAverage1 > 7
  g16 = groupAverage1 <6
  g29 = groupAverage2 < 9.5 & groupAverage2>7
  g26 = groupAverage2 < 6
  
  groupAverage2[g26 & g19] = (apply(exprData[group2,g26 & g19], 2, max))
  groupAverage1[g16 & g29] = (apply(exprData[group2,g16 & g29], 2, max))
  
  groupAverage2[g16 & g29] = (apply(exprData[group2,g16 & g29], 2, min))
  groupAverage1[g26 & g19] = (apply(exprData[group2,g26 & g19], 2, min))
  
  
  group1 = groupAverage1
  group2 = groupAverage2
  
  
  
  g1n = group1<9.5 & group1>7
  g1s = group1<6
  g2n = group2<9.5 & group2>7
  g2s = group2<6
  fold = (groupAverage1 - groupAverage2)
  chosen =  which(fold >= (log(f)/log(2)) | fold <= log(1/f)/log(2) | (g1n & g2s) | (g2n & g1s))
  return(
    data.frame(index = chosen, foldChange = fold[chosen])
  )
}


#Data preperation
####################
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

#disclusion is done below when comparing it to rest.
#discludeSamples = which(design$age<=14)
#exprData = exprData[, -discludeSamples]
#design = design[-discludeSamples, ]
exprData = t(exprData)
#########################################



comparisons =
  "all"

intersections = ''
foldTreshold = 10
#analysis
############################################

realGroups = list(toddler = which(design$age <=7.5 & design$cellKind == 'Neuron'),
                  youngling = which(design$age >7.5 & design$age <=17 & design$cellKind == 'Neuron'),
                  youngAdult = which(design$age > 17 & design$age <= 45 & design$cellKind == 'Neuron'),
                  adult = which(design$age >45 & design$age <=60 & design$cellKind == 'Neuron')
  
  )

conn = textConnection(comparisons)
comp = read.csv(conn, header = F, strip.white = T)

groupNames = names(realGroups)

#uses medians
groupAverages = list()
for (i in realGroups){

  
  groupAverage = i
  
  groupAverages = c(groupAverages, list(groupAverage))
}
names(groupAverages)= groupNames

combinations = combn(1:length(realGroups),2)


for (i in 1:ncol(combinations)){
  compare = combinations[ , i]
  #exit if comparison not requested
  if (!comp[1,1]=='all'){
    daTruth = (comp == groupNames[compare[1]] | comp == groupNames[compare[2]])
    if (!any(apply(daTruth,1,all))){
      next
    }
  }
  
  fileName = paste('Results/',groupNames[compare[1]],'-', groupNames[compare[2]], sep='')
  fChange = foldChange2(groupAverages[[compare[1]]], groupAverages[[compare[2]]], foldTreshold)
  fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fChange$index], geneFoldChange= fChange$foldChange )
  fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]
  print(fileName)
  print(i)
  write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
}

#broken
if ('rest' %in% comp$V1){
  for (i in 1:length(realGroups)){
    fileName = paste('Rest/',names(realGroups)[i],sep ='')
    i = realGroups[[i]]

          
    groupAverage = i
    
    restAverage = 
    fChange = foldChange2(groupAverage, restAverage, foldTreshold)
    fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fChange$index], geneFoldChange= fChange$foldChange )
    fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]
    print(fileName)
    print(i)
    write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
    
  }
}



require(tuneR)
require(cluster)
require(ggplot2)
t = seq(0, 3, 1/8000) #times in seconds if sample for 3 seconds at 8000Hz
u = (2^15-1)*sin(2*pi*440*t) #440 Hz sine wave that lasts t length seconds (here, 3 seconds) 
w = Wave(u, samp.rate = 8000, bit=16) #make the wave variable 

listDepth = function(deList){
  step = 1
  while (T){
    if (typeof(eval( parse(text = paste(c("deList",rep('[[1]]',step)),sep='',collapse = '')))) != "list"){
      return(step)
    }
    step = step +1
  }
}





dpaste = function (...){
  paste(..., sep='')
}

getParent = function(step = 1){
  wd = getwd()
  for (i in 1:step){
    setwd('..')
  }
  parent = getwd()
  setwd(wd)
  return(paste(parent,'/',sep=''))
}

giveGI = function(daGene){
  return(match(daGene, geneData$Gene.Symbol))  
}

#last execution: age 14<= with undesirables discluded

foldChange = function (group1, group2, f = 10){
  g1n = group1<9.5 & group1>8
  g1s = group1<6
  g2n = group2<9.5 & group2>8
  g2s = group2<6
  fold = (group1 - group2)
  chosen =  which(fold >= (log(f)/log(2)) | fold <= log(1/f)/log(2) | (g1n & g2s) | (g2n & g1s))
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
  
  g19 = groupAverage1 < 9.5 & groupAverage1 > 8
  g16 = groupAverage1 <6
  g29 = groupAverage2 < 9.5 & groupAverage2>8
  g26 = groupAverage2 < 6
  
  
  
  tempGroupAv1 = vector(length = length(groupAverage1))
  tempGroupAv2 = vector(length = length(groupAverage2))
  
  tempGroupAv2[g26 & g19] = tryCatch({(apply(exprData[group2, g26 & g19], 2, max))},
                                       error = function(cond){
                                         print('I hate you damn it!')
                                         return(max(exprData[group2, g26 & g19]))
                                       })
   
  tempGroupAv1[g16 & g29] = tryCatch({(apply(exprData[group1,g16 & g29], 2, max))},
                                       error = function(cond){
                                         print('I hate you damn it!')
                                         return(max(exprData[group1,g16 & g29]))
                                       })
   
  tempGroupAv2[g16 & g29] = tryCatch({(apply(exprData[group2,g16 & g29], 2, min))},
                                        error = function(cond){
                                        print('I hate you damn it!')
                                         return(min(exprData[group2,g16 & g29]))
                                       })
   
  tempGroupAv1[g26 & g19] = tryCatch({(apply(exprData[group1,g26 & g19], 2, min))},
                                       error = function(cond){
                                         print('I hate you damn it!')
                                         return(min(exprData[group1,g26 & g19]))
                                       })
  #groupAverage1[5124]
  #groupAverage2[5124]
  

  #groupAverage1[7067]
  #groupAverage2[7067]
  
  add1 = g19 & g26 & tempGroupAv1>tempGroupAv2
  add2 = g29 & g16 & tempGroupAv2>tempGroupAv1 
    

  fold = (groupAverage1 - groupAverage2)
  chosen =  which({(fold >= (log(f)/log(2))) & !(g19 & g26) } | {(fold <= log(1/f)/log(2)) &  !(g29 & g16)}| add1 | add2)
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


groups = list(neurons = list(
  allCholinergicN = list(
    cholinergicN = c(13, 11), 
    #motorCholinN = c(10,12),
    motorCholinStemN = 10,
    motorCholinSpinalN = 12
  ),
  allGlutamergic = list(
    #ubcN = 18,
    glutaN = 53,
    #pyramidalN = c(40, 43, 44,45, 46, 47, 48, 49)
    pyramidalCorticospineN = c(43, 44),
    pyramidalCingulateN = c(45, 47),
    pyramidalSomatosensoryN = c(46),
    pyramidalAmygdalaN = c(48),
    pyramidalHippo = 49
    
  ),
  allGabaergic = list(
    #spinyN = c(15, 16),
    spinyDrd1 = 15,
    spinyDrd2 = 16,
    golgiN = 17,
    purkinjeN = 25,
    granuleN = 20,
    interN = 14,
    #gabaN = c(58, 62, 63, 64, 59, 56, 55, 54, 57),
    gabaPV = c(58, 62, 63, 64),
    gabaReln = 59,
    gabaSSTReln = 57,
    gabaSSTRelnCalb = 56,
    gabaVIPReln = c(54, 55),
    basketN = 19
  ),
  
  mixedN = c(26, 34)
),

glials = list(
  oligoG = list(
    matureOligoG = c(21, 22, 35),
    mixedOligoG = c(23, 24)
  ),
  allAstrocytesG = list(
    bergmanG = 27,
    astrocytesG = c(31, 32)
  )
)

)
# translate indexes into sample nos
realGroups = vector(mode= 'list', length(groups))
names(realGroups) = names(groups)
for (i in 1:length(groups)){
  realGroups[[i]] = vector(mode = 'list', length(groups[[i]]))
  names(realGroups[[i]]) = names(groups[[i]])
  for (j in 1:length(groups[[i]])){
    realGroups[[i]][[j]] = vector(mode = 'list', length(groups[[i]][[j]]))
    names(realGroups[[i]][[j]]) = names(groups[[i]][[j]])
    for (k in 1:length(groups[[i]][[j]])){
      realGroups[[i]][[j]][[k]] = which(design$originalIndex %in% groups[[i]][[j]][[k]])
    }   
  }
}

## failure to automate with 5 hours of sleep. manual stepping it is...
group1 = realGroups[[1]][[3]][[4]]
group2 = realGroups[[1]][[3]][[10]]

group1 = unlist(realGroups[[1]])
group2 = unlist(realGroups[[2]])

fc = foldChange2(group1,group2)
fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fc$index], geneFoldChange= 2^fc$foldChange )
fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]

fileName = 'Results/neuron-glial'
write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)


tempGroups = realGroups[[1]]
foldTreshold = 10
for (tempGroups in realGroups){
  tempAverages = list()
  for (i in tempGroups){
    groupAverage = tryCatch({apply(exprData[unlist(i),], 2, median)},
                            error = function(cond){
                              print('fuu')
                              return(exprData[i,])
                            })
    
    
    tempAverages = c(tempAverages, list(groupAverage))
  }
  tempAverages = tempGroups
  names(tempAverages)= names(tempGroups)
  combinations = combn(1:length(tempGroups),2)
  for (i in 1:ncol(combinations)){
    compare = combinations [ ,i]
    fileName = paste('Results/',names(tempGroups)[compare[1]],'-', names(tempGroups)[compare[2]], sep='')
    fChange = foldChange2(tempAverages[[compare[1]]], tempAverages[[compare[2]]], foldTreshold)
    fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fChange$index], geneFoldChange= 2^fChange$foldChange )
    fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]
    print(fileName)
    print(i)
    write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
  }
  
  for (tempGroups2 in tempGroups){
    tempAverages = list()
    for (j in tempGroups2){
      groupAverage = tryCatch({apply(exprData[unlist(j),], 2, median)},
                              error = function(cond){
                                print('fuu')
                                return(exprData[j,])
                              })
      tempAverages = c(tempAverages, list(groupAverage))
    }
    tempAverages = tempGroups2
    names(tempAverages)= names(tempGroups2)
    combinations = combn(1:length(tempGroups2),2)
    for (i in 1:ncol(combinations)){
      compare = combinations [ ,i]
      fileName = paste('Results/',names(tempGroups2)[compare[1]],'-', names(tempGroups2)[compare[2]], sep='')
      if (fileName == "Results/-"){
        next
      }
      fChange = foldChange2(tempAverages[[compare[1]]], tempAverages[[compare[2]]], foldTreshold)
      fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fChange$index], geneFoldChange= 2^fChange$foldChange )
      fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]
      print(fileName)
      print(i)
      write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
    }
  }
}
#play(w)



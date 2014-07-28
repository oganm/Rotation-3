foldChange = function (group1, group2, f = 10){
  fold = group1 / group2
  chosen =  which(fold <= (1/f) | fold >= f)
  return(
    data.frame(index = chosen, foldChange = fold[chosen])
  )
}


#Data preperation
####################
design = read.table('normalizedDesign.csv',header=T,sep='\t')


design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)

allDataPre = read.csv('allNormalized', header = T)
geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData = 2^exprData
design = design[match(colnames(exprData),design$sampleName,),]
rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<2^5)

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]

#discludeSamples = which(design$age<=14)
#exprData = exprData[, -discludeSamples]
#design = design[-discludeSamples, ]
exprData = t(exprData)
#########################################


#define groups
#####################
#all adults
#mode = 'index'
mode = 'name'
#corr = cor(t(exprData))
#set.seed(3)
#a = heatmap.3(corr,  Rowv = F)
#use the original denrogram
dendro = a$colDendrogram



groups = list(b1_oligodendro = labels(dendro[[1]][[1]][[1]]), 
              b2_astrocytes = labels(dendro[[1]][[1]][[2]]),
              b2_1_oldAstro = labels(dendro[[1]][[1]][[2]][[1]]),
              b2_2_youngAstro = labels(dendro[[1]][[1]][[2]][[1]]),
              b3_oligo = labels(dendro[[1]][[2]][[1]]),
              b3_1_oligoPre = labels(dendro[[1]][[2]][[1]][[1]]),
              b3_2_oligo =labels(dendro[[1]][[2]][[1]][[2]]),
              b4_astro.bergman = labels(dendro[[1]][[2]][[2]][[1]]),
              b4_1_astro = labels(dendro[[1]][[2]][[2]][[1]][[1]]),
              b4_2_astro.bergman = labels(dendro[[1]][[2]][[2]][[1]][[2]]),
              b5_oldOligo = labels(dendro[[1]][[2]][[2]][[2]]),
              b5_1_LoneMatureOligo = labels(dendro[[1]][[2]][[2]][[2]][[1]][[1]]),
              b5_2_mixedOligo = labels(dendro[[1]][[2]][[2]][[2]][[1]][[2]]),
              b5_3_mixedOligoCortex = labels(dendro[[1]][[2]][[2]][[2]][[2]][[1]]),
              b6_purkinje = labels(dendro[[2]][[1]][[1]]),
              b7_ubc = labels(dendro[[2]][[1]][[2]][[1]]),
              b8_granule = labels(dendro[[2]][[1]][[2]][[2]][[1]][[1]]),
              b9_basket = labels(dendro[[2]][[1]][[2]][[2]][[1]][[2]][[1]]),
              b10_golgi = labels(dendro[[2]][[1]][[2]][[2]][[1]][[2]][[2]]),
              b11_cholinergic = labels(dendro[[2]][[1]][[2]][[2]][[2]]),
              b11_1_forebrainCholi = labels(dendro[[2]][[1]][[2]][[2]][[2]][[1]]),
              b11_2_striatum.forebrainCholi = labels(dendro[[2]][[1]][[2]][[2]][[2]][[2]][[1]]),
              b11_3_choli.motor = labels(dendro[[2]][[1]][[2]][[2]][[2]][[2]][[2]]),
              b12_spiny = labels(dendro[[2]][[2]][[1]][[1]][[1]]),
              b13_mixed = labels(dendro[[2]][[2]][[1]][[1]][[2]]),
              b14_pyramidal.inter = labels(dendro[[2]][[2]][[1]][[2]]),
              b15_mixed = labels(dendro[[2]][[2]][[2]][[1]][[1]]),
              b16_youngPyramidal = labels(dendro[[2]][[2]][[2]][[1]][[2]]),
              b17_pyramidal.glutamergic = labels(dendro[[2]][[2]][[2]][[2]][[1]]),
              b17_1_strayYoungPyramidal = labels(dendro[[2]][[2]][[2]][[2]][[1]][[1]][[1]]),
              b17_2_strayGlutamergic = labels(dendro[[2]][[2]][[2]][[2]][[1]][[1]][[2]][[2]][[1]]),
              b17_3_strayOldPyramidal = labels(dendro[[2]][[2]][[2]][[2]][[1]][[1]][[2]][[2]][[2]]),
              b17_4_youngPyramidal = labels(dendro[[2]][[2]][[2]][[2]][[1]][[2]][[1]]),
              b17_5_glutamergic = labels(dendro[[2]][[2]][[2]][[2]][[1]][[2]][[2]][[1]]),
              b17_6_oldPyramidal = labels(dendro[[2]][[2]][[2]][[2]][[1]][[2]][[2]][[2]]),
              b18_gabaergic = labels(dendro[[2]][[2]][[2]][[2]][[2]]),
              b18_1_youngGabaergic = labels(dendro[[2]][[2]][[2]][[2]][[2]][[1]]),
              b18_2_oldGabaergic = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]]),
              b18_3_oldPV = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[1]]), 
              b18_4_VIP. =  labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]),
              b18_5_mixedReln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]),
              b18_6_youngPV  = labels(dendro[[2]][[2]][[2]][[2]][[2]][[1]][[1]][[2]]),
              b18_7_youngSSTReln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[1]][[2]][[2]]),
              b18_8_calb_sst_reln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]),
              b18_9_reln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]),
              b18_10_sst_reln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]])
              )

comparisons =
"all
rest"

intersections = ''
foldTreshold = 10
#analysis
############################################
#convert original index to index in 
if (mode == 'index'){
  realGroups = vector(mode='list',length = length(groups))
  names(realGroups) = names(groups)
  for (i in 1:length(groups)){
    realGroups[[i]] = which(design$originalIndex %in% groups[[i]])
  }
}


if (mode == 'name'){
  realGroups = vector(mode='list',length = length(groups))
  names(realGroups) = names(groups)
  for (i in 1:length(groups)){
    realGroups[[i]] = which(design$sampleName %in% groups[[i]])    
  }
}



conn = textConnection(comparisons)
comp = read.csv(conn, header = F, strip.white = T)
  
groupNames = names(groups)

#uses medians
groupAverages = list()
for (i in realGroups){
  groupAverage = tryCatch({apply(exprData[i,], 2, median)},
                          error = function(cond){
                            print(i)
                            return(exprData[i,])
                          })
  
  groupAverages = c(groupAverages, list(groupAverage))
}
names(groupAverages)= groupNames

combinations = combn(1:length(groups),2)


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
  fChange = foldChange(groupAverages[[compare[1]]], groupAverages[[compare[2]]], foldTreshold)
  fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fChange$index], geneFoldChange= fChange$foldChange )
  fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]
  print(fileName)
  print(i)
  write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
}

if ('rest' %in% comp$V1){
  for (i in 1:length(realGroups)){
    fileName = paste('Rest/',names(realGroups)[i],sep ='')
    i = realGroups[[i]]
    groupAverage = tryCatch({apply(exprData[i,], 2, median)},
                            error = function(cond){
                              print(i)
                              return(exprData[i,])
                            })
    
    restAverage = tryCatch({apply(exprData[-i,], 2, median)},
                           error = function(cond){
                             print(paste('motherFuck!',i))
                             return(exprData[i,])
                           })
    fChange = foldChange(groupAverage, restAverage, foldTreshold)
    fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fChange$index], geneFoldChange= fChange$foldChange )
    fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]
    print(fileName)
    print(i)
    write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
    
  }
}

bbbb[c(length(bbbb)/2,length(bbbb)/2+1)]
aaaa[c(length(aaaa)/2,length(aaaa)/2+1)]


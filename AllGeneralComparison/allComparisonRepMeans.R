loc = getwd()
setwd('..')
parent = getwd()
setwd(loc)


source(paste0(parent,'/ogbox.r'))
require(reshape)
require(cluster)

foldChange2 = function (group1, group2, f = 10){

  
  groupAverage1 = group1
  
  
  
  groupAverage2 = tryCatch({apply(group2, 2, median)},
                           error = function(cond){
                             print('fuu')
                             return(group2)
                           })
  
  g19 = groupAverage1 < 9.5 & groupAverage1 > 8
  g16 = groupAverage1  < 6
  g29 = groupAverage2 < 9.5 & groupAverage2 > 8
  g26 = groupAverage2 < 6
  
  
  
  tempGroupAv2 = vector(length = length(groupAverage2))
  
  tempGroupAv2[g26 & g19] = tryCatch({(apply(group2[, g26 & g19], 2, max))},
                                     error = function(cond){
                                       print('I hate you damn it!')
                                       if (is.null(nrow(group2))){
                                           return(group2[g26 & g19])
                                       }else{ return(max(group2[, g26 & g19]))
                                     }})

  tempGroupAv2[g16 & g29] = tryCatch({(apply(group2[, g16 & g29], 2, mim))},
                                     error = function(cond){
                                         print('I hate you damn it!')
                                         if (is.null(nrow(group2))){
                                             return(group2[g16 & g29])
                                         }else{ return(min(group2[, g26 & g19]))
                                         }})
  

  #groupAverage1[5124]
  #groupAverage2[5124]
  
  
  #groupAverage1[7067]
  #groupAverage2[7067]
  
  add1 = g19 & g26 & groupAverage1>tempGroupAv2
  add2 = g29 & g16 & tempGroupAv2>groupAverage1 
  
  
  fold = (groupAverage1 - groupAverage2)
  chosen =  which({(fold >= (log(f)/log(2))) & !(g19 & g26) } | {(fold <= log(1/f)/log(2)) &  !(g29 & g16)}| add1 | add2)
  return(
    data.frame(index = chosen, foldChange = fold[chosen])
  )
}

giveSilhouette = function(daGeneIndex, groupInfo1, groupInfo2){
    clustering = as.integer(rep(1,nrow(design))*(1:nrow(design) %in% groupInfo1)+1)
    clustering = clustering[1:nrow(design) %in% c(groupInfo1, groupInfo2)]
    data = (exprData[ (1:nrow(design) %in% c(groupInfo1, groupInfo2)),  daGeneIndex])
    cluster = list(clustering = clustering, data = data)    
    silo = silhouette(cluster,dist(data))
    return(mean(silo[,3]))
}


#Data preperation
####################

design = read.table(paste0(parent,'/Data/normalizedDesign.csv'),header=T,sep='\t')
#design = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/Data/normalizedDesignEBENISIKIYIM.csv')
#write.table(design,quote=F,file='normalizedDesign2.csv', sep = "\t", row.names = F)

design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)


allDataPre = read.csv(paste0(parent,'/Data/mostVariableQuantileNormalized'), header = T)
geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData = exprData
design = design[match(colnames(exprData),make.names(design$sampleName),),]
rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<5)

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]

#disclusion is done below when comparing it to rest.
#discludeSamples = which(design$age<=14)
#exprData = exprData[, -discludeSamples]
#design = design[-discludeSamples, ]
exprData = t(exprData)
bok = cbind(design,exprData)[,1:18]

# get replicate means
# a terrible way to preallocate
newExpr = exprData[1:length(unique(design$originalIndex)),]
indexes = unique(design$originalIndex)
for (i in 1:length(indexes)){
  newExpr[i, ] = apply(exprData[design$originalIndex == indexes[i],], 2,mean)
}

newDesign = design[match(indexes,design$originalIndex),]


# deal with groups
nameGroups = list(someNaming = newDesign$someNaming, 
                  ourNamingGeneral = newDesign$ourNamingGeneral,
                  ourNaming = newDesign$ourNaming,
                  cellKind = newDesign$cellKind,
                  justGaba = newDesign$justGaba,
                  justGabaPV = newDesign$justGabaPV,
                  Cortex = newDesign$Cortex,
                  Amygdala = newDesign$Amygdala,
                  Striatum = newDesign$Striatum,
                  Midbrain = newDesign$Midbrain,
                  CortexWGaba = newDesign$CortexWGaba,
                  CortexWInter = newDesign$CortexWInter,
                  Cortex.SingleGroupInterGaba = newDesign$Cortex.SingleGroupInterGaba,
                  CortexWGabaDetailed = newDesign$CortexWGabaDetailed)

stepi = 1
for (i in nameGroups){
    groupNames = trimNAs(unique(i))
    realGroups = vector(mode = 'list', length = length(groupNames))
    names(realGroups) = groupNames
    for (j in 1:length(groupNames)){
        realGroups[[j]] = which(i == groupNames[j])           
    }
    groupAverages = list()
    
    
    for (j in realGroups){
        groupAverage = tryCatch({apply(newExpr[j,], 2, mean)},
                                error = function(cond){
                                    print('fuu')
                                    return(newExpr[j,])
                                })
        groupAverages = c(groupAverages, list(groupAverage))
    }
    
    names(groupAverages)= groupNames
    groupAverages = t(as.data.frame(groupAverages))
    
    dir.create('Rest/InBetween/')
    dir.create("Rest/InBetween/Marker")
    for (j in 1:nrow(groupAverages)){
        dir.create('Rest/InBetween/' + names(nameGroups)[stepi] + '/')
        fileName = 'Rest/InBetween/' + names(nameGroups)[stepi] + '/' +  names(realGroups)[j]
        dir.create('Rest/InBetween/Marker/' + names(nameGroups)[stepi] + '/')
        fileName2 = 'Rest/InBetween/Marker/' + names(nameGroups)[stepi] + '/'+ names(realGroups)[j]
        
        #find markers
        isMarker = vector(length = ncol(groupAverages))
        for (t in 1:ncol(groupAverages)){
            isMarker[t] = all(groupAverages[-j, t] + log(10, base=2) < groupAverages[j,t])
        }
        fMarker = data.frame(geneData$Gene.Symbol[isMarker], groupAverages[j,isMarker], tryCatch({apply(groupAverages[-j,isMarker],2,max)}, error = function(e){groupAverages[-j,isMarker]}))
        fChange = foldChange2(groupAverages[j, ], groupAverages[-j,] )
        fChangePrint = data.frame(geneNames = geneData$Gene.Symbol[fChange$index], geneFoldChange= fChange$foldChange )
        fChangePrint = fChangePrint[order(fChangePrint$geneFoldChange, decreasing=T) ,]
        
        #silhouette
        groupInfo1 = which(design[,names(nameGroups)[stepi]] == names(realGroups)[j])
        groupInfo2 = which(design[,names(nameGroups)[stepi]] != names(realGroups)[j] & !is.na(design[,names(nameGroups)[stepi]]))
        
        silo = vector(length = nrow(fChangePrint))
        for (t in 1:nrow(fChangePrint)){
           silo[t] = giveSilhouette(which(geneData$Gene.Symbol == fChangePrint$geneNames[t]),
                                    groupInfo1,
                                    groupInfo2)
        }
        
        fChangePrint = cbind(fChangePrint, silo)
        
        print(fileName)
        print(i)
        write.table(fChangePrint, quote = F, row.names = F, col.names = F, fileName)
        write.table(fMarker, quote = F, row.names = F, col.names = F, fileName2)
        
    }
    stepi = stepi + 1
}
    
    

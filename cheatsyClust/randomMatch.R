require(fields)
require(cluster)

trimElement = function (aVector,e){
  return(aVector[aVector!=e])
}

giveGI = function(daGene){
  return(match(daGene, geneData$Gene.Symbol))  
}

#filenames = list.files("Results",include.dirs = FALSE)
filenames2 = list.files("Results/min8FinerGroups",include.dirs = FALSE)

#ldf <- lapply(paste('Results/', filenames, sep = ''), read.table)
ldf2 <-lapply(paste('Results/min8FinerGroups/', filenames2, sep = ''), read.table)
                  
  
                  


geneList2 = vector(mode = 'character', length = 68000)
step = 1
for (i in ldf2){
  length(i$V1)
  geneList2[step:(step + length(i$V1)-1)] = as.character(i$V1)
  step = step+length(i$V1)
}

geneList2 = trimElement(geneList2,'')
length(unique(geneList2))
geneList2 = unique(geneList2)
geneList = geneList2

step = 1
#for (i in ldf1){
#  length(i$V1)
#  geneList[step:(step + length(i$V1)-1)] = as.character(i$V1)
#  step = step+length(i$V1)
#}

#geneList = trimElement(geneList,'')
#geneList = unique(geneList)
#length(geneList)

#Matched Random
################################################
design = read.table('C:/Users/Ogan/Dropbox/Rotation 3/Data/normalizedDesign.csv',header=T,sep='\t')


design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)

allDataPre = read.csv(dpaste(parent,'Data/mostVariableQuantileNormalized'), header = T)
geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData = exprData
design = design[match(colnames(exprData),design$sampleName,),]
rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<5)


exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]
selectGenes = giveGI(geneList)
selectGeneData = geneData[selectGenes, ]
selectExprData = exprData[selectGenes, ]

#Defining groups again
##########################
groups = list(#pyramidalN = c(40, 43, 44,45, 46, 47, 48, 49), 
              pyramidalCorticospineN = c(43, 44),
              pyramidalCingulateN = c(45, 47),
              pyramidalSomatosensoryN = c(46),
              pyramidalAmygdalaN = c(48),
              pyramidalHippo = 49,
              cholinergicN = c(13, 11), 
              spinyN = c(15, 16), 
              #gabaN = c(58, 62, 63, 64, 59, 56, 55, 54),
              gabaPV = c(58, 62, 63, 64),
              gabaReln = 59,
              gabaSSTReln = 57,
              gabaSSTRelnCalb = 56,
              gabaVIPReln = c(54, 55),
              glutaN = 53,
              golgiN = 17,
              granuleN = 20,
              interN = 14,
              #motorCholinN = c(10,12),
              motorCholinStemN = 10,
              motorCholinSpinalN = 12,
              mixedN = c(26, 34),
              basketN = 19,
              ubcN = 18,
              purkinjeG = 25,
              
              
              bergmanG = 27,  
              #oligoG = c(21, 22, 23, 24, 35), 
              matureOligoG = c(21, 22, 35),
              mixedOligoG = c(23, 24),
              astrocytesG = c(31, 32)
)

  realGroups = vector(mode='list',length = length(groups))
  names(realGroups) = names(groups)
  for (i in 1:length(groups)){
    realGroups[[i]] = which(design$originalIndex %in% groups[[i]])
  }

##########################
#silhouette with the selected genes. realGroups is different than the one in cheatsyClust.r
clustering = vector(mode = 'integer', length = max(unlist(realGroups)))
for (i in 1:length(realGroups)){
  clustering[realGroups[[i]]] = i
}
data = t(exprData[selectGenes, clustering > 0])
cluster = list(clustering = clustering[clustering > 0], data = data)
silo = silhouette(cluster, dist(data))
chosenSilo = mean(silo[,3])




selectMean = apply(selectExprData, 1, var)
selectVar = apply(selectExprData, 1, mean)
selectStat = apply(selectExprData, 1, var) / apply(selectExprData, 1, mean)
allStat = apply(exprData, 1, var) / apply(exprData, 1, mean)

binSep = seq(0, 2, 0.05)
bin = vector(mode = 'list', length = 40)

a=hist(selectStat,seq(0,2,0.05))
a$counts

count = 0
randomSilos = vector(mode = 'double', length = 1000)
for (n in 1:1000){
  randomGenes = exprData[sample(1:nrow(exprData),nrow(selectExprData)),clustering > 0]
  data = t(randomGenes)
  cluster = list(clustering = clustering[clustering > 0], data = data)
  silo = silhouette(cluster, dist(data))
  randomSilos[n] = mean(silo[,3])
  if (mean(silo[,3])<chosenSilo){
    count = count+1
  }
}




# count = 0
# randomSilos = vector(mode = 'double', length = 1000)
# for (n in 1:1000){
#   for (i in 1:length(a$counts)){
#     qualify = which(allStat >= binSep[i] & allStat <= binSep[i+1])
#     bin[[i]] = as.integer(sample(qualify, a$counts[i]))
#   }
#   
#   randomExpr = exprData[unlist(bin), ]
#   randomGeneData = geneData [unlist(bin), ]
#   
#   #random silhouette 
#   data = t(exprData[unlist(bin), clustering > 0])
#   cluster = list(clustering = clustering[clustering > 0], data = data)
#   silo = silhouette(cluster, dist(data))
#   randomSilos[n] = mean(silo[,3])
#   if (mean(silo[,3])<chosenSilo){
#     count = count+1
#   }
#   
#   
# }

hist(trimElement(randomSilos,0))

require(ggplot2)

( chosenSilo - mean(randomSilos))/sd(randomSilos)





#Heatmap
#######################################################
dpaste = function (...){
  paste(..., sep='')
}
getParent = function(step = 1){
  wd = getwd()
  setwd('..')
  parent = getwd()
  setwd(wd)
  return(parent)
}
parent = getParent()
parent = dpaste(parent, '/')

allDataPre = read.csv(dpaste(parent,'Data/mostVariableQuantileNormalized' ), header = T)

set.seed(3)
useColors = c('darkgreen',
              'black',
              'blue',
              'orange',
              'brown',
              'red',
              'green',
              'yellow',
              'darkgray',
              'cyan',
              'darkseagreen',
              'deeppink',
              'aquamarine4',
              'darkkhaki',
              'darkolivegreen',
              'red4',
              'tan2',
              'thistle',
              'slateblue4',
              'powderblue',
              'lightslategray',
              'maroon',
              'thistle4'
)



#design = read.csv('normalizedDesignNonOrdered')
design = read.table('normalizedDesign.csv',header=T,sep='\t')


design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)


geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData = 2^exprData
design = design[match(colnames(exprData),design$sampleName,),]

rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<2^5)

exprData = exprData[-discludeGenes,] #order here is important as these were discluded when calculating silhouettes as well
geneData = geneData[-discludeGenes,]
selectGenes = giveGI(geneList)

exprData = exprData[selectGenes, ]
geneData = geneData[selectGenes, ]

discludeSamples = which(
  design$age < 14 | 
    design$cellType2 == 'Astroglia' |
    design$originalIndex %in% c(6,7,8)
)

exprData = exprData[, -discludeSamples]
design = design [-discludeSamples, ]

corr = cor(exprData)
######################################
refColor = rep('white', ncol(exprData))

reference = as.character(design$reference)



refNames = unique(reference)

for (i in 1:length(refNames)){
  refColor[grepl(refNames[i], reference)] = useColors[i]
}
legend("bottom", legend= refNames,
       fill = c(useColors))

#############################
cellTypeColor = rep('white', ncol(exprData))
cellType = as.character(design$cellType2)
cellType[grepl('^Neurons', cellType)] = 'Mixed Neurons'
cellType[grepl('Medium Spiny Neurons', cellType)] = 'Medium Spiny Neurons'
cellType[grepl('^Cholinergic', cellType)] = 'Cholinergic'
cellType[grepl('Dopaminergic', cellType)] = 'Dopaminergic'
cellType[grepl('Motor', cellType)] = 'Motor+Cholinergic'
cellType[grepl('Pyramidal', cellType)] = 'Pyramidal'
cellType[grepl('Glutamatergic', cellType)] = 'Glutamatergic'
cellTypeNames = unique(cellType)

cellTypeNames = c('Oligodendrocytes','Mature Oligodendrocytes', 'Oligodendrocyte Precursors',"Mixed Oligodendrocytes", 'Astroglia', 'Astrocytes', 'Bergman Glia', 'Motor+Cholinergic', 'Cholinergic', 'Medium Spiny Neurons', 'Glutamatergic', 'Stellate Basket Cells','Golgi Cells', 'Pyramidal', 'Purkinje Cells', 'Interneurons', 'Granule Cells', 'Unipolar Brush cells', 'Mixed Neurons', 'GABAergic Interneurons, PV+', "GABAergic Interneurons, SST+ Reln+ Calb1+", "GABAergic Interneurons, SST+ Reln+", "GABAergic Interneurons, Reln+", "GABAergic Interneurons, VIP+ Reln+")
typeColors =    c('darkgreen'       ,'forestgreen'            , 'greenyellow'               ,'green'                 , 'yellow'   , 'tan'       , 'palegreen'   , 'red'              , 'darkorange' , 'blanchedalmond'      , 'slategray'    , 'mediumpurple4'        ,'orchid'     , 'turquoise', 'purple'        , 'pink'        , 'thistle'      , 'powderblue'          , 'black'        , colorRampPalette(c("firebrick4",'coral'))(n = 5))






for (i in 1:length(cellTypeNames)){
  cellTypeColor[cellType %in% cellTypeNames[i]] = typeColors[i]
}

legend("bottomleft", legend= cellTypeNames,
       fill = c(typeColors))

####################################################



#####################################################
ageUseColors = colorRampPalette(c("red",'white', "blue"))(n = 12)
ageColor = rep('white', ncol(exprData))
age = as.numeric(design$age)
ageNames = sort(unique(age))

for (i in 1:length(ageNames)){
  ageColor[grepl(ageNames[i], age)] = ageUseColors[i]
}

legend("bottomright", legend= ageNames,
       fill = c(ageUseColors), cex=0.8)
#############################################################


######################################
sexColors = rep('white', ncol(exprData))



####################################

allColors=cbind(reference=refColor, age= ageColor, cellType=cellTypeColor)


source('C:\\Users\\Ogan\\Dropbox\\Rotation 3\\heatmap.3.r')
#######################################################################
palette = colorRampPalette(c("navyblue", "blanchedalmond","firebrick4"))(n = 1000)
palette = colorRampPalette(c("#C7BC69","#59311C"))(n = 1000)

colnames(corr) = paste(design$originalIndex, rownames(corr),sep='.')
a = heatmap.3(corr, trace = "none", Rowv = F, Colv = T, 
              col = palette, ColSideColors = (as.matrix(allColors)), cexCol=1,margins = c(7,5))
legend("bottomleft", legend= cellTypeNames,
       fill = c(typeColors), cex=0.8)

samplesss = as.character(design$sampleName[which(grepl('GABA',design$cellType))])
sampless = which(colnames(exprData) %in% samplesss)
geness = which(geneData$Gene.Symbol %in% c('Reln','Pvalb','Sst','Calb1','Vip'))
genessNames = geneData$Gene.Symbol[geness]

table=exprData[geness,sampless]
rownames(table)=genessNames

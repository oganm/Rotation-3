#Startup
###############

allSilhouette = read.csv('Results/allSilhouette',header = T)
tsilly = allSilhouette[, 4:ncol(allSilhouette)]
plot(density(unlist(tsilly)))

hist(unlist(tsilly))

#selectGenes = (apply(tsilly>0.5, 1, any))
selectGenes = apply(tsilly,2,order,decreasing = T)

selectGenes = selectGenes[1:50,]
dim(selectGenes) = NULL
selectGenes = unique(selectGenes)


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


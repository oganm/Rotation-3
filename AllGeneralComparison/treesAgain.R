#Startup
###############
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
design = read.table(parent +'/Data/normalizedDesignSex.csv',header=T,sep='\t')


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

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]

#random = sample(1:nrow(exprData),10)
#exprData = exprData[random,]
#geneData = geneData[random,]
#exprData = exprData[discludeGenes,]
#geneData = geneData[discludeGenes,]

# discludeSamples = which(
#   design$age < 14 | 
#     design$cellType2 == 'Astroglia' |
#   design$originalIndex %in% c(6,7,8)
#   )

discludeSamples = which(is.na(design$ourNamingGeneral))
exprData = exprData[, -discludeSamples]
design = design [-discludeSamples, ]

filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/ourNamingGeneral",include.dirs = FALSE)
fileContents = lapply("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/ourNamingGeneral"+'/'+ filenames, read.table)
geneList = vector(mode = 'list', length = length(fileContents))
names(geneList) = filenames
for (i in 1:length(fileContents)){
    geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>0])
}

puristList = vector(mode = 'list', length = length(geneList))
for (i in 1:length(geneList)){
    puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
}

names(puristList) = names(geneList)

foldGenes = puristList

exprData=exprData[geneData$Gene.Symbol %in% unlist(puristList),]
geneData=geneData[geneData$Gene.Symbol %in% unlist(puristList),]

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
cellType = as.character(design$ourNamingGeneral)
#cellType[grepl('^Neurons', cellType)] = 'Mixed Neurons'
#cellType[grepl('Medium Spiny Neurons', cellType)] = 'Medium Spiny Neurons'
#cellType[grepl('^Cholinergic', cellType)] = 'Cholinergic'
#cellType[grepl('Dopaminergic', cellType)] = 'Dopaminergic'
#cellType[grepl('Motor', cellType)] = 'Motor+Cholinergic'
#cellType[grepl('Pyramidal', cellType)] = 'Pyramidal'
#cellType[grepl('Glutamatergic', cellType)] = 'Glutamatergic'
cellTypeNames = unique(cellType)

cellTypeNames = c('Oligodendrocytes','Mature Oligodendrocytes', 'Oligodendrocyte Precursors',"Mixed Oligodendrocytes", 'Astroglia', 'Astrocytes', 'Bergman Glia', 'Motor+Cholinergic', 'Cholinergic', 'Medium Spiny Neurons', 'Glutamatergic', 'Stellate Basket Cells','Golgi Cells', 'Pyramidal', 'Purkinje Cells', 'Interneurons', 'Granule Cells', 'Unipolar Brush cells', 'Mixed Neurons', 'Microglia','GABAergic Interneurons, PV+', "GABAergic Interneurons, SST+ Reln+ Calb1+", "GABAergic Interneurons, SST+ Reln+", "GABAergic Interneurons, Reln+", "GABAergic Interneurons, VIP+ Reln+")
typeColors =    c('darkgreen'       ,'forestgreen'            , 'greenyellow'               ,'green'                 , 'yellow'   , 'tan'       , 'palegreen'   , 'red'              , 'darkorange' , 'blanchedalmond'      , 'slategray'    , 'mediumpurple4'        ,'orchid'     , 'turquoise', 'purple'        , 'pink'        , 'thistle'      , 'powderblue'          , 'black'        , 'white',colorRampPalette(c("firebrick4",'coral'))(n = 5))
typeColors = typeColors[1:25]

cellTypeNames = c('Oligo',      'Astro'    , 'Bergman', 'MotorCholin', 'Cholinergic', 'Spiny',                'Gluta',    'Basket',          'Golgi', 'Pyramidal', 'Purkinje', 'Inter', 'Granule', 'Ubc',        'Microglia','Gaba')
typeColors =    c('darkgreen', 'yellow'   , 'palegreen'   , 'red'   , 'darkorange' , 'blanchedalmond'      , 'slategray', 'mediumpurple4'  ,'orchid', 'turquoise', 'purple',   'pink' , 'thistle', 'powderblue', 'white',"firebrick4")



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
       fill = c(ageUseColors), cex=0.5)
#############################################################


######################################
sexColors = rep('white', ncol(exprData))



####################################

allColors=cbind(reference=refColor, age= ageColor, cellType=cellTypeColor)

w
source('C:\\Users\\Ogan\\Dropbox\\Rotation 3\\heatmap.3.r')
#######################################################################
palette = colorRampPalette(c("orange", "white"))(n = 1000)
palette = colorRampPalette(c("#C7BC69","#59311C"))(n = 1000)

colnames(corr) = paste(design$originalIndex, rownames(corr),sep='.')
a = heatmap.3(corr, trace = "none", Rowv = T, Colv = T, 
              col = palette, ColSideColors = (as.matrix(allColors)), cexCol=1,margins = c(7,5),dendrogram = 'column')
legend("bottomleft", legend= cellTypeNames,
       fill = c(typeColors), cex=0.6)
legend("bottomright", legend= ageNames,
       fill = c(ageUseColors), cex=0.6)

 samplesss = as.character(selectDes$sampleName[which(grepl('GABA',design$cellType))])
sampless = which(colnames(exprData) %in% samplesss)
geness = which(geneData$Gene.Symbol %in% c('Reln','Pvalb','Sst','Calb1','Vip'))
genessNames = geneData$Gene.Symbol[geness]

table=exprData[geness,sampless]
rownames(table)=genessNames


a = heatmap.3(corr, trace = "none", Rowv = T, Colv = T, 
              col = palette, ColSideColors = (as.matrix(allColors)), cexCol=1,margins = c(7,5),dendrogram = 'column')
legend("bottomleft", legend= cellTypeNames,
       fill = c(typeColors), cex=0.6)


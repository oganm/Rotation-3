library(rgl)

source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')
parent = getParent()
allDataPre = read.csv(parent + 'Data/mostVariableQuantileNormalized', header = T)
require(ggplot2)

design = read.table(parent + 'Data/normalizedDesign.csv',header=T,sep='\t')

design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)




geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]

rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<5)

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]

design = design[match(colnames(exprData),design$sampleName,),]

filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/justGabaPV",include.dirs = FALSE)
fileContents = lapply("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/justGabaPV"+'/'+ filenames, read.table)
geneList = vector(mode = 'list', length = length(fileContents))
names(geneList) = filenames
for (i in 1:length(fileContents)){
  geneList[[i]] = as.character(fileContents[[i]]$V1)
}
geneList


#filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/ourNamingGeneral",include.dirs = FALSE)
#fileContents = lapply("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/ourNamingGeneral"+'/'+ filenames, read.table)
#geneList = vector(mode = 'list', length = length(fileContents))
#names(geneList) = filenames
#for (i in 1:length(fileContents)){
#  geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>0])
#}

puristList = vector(mode = 'list', length = length(geneList))
for (i in 1:length(geneList)){
  puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
}



names(puristList) = names(geneList)


exprData=exprData[geneData$Gene.Symbol %in% unlist(puristList[[2]]),]
geneData=geneData[geneData$Gene.Symbol %in% unlist(puristList[[2]]),]


discludeSamples = which(is.na(design$ourNamingGeneral))
exprData = exprData[, -discludeSamples]
design = design [-discludeSamples, ]



exprData = t(exprData)
#exprData=exprData[,geneData$Gene.Symbol %in% c('Mobp', 'Slc1a3', 'Slc1a2', 'Aqp4', 'Mag', 'Mog', 'Tac1', 'Neurod6') ]
#geneData = geneData[geneData$Gene.Symbol %in% c('Mobp', 'Slc1a3', 'Slc1a2', 'Aqp4', 'Mag', 'Mog', 'Tac1', 'Neurod6'),


#variance = apply(exprData, 1, var)

#variance = apply(expData, 2, var)

boxy=boxplot(t(exprData))
#variance = boxy$stats[4,]

#boxplot(t(exprData[order(variance),]))
#transition is smooth. adding it as a numerical variable rather than categorical

#pcaData = rbind(exprData, variance)
#pcaData = rbind(expData, variance)

#prcs = prcomp((pcaData), scale=T)
prcs = prcomp(exprData, scale = T)
a=prcs$rotation

cellTypeColor = rep('white', nrow(exprData))
cellType = as.character(design$someNaming)
cellTypeNames = c('Oligo',      'Astro'    , 'Bergman', 'MotorCholin', 'Cholinergic', 'Spiny',                'Gluta',    'Basket',          'Golgi', 'Pyramidal', 'Purkinje', 'Inter', 'Granule', 'Ubc',        'Microglia','Gaba', 'GabaPV')
typeColors =    c('darkgreen', 'yellow'   , 'palegreen'   , 'red'   , 'darkorange' , 'blanchedalmond'      , 'slategray', 'mediumpurple4'  ,'orchid', 'turquoise', 'purple',   'pink' , 'thistle', 'powderblue', 'black',"firebrick4", 'blue')
for (i in 1:length(cellTypeNames)){
  cellTypeColor[cellType %in% cellTypeNames[i]] = typeColors[i]
11}
frame = data.frame(prcs$x[,1], prcs$x[,2], prcs$x[,3],  ty = trimNAs(design$someNaming))

plot3d(frame, col = cellTypeColor, size = 10)                                                                                                                                                                                  # panel.grid.minor = theme_blank(),



#which(abs(a[21086,]) %in% max(abs(a[21086,])))

max(abs(a[,1]))
Xist=which(geneData$Gene.Symbol %in% 'Xist')[1]
sex=exprData[,Xist]>5.6
frame = data.frame(prcs$x[,2], prcs$x[,1], prcs$x[,3],  ty = design$someNaming)


cbPalette <- typeColors
(ggplot(frame, aes(y=prcs.x...2., x= prcs.x...1. , color = ty)))+geom_point(size=7)+ scale_colour_manual(values=cbPalette)#+scale_shape_manual(values=unlist(lapply(c('O','I'), utf8ToInt)))+theme_bw()+opts(axis.line = theme_segment(colour = "black"),
                                                                                                                                                                                  # panel.grid.major = theme_blank(),
                                                                                                                                                                                  # panel.border = theme_blank(),
                                                                                                                                                                                  # panel.background = theme_blank()))


#(ggplot(frame, aes(y=prcs.x...2., x= prcs.x...1. , color = ty))+geom_point(size=5))

#t(pcaData)

#t(as.matrix(pcaData[,1])) %*% as.matrix(prcs$rotation[,1])

#View(prcs$x)
#dim(t(as.matrix(exprData)) %*% as.matrix(prcs$rotation))
  
which(geneData$Gene.Symbol[order(abs(a[,2]),decreasing=T)] %in% 'Xist')
a[Xist,2]
max(abs(a[,2]))

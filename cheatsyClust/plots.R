require(scatterplot3d)
require(plot3D)
library(rgl)
#non expressed. reload first. 
#toPlot = exprData[discludeGenes,(design$age>=14)&(!design$cellType2=='Astroglia')&(!design$originalIndex %in% c(6,7,8))]
#random
#toPlot = exprData[sample(1:nrow(exprData),nrow(selectExprData)),(design$age>=14)&(!design$cellType2=='Astroglia')&(!design$originalIndex %in% c(6,7,8))]
#selected         
toPlot = selectExprData[,(design$age>=14)&(!design$cellType2=='Astroglia')&(!design$originalIndex %in% c(6,7,8))]
desSele = design[(design$age>=14)&(!design$cellType2=='Astroglia')&(!design$originalIndex %in% c(6,7,8)),]


toPlot = toPlot[selectGeneData$Gene.Symbol %in% c('Mobp', 'Slc1a3', 'Slc1a2', 'Aqp4', 'Mag', 'Mog', 'Tac1', 'Neurod6'),]



prcs = prcomp(t(toPlot), scale = T)
daFrame = data.frame(prcs$x[,3], prcs$x[,1], prcs$x[,2],prcs$x[,4], 
                   #index = design$originalIndex[c(realGroups[[1]], realGroups[[9]])], 
                   #anatomy = design$anatomical[c(realGroups[[1]], realGroups[[9]])],
                   cellType = desSele$cellType2,
                   age = desSele$age
)
scatterplot3d(x=  prcs$x[,1], y =  prcs$x[,2], z=  prcs$x[,3], color= cellTypeColor, pch= 16)
cellType = as.character(desSele$cellType2)
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

legend("topleft", legend= cellTypeNames,
       fill = c(typeColors), cex=0.6)

plot3d(x=  prcs$x[,1], y =  prcs$x[,2], z=  prcs$x[,3],size= 8, col =cellTypeColor )


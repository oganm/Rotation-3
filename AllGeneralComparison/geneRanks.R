foldChangeOut = function (group1, group2){
  return(fold = uselog^(group1 - group2))
}


#Data preperation
####################
allDataPre = read.csv('allNormalized', header = T)

design = read.table('normalizedDesign.csv',header=T,sep='\t')


design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)
design = design[match(colnames(allDataPre),design$sampleName,),]

geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData[exprData<6] = 6
exprData = 2^exprData
uselog=5
exprData = log(exprData, base = uselog)
exprData = round(exprData)

design = design[match(colnames(exprData),design$sampleName,),]


#discludeSamples = which(design$age<=14)
#exprData = exprData[, -discludeSamples]
#design = design[-discludeSamples, ]
exprData = t(exprData)
#########################################


#define groups
#####################

dendro = a$colDendrogram

#groupsCluster
####################
#groups = list(b1_oligodendro = labels(dendro[[1]][[1]][[1]]), 
#              b2_astrocytes = labels(dendro[[1]][[1]][[2]]),
#              b2_1_oldAstro = labels(dendro[[1]][[1]][[2]][[1]]),
#              b2_2_youngAstro = labels(dendro[[1]][[1]][[2]][[1]]),
#              b3_oligo = labels(dendro[[1]][[2]][[1]]),
#              b3_1_oligoPre = labels(dendro[[1]][[2]][[1]][[1]]),
#              b3_2_oligo =labels(dendro[[1]][[2]][[1]][[2]]),
#              b4_astro.bergman = labels(dendro[[1]][[2]][[2]][[1]]),
#              b4_1_astro = labels(dendro[[1]][[2]][[2]][[1]][[1]]),
#              b4_2_astro.bergman = labels(dendro[[1]][[2]][[2]][[1]][[2]]),
#              b5_oldOligo = labels(dendro[[1]][[2]][[2]][[2]]),
#              b5_1_LoneMatureOligo = labels(dendro[[1]][[2]][[2]][[2]][[1]][[1]]),
#               b5_2_mixedOligo = labels(dendro[[1]][[2]][[2]][[2]][[1]][[2]]),
#               b5_3_mixedOligoCortex = labels(dendro[[1]][[2]][[2]][[2]][[2]][[1]]),
#               b6_purkinje = labels(dendro[[2]][[1]][[1]]),
#               b7_ubc = labels(dendro[[2]][[1]][[2]][[1]]),
#               b8_granule = labels(dendro[[2]][[1]][[2]][[2]][[1]][[1]]),
#               b9_basket = labels(dendro[[2]][[1]][[2]][[2]][[1]][[2]][[1]]),
#               b10_golgi = labels(dendro[[2]][[1]][[2]][[2]][[1]][[2]][[2]]),
#               b11_cholinergic = labels(dendro[[2]][[1]][[2]][[2]][[2]]),
#               b11_1_forebrainCholi = labels(dendro[[2]][[1]][[2]][[2]][[2]][[1]]),
#               b11_2_striatum.forebrainCholi = labels(dendro[[2]][[1]][[2]][[2]][[2]][[2]][[1]]),
#               b11_3_choli.motor = labels(dendro[[2]][[1]][[2]][[2]][[2]][[2]][[2]]),
#               b12_spiny = labels(dendro[[2]][[2]][[1]][[1]][[1]]),
#               b13_mixed = labels(dendro[[2]][[2]][[1]][[1]][[2]]),
#               b14_pyramidal.inter = labels(dendro[[2]][[2]][[1]][[2]]),
#               b15_mixed = labels(dendro[[2]][[2]][[2]][[1]][[1]]),
#               b16_youngPyramidal = labels(dendro[[2]][[2]][[2]][[1]][[2]]),
#               b17_pyramidal.glutamergic = labels(dendro[[2]][[2]][[2]][[2]][[1]]),
#               b17_1_strayYoungPyramidal = labels(dendro[[2]][[2]][[2]][[2]][[1]][[1]][[1]]),
#               b17_2_strayGlutamergic = labels(dendro[[2]][[2]][[2]][[2]][[1]][[1]][[2]][[2]][[1]]),
#               b17_3_strayOldPyramidal = labels(dendro[[2]][[2]][[2]][[2]][[1]][[1]][[2]][[2]][[2]]),
#               b17_4_youngPyramidal = labels(dendro[[2]][[2]][[2]][[2]][[1]][[2]][[1]]),
#               b17_5_glutamergic = labels(dendro[[2]][[2]][[2]][[2]][[1]][[2]][[2]][[1]]),
#               b17_6_oldPyramidal = labels(dendro[[2]][[2]][[2]][[2]][[1]][[2]][[2]][[2]]),
#               b18_gabaergic = labels(dendro[[2]][[2]][[2]][[2]][[2]]),
#               b18_1_youngGabaergic = labels(dendro[[2]][[2]][[2]][[2]][[2]][[1]]),
#               b18_2_oldGabaergic = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]]),
#               b18_3_oldPV = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[1]]), 
#               b18_4_VIP.Reln =  labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]),
#               b18_5_mixedReln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]),
#               b18_6_youngPV  = labels(dendro[[2]][[2]][[2]][[2]][[2]][[1]][[1]][[2]]),
#               b18_7_youngSSTReln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[1]][[2]][[2]]),
#               b18_8_calb_sst_reln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]),
#               b18_9_reln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]),
#               b18_10_sst_reln = labels(dendro[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]])
# )
# 
# groupsCluster = vector(mode='list',length = length(groups))
# names(groupsCluster) = names(groups)
# for (i in 1:length(groups)){
#   groupsCluster[[i]] = which(design$sampleName %in% groups[[i]])    
# }

###########


#groupsOfficial
######################
groups = list(pyramidalN = c(40, 43, 44, 45, 46, 47, 48, 49),     
              cholinergicN = c(13, 11), 
              spinyN = c(15, 16), 
             # gabaN = c(58, 62, 63, 64, 59, 56, 55, 54),
              gabaPV = c(58, 62, 63, 64),
              gabaReln = 59,
              gabaSSTReln = 57,
              gabaSSTRelnCalb = 56,
              gabaVIPReln = c(54, 55),
              glutaN = 53,
              golgiN = 17,
              granuleN = 20,
              interN = 14,
              motorCholinN = c(10,12),
           #   mixedN = c(26, 34),
              basketN = 19,
              ubcN = 18,
             # neurons = c(40, 43, 44, 45, 46, 47, 48, 49, 13, 11, 15, 16, 58, 62, 63, 64, 59, 56, 55, 54, 58, 62, 63, 64, 59, 57, 56, 55, 54, 53, 17, 20, 14, 10, 12, 26, 34, 19, 18),
              
              
            #  astroG = c(28, 29),
              bergmanG = 27,  
              oligoG = c(21, 22, 23, 24, 35), 
              astrocytesG = c(31, 32), 
              purkinjeG = 25
             # glials = c(28,29,27,21,22,35,31,32 ,25, 23, 24)
)
groupsOfficial = vector(mode='list',length = length(groups))
names(groupsOfficial) = names(groups)
for (i in 1:length(groups)){
  groupsOfficial[[i]] = which(design$originalIndex %in% groups[[i]])
}

###############



#group medians
##############
groupAvO = vector(mode = 'double', length = length(groupsOfficial) * ncol(exprData))
dim(groupAvO) = c(length(groupsOfficial), ncol(exprData))

for (i in 1:length(groupsOfficial)){
  groupAverage = tryCatch({apply(exprData[groupsOfficial[[i]],], 2, median)},
                          error = function(cond){
                            print(i)
                            return(exprData[groupsOfficial[[i]],])
                          })
  
  groupAvO[i,] = groupAverage
}

# groupAvC = vector(mode = 'double', length = length(groupsCluster) * ncol(exprData))
# dim(groupAvC) = c(length(groupsCluster), ncol(exprData))

for (i in 1:length(groupsCluster)){
  groupAverage = tryCatch({apply(exprData[groupsCluster[[i]],], 2, median)},
                          error = function(cond){
                            print(i)
                            return(exprData[groupsCluster[[i]],])
                          })
  
  groupAvC[i,] =groupAverage
}

# ranksC = vector(mode = 'double', length = nrow(groupAvC) * ncol(groupAvC))
ranksO = vector(mode = 'double', length = nrow(groupAvO) * ncol(groupAvO))




#####################
#groupAvC
#groupAvO


#get fold changes
####################
fChangeO = vector(mode = 'double', length = length(groupsOfficial) * ncol(exprData))
# fChangeC = vector(mode = 'double', length = length(groupsCluster)  * ncol(exprData))
# dim(fChangeC) = c(length(groupsCluster), ncol(exprData))
dim(fChangeO) = c(length(groupsOfficial), ncol(exprData))

# n=1
# for (i in groupsCluster){
#   
#   groupAverage = tryCatch({apply(exprData[i,], 2, median)},
#                           error = function(cond){
#                             print(i)
#                             return(exprData[i,])
#                           })
#   
#   restAverage = tryCatch({apply(exprData[-i,], 2, median)},
#                          error = function(cond){
#                            print(paste('motherFuck!',i))
#                            return(exprData[i,])
#                          })
#   
#   fChange = foldChangeOut(groupAverage, restAverage)
#   fChangeC[n,] = fChange
#   n = n + 1  
# }

n=1
for (i in groupsOfficial){
  
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
  
  fChange = foldChangeOut(groupAverage, restAverage)
  fChangeO[n,] = fChange
  n = n + 1  
}
#################################






#get ranked 1 gene list
################
# ranksC = apply(-groupAvC,2,rank)
ranksO = apply(-groupAvO,2,rank)

# ranksC = ranksC == 1
ranksO = ranksO == 1



# EpicGenesC = list()
# EpicFoldC = list()
EpicGenesO = list()
EpicFoldO = list()
#also filters for a min fold change =3
# for (i in 1:nrow(ranksC)){
#   EpicGenesC[[names(groupsCluster)[i]]] = geneData$Gene.Symbol[ranksC[i,]&fChangeC[i,]>3]
#   EpicFoldC[[names(groupsCluster)[i]]] = fChangeC[i,ranksC[i,]&fChangeC[i,]>3]
# }

for (i in 1:nrow(ranksO)){
  EpicGenesO[[names(groupsOfficial)[i]]] = geneData$Gene.Symbol[ranksO[i,]&fChangeO[i,]>3]
  EpicFoldO[[names(groupsOfficial)[i]]] = fChangeO[i,ranksO[i,]&fChangeO[i,]>3]
  
}

######################

#for (i in 1:length(EpicGenesC)){
#  filename = paste('RankedResults/cluster/', names(EpicGenesC)[i], sep = '')
#  frame = data.frame(genes = EpicGenesC[[i]], fold = EpicFoldC[[i]])
#  frame = frame[order(frame$fold, decreasing = T),]
#  frame = frame[match(unique(frame$genes), frame$genes),]
#  print(i)
#  write.table(frame, quote = F, row.names = F, col.names = F, file = filename)
#}

for (i in 1:length(EpicGenesO)){
  filename = paste('RankedResults/original/', names(EpicGenesO)[i], sep = '')
  frame = data.frame(genes = EpicGenesO[[i]], fold = EpicFoldO[[i]])
  frame = frame[order(frame$fold, decreasing = T),]
  frame = frame[match(unique(frame$genes), frame$genes),]
  print(i)
  write.table(frame, quote = F, row.names = F, col.names = F, file = filename)
}
m=0
play(w)

heatNames = vector()
for (i in 1:length(EpicGenesO)){
  heatNames= c(heatNames, as.character(EpicGenesO[[i]]))
  
}
heatNames=unique(heatNames)

exprData2 = t(allDataPre[,4:ncol(allDataPre)])

heatData = exprData2[,apply(ranksO,2,any)]
heatData = t(exprData2[,geneData$Gene.Symbol %in% heatNames])
typeColors =    c('darkgreen'       ,'forestgreen'            , 'greenyellow'               ,'green'                 , 'yellow'   , 'tan'       , 'palegreen'   , 'red'              , 'darkorange' , 'blanchedalmond'      , 'slategray'    , 'mediumpurple4'        ,'orchid'     , 'turquoise', 'purple'        , 'pink'        , 'thistle'      , 'powderblue'          , 'black'        , colorRampPalette(c("firebrick4",'coral'))(n = 6))

groupColors = rep('white',165)
for (i in 1:length(groupsOfficial)){
  if (names(groupsOfficial)[i] %in% c('glials', 'neurons', 'gabaN')){
    next
  }
  groupColors[groupsOfficial[[i]]]=useColors[i]
  
}
  
matr = (as.matrix(groupColors))

heatmap.3(heatData,dendrogram='col',ColSideColors = matr,cexCol=1)

names(groupsOfficial)[!(names(groupsOfficial) %in% c('glials', 'neurons', 'gabaN'))]

bokbokbok = names(groupsOfficial)[!(names(groupsOfficial) %in% c('glials', 'neurons', 'gabaN'))]
okok = useColors[!(names(groupsOfficial) %in% c('glials', 'neurons', 'gabaN'))]
legend("bottomleft", legend = bokbokbok,
              fill = okok, cex=0.8)

i=4
require(gplots)

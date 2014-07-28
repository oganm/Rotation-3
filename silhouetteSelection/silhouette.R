source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')

require(tuneR)
require(cluster)
require(ggplot2)
t = seq(0, 3, 1/8000) #times in seconds if sample for 3 seconds at 8000Hz
u = (2^15-1)*sin(2*pi*440*t) #440 Hz sine wave that lasts t length seconds (here, 3 seconds) 
w = Wave(u, samp.rate = 8000, bit=16) #make the wave variable 


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

plotGene = function(daGene, groupInfo){
  group = as.integer(rep(1,nrow(design))*(1:nrow(design) %in% groupInfo)+1)
  daIndex = which(geneData$Gene.Symbol %in% daGene)
  fr = data.frame(expr = as.double(exprData[daIndex,]), group = group)
  (ggplot(fr, aes(x=expr,y=group))
   +geom_point(alpha=0.3))
  
}

plotGenes = function (daGenes, groupInfo){
  group = as.integer(rep(1,nrow(design))*(1:nrow(design) %in% groupInfo)+1)
  daIndex = match(daGenes, geneData$Gene.Symbol)
  fr = data.frame(as.double(exprData[daIndex,][1,]),as.double(exprData[daIndex,][2,]), group = factor(group))
  colnames(fr) = c('x','y','group')
  (ggplot(fr, aes(x=x,y=y))
   +geom_point(aes(color=group))
   +xlab(daGenes[1])
   +ylab(daGenes[2])
   )
}

giveGI = function(daGene){
  return(match(daGene, geneData$Gene.Symbol))  
}

giveSilhouette = function(daGeneIndex, groupInfo){
  clustering = as.integer(rep(1,nrow(design))*(1:nrow(design) %in% groupInfo)+1)
  data = t(exprData[daGeneIndex,])
  cluster = list(clustering = clustering, data = data)    
  silo = silhouette(cluster,dist(data))
  return(mean(silo[,3]))
}

parent = getParent()

#allDataPre = read.csv(dpaste(parent,'Data/mostVariableExp' ), header = T)
allDataPre = read.csv(dpaste(parent,'Data/mostVariableQuantileNormalized' ), header = T)

design = read.table('C:/Users/Ogan/Dropbox/Rotation 3/Data/normalizedDesign.csv',header=T,sep='\t')

design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)


geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]

design = design[match(colnames(exprData),design$sampleName,),]
rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<5)

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]


discludeSamples = discludeSamples = which(is.na(design$ourNaming))


exprData = exprData[, -discludeSamples]
design = design [-discludeSamples, ]

# groups = list(pyramidalN = c(40, 43, 44,45, 46, 47, 48, 49),     
#               cholinergicN = c(13, 11), 
#               spinyN = c(15, 16), 
#               gabaN = c(58, 62, 63, 64, 59, 56, 55, 54),
#               gabaPV = c(58, 62, 63, 64),
#               gabaReln = 59,
#               gabaSSTReln = 57,
#               gabaSSTRelnCalb = 56,
#               gabaVIPReln = c(54, 55),
#               glutaN = 53,
#               golgiN = 17,
#               granuleN = 20,
#               interN = 14,
#               motorCholinN = c(10,12),
#               mixedN = c(26, 34),
#               basketN = 19,
#               ubcN = 18,
#               purkinjeN = 25,
#               
#               neurons = c(40, 43, 44, 45, 46, 47, 48, 49, 13, 11, 15, 16, 58, 62, 63, 64, 59, 56, 55, 54, 58, 62, 63, 64, 59, 57, 56, 55, 54, 53, 17, 20, 14, 10, 12, 26, 34, 19, 18, 25),
#               
#               
#               bergmanG = 27,  
#               oligoG = c(21, 22, 23, 24, 35), 
#               astrocytesG = c(31, 32), 
#               glials = c(27,21,22,35,31,32, 23, 24)
# )


# realGroups = vector(mode='list',length = length(groups))
# names(realGroups) = names(groups)
# for (i in 1:length(groups)){
#   realGroups[[i]] = which(design$originalIndex %in% groups[[i]])
# }


# grouping using ourNaming fields in design
groupNames = with(design, unique(c(as.character(ourNamingGeneral), as.character(ourNaming), as.character(cellKind))))
groupNames = trimNAs(groupNames)
realGroups = vector(mode = 'list', length = length(groupNames))
names(realGroups) = groupNames
for (i in 1:length(groupNames)){
  realGroups[[i]] = which(!is.na(design$ourNaming) & 
                            (design$ourNaming == groupNames[i] | design$ourNamingGeneral == groupNames[i] | design$cellKind == groupNames[i]))
}



progy = txtProgressBar(min=1,max=length(realGroups)*nrow(exprData), style=3)
step = 0
silly = vector(mode = 'double', length = length(realGroups)*nrow(exprData))
dim(silly) = c(length(realGroups), nrow(exprData))
#for single group trial
#progy = txtProgressBar(min=1,max=nrow(exprData))
for (i in 1:length(realGroups)){
  for (j in 1:nrow(exprData)){
    members = realGroups [[i]]
    clustering = as.integer(rep(1,nrow(design))*(1:nrow(design) %in% members)+1)
    data = t(exprData[j,])
    cluster = list(clustering = clustering, data = data)    
    silo = silhouette(cluster,dist(data))
    silly[i, j] = mean(silo[,3])
    setTxtProgressBar(progy,step)
    step = step+1
  }
}
close(progy)


allSilhouette = cbind(geneData,t(silly))
colnames(allSilhouette) = c(colnames(allSilhouette)[1:3], names(realGroups))
write.csv(allSilhouette, quote = T, row.names = F, col.names = T, 'results/allSilhouette')
save(allSilhouette,file = 'results/allSilhouette.rData')
 for (i in 1:length(realGroups)){
   filename = dpaste('results/', names(realGroups)[i])
   geneList = geneData$Gene.Symbol[order(silly[i,], decreasing = T)]
   silhouettes = silly[i,order(silly[i,], decreasing = T)]

  
  
  df = data.frame (genes = geneList, silhouettes = silhouettes)
  write.table(df, quote = F, row.names = F, col.names= T, filename )
  
}
  
play(w)
# 
# geneData$Gene.Symbol[order(silly[1,], decreasing = T)][1:50]
# 
# plotGene('Pank1', realGroups[[1]])
# plotGene('Gad1', realGroups[[4]])
# 
# j=giveGeneIndex('A630007B06Rik')
# j=giveGeneIndex('Gad1')
# 
# toFile = as.data.frame(t(silly))
# colnames(toFile) = names(realGroups)
# 
# 
# for (i in 1:length(realGroups))

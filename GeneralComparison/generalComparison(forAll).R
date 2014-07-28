require(gdata)
#startup
######################
dpaste = function (...){
  paste(..., sep='')
}
giveName = function(numList){
  namesANDindex = data.frame(geneNames = geneNames, geneIndex = geneIndex)
  nAi = namesANDindex[order(namesANDindex$geneIndex),]
  
  matches = which(nAi$geneIndex %in% numList)
  gNames = nAi$geneNames[matches]
  indexes = which(numList %in% nAi$geneIndex)
  
  return(data.frame(geneNames = gNames, index = indexes ))
}

giveIndex = function(geneList){
  indexes = vector(length = length(geneList))
  for (i in 1:length(geneList)){
  indexes[i] = geneIndex[which(geneNames %in% geneList[i])]
  }
  return(indexes)
}

des = function(sb){
  return(sb[, (ncol(sb)-ncol(transcriptomeSources)+2):ncol(sb)])
}

dat = function(sb){
  return(sb[, 1:(ncol(sb)-ncol(transcriptomeSources))])
}

foldChange = function (group1, group2, f = 10){
  fold = group1 / group2
  chosen =  which(fold <= (1/f) | fold >= f)
  return(
    data.frame(index = chosen, foldChange = fold[chosen])
  )
}
geneIndex = as.numeric(read.csv('colsToUse', header=FALSE))
geneNames = as.character(read.table('usedGenes')$V1)

transcriptomeSources = read.xls('journal.pone.0016493.s010.xls')#perl='C:/Perl64/bin/perl.exe')
transcriptomeSources$SampleNo=1:64
colnames(transcriptomeSources)=c('officialName','cellType','anatomical','age','method','isolation','amplification','RNAamount','platform','reference','sampleNo')

transcriptomeSources$age = gsub('~', '', transcriptomeSources$age)
transcriptomeSources$age = gsub('P', '', transcriptomeSources$age)
transcriptomeSources$age = gsub('7-8', '7.5', transcriptomeSources$age)
transcriptomeSources$age[grepl('(precise age not given)',transcriptomeSources$age)] = 60
transcriptomeSources$age = as.numeric(transcriptomeSources$age)

transcriptomeSources$anatomical = gsub('Corpus Striatum', 'Striatum', transcriptomeSources$anatomical)
transcriptomeSources$anatomical = gsub('Motor Cortex', 'Motor cortex', transcriptomeSources$anatomical)
transcriptomeSources$anatomical = gsub('Layer 5A Cortex', 'Cortex', transcriptomeSources$anatomical)
transcriptomeSources$anatomical = gsub('Layer 5B Cortex', 'Cortex', transcriptomeSources$anatomical)
transcriptomeSources$anatomical = gsub('Layer 6 Cortex', 'Cortex', transcriptomeSources$anatomical)

transcriptomeSources$cellType[grepl('^Neurons', transcriptomeSources$cellType)] = 'Mixed Neurons'
transcriptomeSources$cellType=gsub(', P.*?$','',transcriptomeSources$cellType)
allDataPre = read.csv('allData', header = F)
allData=2^allDataPre
allDataD = cbind(allData, transcriptomeSources)
###################################################################
# transcriptomeSources
# geneNames
# geneIndex
# allData
# usedData
# giveName (f)
# giveIndex (f)
# des (f)
# dat (f)
# foldChange (f)


#Group Select
################

#SEND IN INDEXES Of allDataD. conamination removal should be done during grouping
#two lines of groups are there to allow combining groups
#look for data here using the parameters
#des(allDataD)

#manually looking for contamination with markers
giveIndex('')
allDataD[c('')]


#groups = list(cortexOligo = 35, cortexAstro = c(31,32), cortexPyramidal = c(7,8,40,43,44), cortexMix = c(9,26,34), cortexInter= 14)
#groups = c(groups, with(groups, list(neurons = c(cortexPyramidal, cortexMix, cortexInter), glials = c(cortexOligo, cortexAstro))))

#manually looking for contamination
markers = c('Eno2','Aqp4','Adamts4','Gad1','Slc32a1','Th')
markerI = giveIndex(markers)
markerCol = dpaste('V',markerI)
a=allDataD[,c(markerCol,'cellType','age')]
colnames(a)= c('Eno2-Neur', 'Aqp4-astro', 'Adamts4-oligo','Gad1-gaba','Slc32a1-gaba', 'Th-dopa','cellType','age')

des(allDataD[-c(1, 6,9, 21,  22, 24, 28, 29),][allDataD$age[-c(1, 6,9, 21,  22, 24, 28, 29)]>=14,])
des(allDataD[-c(1, 6,9, 21,  22, 24, 28, 29,2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 26, 34, 40, 43, 44, 45:50, 51, 52, 53,  54:64),][allDataD$age[-c(1, 6,9, 21,  22, 24, 28, 29,2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 26, 34, 40, 43, 44, 45:50, 51, 52, 53,  54:64)]>=14,])

groups = list(neurons = c(2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 26, 34, 40, 43, 44, 45:50, 51, 52, 53,  54:64), glial = c(23,27,31,32,35), oligo = c(23, 35), astro = c(31,32))
              
              
#make 'all' do do every comparison
#comparisons = 'all'
#write as is

comparisons =
"neurons, glial
neurons, oligo
neurons, astro
oligo, astro"

#comparisons =
#"neurons, glials
#cortexPyramidal, cortexInter
#cortexAstro, cortexOligo"

#write names in the order in groups. if no intersection. manual removal is required

intersections = ''

#intersections = 
#"neurons-glials, cortexOligo-cortexAstro"


foldTreshold = 10
#########################


#Analysis
#############
conn = textConnection(comparisons)
comp = read.csv(conn, header = F, strip.white = T)

groupNames = names(groups)

groupAverages = list()
for (i in groups){
  groupAverage = apply(dat(allDataD[i,]), 2, mean)
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
  temp = giveName(fChange$index)
  fChangeKnown = data.frame(geneNames = temp$geneNames, geneFoldChange = fChange$foldChange[temp$index])
  fChangeKnown = fChangeKnown[order(fChangeKnown$geneFoldChange, decreasing=T) ,]
  print(fileName)
  print(i)
  write.table(fChangeKnown, quote = F, row.names = F, col.names = F, fileName)
}


#intersection between groups
conn = textConnection(intersections)
comp = read.csv(conn, header = F, strip.white = T)
for (i in 1:nrow(comp)){
  a = read.table(dpaste('Results/',comp[i,1]))
  b = read.table(dpaste('Results/',comp[i,2]))
  fileName = dpaste('Intersections/', comp[i,1],' & ', comp[i,2])
  write.table(intersect(a$V1,b$V1), quote = F, row.names = F, col.names = F, fileName)
}


#########################



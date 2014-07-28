#last execution: age 14<= with undesirables discluded

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

#disclusion is done below when comparing it to rest.
#discludeSamples = which(design$age<=14)
#exprData = exprData[, -discludeSamples]
#design = design[-discludeSamples, ]
exprData = t(exprData)
#########################################


#define groups
#####################
#all adults
fromIndexes = T
groups = list(pyramidalN = c(40, 43, 44,45, 46, 47, 48, 49),     
              cholinergicN = c(13, 11), 
              spinyN = c(15, 16), 
              gabaN = c(58, 62, 63, 64, 59, 56, 55, 54),
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
              mixedN = c(26, 34),
              basketN = 19,
              ubcN = 18,
              purkinjeN = 25,
              neurons = c(40, 43, 44, 45, 46, 47, 48, 49, 13, 11, 15, 16, 58, 62, 63, 64, 59, 56, 55, 54, 58, 62, 63, 64, 59, 57, 56, 55, 54, 53, 17, 20, 14, 10, 12, 26, 34, 19, 18, 25),
              
              
              bergmanG = 27,  
              oligoG = c(21, 22, 23, 24, 35), 
              astrocytesG = c(31, 32), 
              glials = c(27,21,22,35,31,32, 23, 24)
)

comparisons =
"all
rest"

intersections = ''
foldTreshold = 10
#analysis
############################################
#convert original index to index in 
if (fromIndexes == T){
  realGroups = vector(mode='list',length = length(groups))
  names(realGroups) = names(groups)
  for (i in 1:length(groups)){
    realGroups[[i]] = which(design$originalIndex %in% groups[[i]])
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
                            print('fuu')
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
                              print('fuu')
                              return(exprData[i,])
                            })
    
    restAverage = tryCatch({apply(exprData[-c(i,which(design$age<14), which(design$cellType2 == 'Astroglia'),which(design$originalIndex %in% c(6,7,8))),], 2, median)},
                           error = function(cond){
                             print('fuuiuu')
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

log(exprData[c('GSM63033',
               'GSM63034',
               'GSM63035',
               'GSM63039',
               'GSM63040',
               'GSM63041',
               'GSM63036',
               'GSM63037',
               'GSM63038',
               'GSM63045',
               'GSM63046',
               'GSM63047',
               'GSM63042',
               'GSM63043',
               'GSM63044'
               ),c(2753,11297)])/log(2)
log(exprData[,c(2753,11297)])/log(2)

'GSM63048',
'GSM63049',
'GSM63050'


log(exprData[c('GSM63048',
               'GSM63049',
               'GSM63050'
),c(2753,11297)])/log(2)







'GSM337773',
'GSM337774',
'GSM337775',
'GSM337779',
'GSM337780',
'GSM337781',
'GSM337785',
'GSM337786',
'GSM337787'


log(exprData[c('GSM337773',
               'GSM337774',
               'GSM337775',
               'GSM337779',
               'GSM337780',
               'GSM337781',
               'GSM337785',
               'GSM337786',
               'GSM337787'
),5105])/log(2)


log(exprData[
,5105])/log(2)

GSM63033
GSM63034
GSM63035
GSM63039
GSM63040
GSM63041
GSM63036
GSM63037
GSM63038
GSM63045
GSM63046
GSM63047
GSM63042
GSM63043
GSM63044



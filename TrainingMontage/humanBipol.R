#initial loading
load('GSE12649_GSE5388_expression.RData')


source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')
insist(R.matlab)
insist(biomaRt)
insist(ggplot2)
insist(plyr)


Matlab$startServer()
matlab = Matlab()
isOpen = open(matlab)


predCon = function(samples, types){
  selectGroupData = types
  if (!isOpen)
    throw("MATLAB server is not running")
  
  random = samples
  selectGroupData = types
  
  predictions = matrix(0, nrow = ncol(selectGroupData), ncol = ncol(random))
  
  MTakenGeneByType = as.matrix(selectGroupData)
  setVariable(matlab, MTakenGeneByType=MTakenGeneByType)
  
  diagNegTypes = diag(x = -1, nrow = ncol(selectGroupData), ncol = ncol(selectGroupData))
  refZeroTypes = matrix(0, ncol(selectGroupData), ncol = 1)
  setVariable(matlab, refZeroTypes=refZeroTypes)
  setVariable(matlab, diagNegTypes=diagNegTypes)
  
  
  for (i in 1:ncol(random)){
    voxelData = as.matrix(random[, i])
    setVariable(matlab, voxelData=voxelData)
    evaluate(matlab, 'fitVoxelsToTypesBis= lsqlin( MTakenGeneByType, voxelData, diagNegTypes, refZeroTypes );')
    predictions [, i] =  getVariable(matlab, 'fitVoxelsToTypesBis')[[1]]
    
  } 
  
  predictions[abs(predictions) <10e-15] = 0
  return(predictions)
}

giveGI = function(daGene){
  return(match(daGene, geneData$Gene.Symbol))  
}

# loading files ----------------


design = read.table('C:/Users/Ogan/Dropbox/Rotation 3/Data/normalizedDesign.csv',header=T,sep='\t')
design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)
allDataPre = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/Data/mostVariableQuantileNormalized', header = T)
geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData = exprData
design = design[match(colnames(exprData),design$sampleName,),]
rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<5)

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]
exprData = 2^exprData

# Get the list and indexes of final selected genes. Also create the selected matrix --------------------
filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/ourNamingGeneral",include.dirs = FALSE)
fileContents = lapply(paste('C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/ourNamingGeneral/', filenames, sep = ''), read.table)
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
geneList = puristList


length(unique(geneList))
geneList = unique(geneList)
#geneList = c('Mobp', 'Slc1a3', 'Slc1a2', 'Aqp4', 'Mag', 'Mog', 'Tac1', 'Neurod6')
geneIndexes = trimNAs(giveGI(unlist(geneList))) #note that the data is already trimmed for low expressed ones so indexes arent from original exprData

# also remove the non desirables
selectGeneData = geneData[geneIndexes, ]
#old criteria
#selectExprData = exprData[geneIndexes, (design$age>=14)&(!design$cellType2=='Astroglia')&(!design$originalIndex %in% c(6,7,8))]
selectGeneData = geneData[geneIndexes, ]
selectExprData = exprData[geneIndexes, !is.na(design$ourNamingGeneral)]
selectDes = design[!is.na(design$ourNamingGeneral),]


selectGroupData = matrix(0, nrow = nrow(selectExprData), ncol = length(unique(selectDes$ourNamingGeneral)))
groupNames = unique(selectDes$ourNamingGeneral)
colnames(selectGroupData) = groupNames
for (i in 1:length(groupNames)){
  groupIndex = selectDes$ourNamingGeneral == groupNames[i]
  groupAverage =  tryCatch({apply(selectExprData[, groupIndex ], 1, mean)},
                           error = function(cond){
                             print('fuu')
                             return(selectExprData[, groupIndex ])
                           })
  selectGroupData[, i] = groupAverage 
}  


#loading up the bipolar data. finding mouse homologues to match expressions 

bipolarGeneData1 = aned_high_GSE12649[, 1:3]
bipolarGeneData2 = aned_high_GSE5388[, 1:3]

bipolarExpr1 = 2^aned_high_GSE12649[,4:ncol(aned_high_GSE12649)]
bipolarExpr2 = 2^aned_high_GSE5388[,4:ncol(aned_high_GSE5388)]

orthoInfo = read.csv('HG-U133A.na33.ortholog.csv')
orthoInfo = orthoInfo[orthoInfo$Ortholog.Array == 'MOE430A',1:5 ]



usePairs1 = orthoInfo[{tolower(orthoInfo$Ortholog.Probe.Set) %in% tolower(selectGeneData$Probe)}&
{tolower(orthoInfo$Probe.Set.ID) %in% tolower(bipolarGeneData1$Probe)}
,c(1,3)]
dim(aned_high_GSE12649)

usePairs1 = usePairs1[-which(duplicated(usePairs1[,1])),]
usePairs1 = usePairs1[-which(duplicated(usePairs1[,2])),]
mouseRegionData = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/mouseRegions/mouseRegions')
mouseRegionExpr = mouseRegionData[,4:ncol(mouseRegionData)]
mouseRegionGenes = mouseRegionData[,1:3]

useMouseGeneData = exprData
useMouseExp =geneData

useMouseExp = mouseRegionExpr[ - which( !mouseRegionGenes$Probe %in% selectGeneData$Probe ), ]
useMouseGeneData = mouseRegionGenes[- which( !mouseRegionGenes$Probe %in% selectGeneData$Probe ), ]

useMouseExp = useMouseExp[match(selectGeneData$Probe, useMouseGeneData$Probe), ] 
useMouseGeneData = useMouseGeneData[match(selectGeneData$Probe, useMouseGeneData$Probe), ] 


selectGroupData = matrix(0, nrow = nrow(selectExprData), ncol = length(unique(selectDes$ourNamingGeneral)))
groupNames = unique(selectDes$ourNamingGeneral)
colnames(selectGroupData) = groupNames
for (i in 1:length(groupNames)){
  groupIndex = selectDes$ourNamingGeneral == groupNames[i]
  groupAverage =  tryCatch({apply(selectExprData[, groupIndex ], 1, mean)},
                           error = function(cond){
                             print('fuu')
                             return(selectExprData[, groupIndex ])
                           })
  selectGroupData[, i] = groupAverage 
}  


contribs = predCon(samples=useMouseExp, types=selectGroupData)
colnames(contribs) = Samples_GSE12649[1,]
rownames(contribs) = groupNames
write.csv(contribs, file = 'contribs_GSE12649')




useBipolarGeneData1 = bipolarGeneData1[tolower(bipolarGeneData1$Probe) %in% tolower(usePairs1$Probe.Set.ID),]
useBipolarExpr1 = bipolarExpr1[tolower(bipolarGeneData1$Probe) %in% tolower(usePairs1$Probe.Set.ID),]

useSelectGeneData1 = selectGeneData[tolower(selectGeneData$Probe) %in% tolower(usePairs1$Ortholog.Probe.Set),]
useSelectExprData1 = selectExprData[tolower(selectGeneData$Probe) %in% tolower(usePairs1$Ortholog.Probe.Set),]


useSelectGeneData1 = useSelectGeneData1[match(tolower(usePairs1$Ortholog.Probe.Set), tolower(useSelectGeneData1$Probe)),]
useBipolarGeneData1 = useBipolarGeneData1[match(tolower(usePairs1$Probe.Set.ID), tolower(useBipolarGeneData1$Probe)),]

useBipolarExpr1 = useBipolarExpr1[match(tolower(usePairs1$Probe.Set.ID), tolower(useBipolarGeneData1$Probe)),]
useSelectExprData1 = useSelectExprData1[match(tolower(usePairs1$Ortholog.Probe.Set), tolower(useSelectGeneData1$Probe)),]


selectGroupData = matrix(0, nrow = nrow(useSelectExprData1), ncol = length(unique(selectDes$ourNamingGeneral)))
groupNames = unique(selectDes$ourNamingGeneral)
colnames(selectGroupData) = groupNames
for (i in 1:length(groupNames)){
  groupIndex = selectDes$ourNamingGeneral == groupNames[i]
  groupAverage =  tryCatch({apply(useSelectExprData1[, groupIndex ], 1, mean)},
                           error = function(cond){
                             print('fuu')
                             return(useSelectExprData1[, groupIndex ])
                           })
  selectGroupData[, i] = groupAverage 
}  


contribs = predCon(samples=useBipolarExpr1, types=selectGroupData)
colnames(contribs) = Samples_GSE12649[1,]
rownames(contribs) = groupNames
write.csv(contribs, file = 'contribs_GSE12649')

#remove non cortex stuff
selectGroupData = matrix(0, nrow = nrow(useSelectExprData1), ncol = length(unique(selectDes$ourNamingGeneral[selectDes$anatomical %in% 'Cortex'])))
selectGroupData = matrix(0, nrow = nrow(useSelectExprData1), ncol = length(unique(selectDes$ourNamingGeneral[selectDes$ourNamingGeneral %in% c('Oligo','Astro','MotorCholin','Inter','Pyramidal', 'Inter','Gaba','Microglia')])))
cortexDes = selectDes[selectDes$anatomical %in% 'Cortex',]
cortexDes = selectDes[selectDes$ourNamingGeneral %in% c('Oligo','Astro','MotorCholin','Inter','Pyramidal', 'Inter','Gaba','Microglia'),]
groupNames = unique(cortexDes$ourNamingGeneral)
colnames(selectGroupData) = groupNames
for (i in 1:length(groupNames)){
  groupIndex = selectDes$ourNamingGeneral == groupNames[i]
  groupAverage =  tryCatch({apply(useSelectExprData1[, groupIndex ], 1, mean)},
                           error = function(cond){
                             print('fuu')
                             return(useSelectExprData1[, groupIndex ])
                           })
  selectGroupData[, i] = groupAverage 
}  


contribsNoCortex = predCon(samples=useBipolarExpr1, types=selectGroupData)
colnames(contribsNoCortex) = Samples_GSE12649[1,]
rownames(contribsNoCortex) = groupNames
write.csv(contribs, file = 'Cortex_contribs_GSE12649')









# second bipolar data --------------
usePairs2 = orthoInfo[{tolower(orthoInfo$Ortholog.Probe.Set) %in% tolower(selectGeneData$Probe)} & 
{tolower(orthoInfo$Probe.Set.ID) %in% tolower(bipolarGeneData2$Probe)} 
,c(1,3)]

usePairs2 = usePairs2[-which(duplicated(usePairs2[,1])),]
usePairs2 = usePairs2[-which(duplicated(usePairs2[,2])),]


useBipolarGeneData2 = bipolarGeneData2[tolower(bipolarGeneData2$Probe) %in% tolower(usePairs2$Probe.Set.ID),]
useBipolarExpr2 = bipolarExpr2[tolower(bipolarGeneData2$Probe) %in% tolower(usePairs2$Probe.Set.ID),]

useSelectGeneData2 = selectGeneData[tolower(selectGeneData$Probe) %in% tolower(usePairs2$Ortholog.Probe.Set),]
useSelectExprData2 = selectExprData[tolower(selectGeneData$Probe) %in% tolower(usePairs2$Ortholog.Probe.Set),]


useSelectGeneData2 = useSelectGeneData2[match(tolower(usePairs2$Ortholog.Probe.Set), tolower(useSelectGeneData2$Probe)),]
useBipolarGeneData2 = useBipolarGeneData2[match(tolower(usePairs2$Probe.Set.ID), tolower(useBipolarGeneData2$Probe)),]

useBipolarExpr2 = useBipolarExpr2[match(tolower(usePairs2$Probe.Set.ID), tolower(useBipolarGeneData2$Probe)),]
useSelectExprData2 = useSelectExprData2[match(tolower(usePairs2$Ortholog.Probe.Set), tolower(useSelectGeneData2$Probe)),]



selectGroupData = matrix(0, nrow = nrow(useSelectExprData2), ncol = length(unique(selectDes$ourNamingGeneral)))
groupNames = unique(selectDes$ourNamingGeneral)
colnames(selectGroupData) = groupNames
for (i in 1:length(groupNames)){
  groupIndex = selectDes$ourNamingGeneral == groupNames[i]
  groupAverage =  tryCatch({apply(useSelectExprData2[, groupIndex ], 1, mean)},
                           error = function(cond){
                             print('fuu')
                             return(useSelectExprData2[, groupIndex ])
                           })
  selectGroupData[, i] = groupAverage 
}  


contribs = predCon(samples=useBipolarExpr2, types=selectGroupData)
colnames(contribs) = Samples_GSE5388[1,]
rownames(contribs) = groupNames
write.csv(contribs, file = 'contribs_GSE5388')

#remove non cortex stuff
selectGroupData = matrix(0, nrow = nrow(useSelectExprData2), ncol = length(unique(selectDes$ourNamingGeneral[selectDes$ourNamingGeneral %in% c('Oligo','Astro','MotorCholin','Inter','Pyramidal', 'Inter','Gaba','Microglia')])))
cortexDes = selectDes[selectDes$anatomical %in% 'Cortex',]
cortexDes = selectDes[selectDes$ourNamingGeneral %in% c('Oligo','Astro','MotorCholin','Inter','Pyramidal', 'Inter','Gaba','Microglia'),]
groupNames = unique(cortexDes$ourNamingGeneral)
colnames(selectGroupData) = groupNames
for (i in 1:length(groupNames)){
  groupIndex = selectDes$ourNamingGeneral == groupNames[i]
  groupAverage =  tryCatch({apply(useSelectExprData2[, groupIndex ], 1, mean)},
                           error = function(cond){
                             print('fuu')
                             return(useSelectExprData2[, groupIndex ])
                           })
  selectGroupData[, i] = groupAverage 
}  


contribsNoCortex = predCon(samples=useBipolarExpr2, types=selectGroupData)
colnames(contribsNoCortex) = Samples_GSE5388[1,]
rownames(contribsNoCortex) = groupNames
write.csv(contribsNoCortex, file = 'Cortex_contribs_GSE5388')


#
#useless bullpoo ---------------
mmart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#mmusculus_homolog_orthology_type
ensemblIDs = getBM(attributes = c('mgi_symbol', 'ensembl_gene_id'),filters ='mgi_symbol' ,
                   values = as.character(selectGeneData$Gene.Symbol),
                   mart = mmart)
orthologues = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene', 'mmusculus_paralog_perc_id_r1'), filters = "ensembl_gene_id", 
                    values = ensemblIDs$ensembl_gene_id,
                    mart = mmart)


hmart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

humanOrtho = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id',
                   values = orthologues$hsapiens_homolog_ensembl_gene,
                   mart = hmart)


homoFrame = cbind(orthologues, 
                  ensemblIDs[ match(orthologues$ensembl_gene_id, ensemblIDs$ensembl_gene_id), ], 
                  humanOrtho[ match(orthologues$hsapiens_homolog_ensembl_gene, humanOrtho$ensembl_gene_id) ,])

relevantFrame = unique(homoFrame[,c('mgi_symbol', 'hgnc_symbol','mmusculus_paralog_perc_id_r1' )])

relevantFrame = relevantFrame[!is.na(relevantFrame$hgnc_symbol) & !relevantFrame$hgnc_symbol == '', ]

relevantFrame2 = frameCollapse(relevantFrame, c('mgi_symbol', 'hgnc_symbol'), max)


bipolarGeneData1[which(relevantFrame2 %in% bipolarGeneData1$name), ]
bipolarGeneData1[which(relevantFrame2 %in% bipolarGeneData2$name), ]




aggregate(relevantFrame,by=list(key=relevantFrame$x),FUN=function(x){x[1]})


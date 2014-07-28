insist = function(name){
  name = substitute(name) 
  name = as.character(name) 
  if (!require(name, character.only = T)) {
    install.packages(name)
    library(name, character.only = T)
  }
}

insist(R.matlab)
insist(biomaRt)
insist(ggplot2)

Matlab$startServer()
matlab = Matlab()
isOpen = open(matlab)


# Functions for the aftermath ---------

trimNAs = function(aVector) {
  return(aVector[!is.na(aVector)])
}

checkPred = function(truth, prediction){
  #make all sums 1 so only proportions will effect the results
  sumPre = apply(prediction, 2, sum)
  sumTru = apply(truth, 2, sum)
  objPre = prediction / repRow(sumPre, nrow(prediction))
  objTru = truth / repRow(sumTru, nrow(truth))
  return(diag(cor(truth, prediction, method = 'pearson')))
}

checkPred2 = function(truth,prediction){
  sumPre = apply(prediction, 2, sum)
  sumTru = apply(truth, 2, sum)
  objPre = prediction / repRow(sumPre, nrow(prediction))
  objTru = truth / repRow(sumTru, nrow(truth))
  return(diag(cor(t(truth), t(prediction), method = 'pearson')))
}

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
  
  predictions[predictions <10e-15] = 0
  return(predictions)
}

#this function selects a sample to add as noise according to noiseMax and if maxMix = 0, it outputs just the noisyOne + noiseMax noise as in the first form of the function
generateSamplesSetNoise = function(noiseMax, randomNo = 30, maxMix = 5,mode = 'constant', noisyOne){
  random =  matrix(0, nrow = length(geneIndexes), ncol = randomNo)
  contMatrix = matrix (0, nrow = ncol(selectGroupData), ncol = randomNo)
  
  if (mode == 'constant'){
    noise = rep(noiseMax,randomNo)
    noOfSamples = maxMix
  } else if (mode == 'random'){
    noise = runif(randomNo, 0, noiseMax)
    noOfSamples = floor(runif(1,2,maxMix))
  }
  for (i in 1:randomNo){
    #select how many samples to mix and which samples to mix
    sampleWeights = runif(noOfSamples, 0, 1)
    whichSamples = sample(ncol(selectGroupData) -1, noOfSamples )
    # prevent noisyOne from being selected. throw it to the back. if its the last one it will never be equal anyway.
    
    whichSamples[whichSamples == noisyOne] = noOfSamples
    
    #adding noise based on one of the samples but with shuffled titles
    
    if (maxMix == 0){
      noice = selectGroupData [sample(nrow(selectGroupData), nrow(selectGroupData)), sample(ncol(selectGroupData), 1)] * noise[i]
      random[, i]= selectGroupData [, noisyOne] * (1 - noise) + noice * noise
    }
    
    addNoise = selectGroupData [, noisyOne] * noise[i]
    contWrite = sampleWeights/sum(sampleWeights)
    sampleWeights = {sampleWeights-sampleWeights*noise[i]}/sum(sampleWeights)
    contMatrix [whichSamples, i] = contWrite
    sampleMult = lapply(sampleWeights, rep , nrow(selectGroupData))
    sampleMult = unlist(sampleMult)
    if (!maxMix == 0){
    dim(sampleMult) = c(nrow(selectGroupData),noOfSamples)
    random[, i] = apply(selectGroupData[, whichSamples] * sampleMult , 1, sum) +addNoise
    }
  }
  outlist = list (random = random, noise = noise, contMat = contMatrix[-noisyOne, ], sourceSemps = selectGroupData[, -noisyOne])
  return(outlist)
}


generateSamples = function(noiseMax, randomNo = 30, maxMix = 5,mode = 'constant', noisyGenes = NA, geneNoiseMean = 0, geneNoiseSD = 1, mixSamp = NA){
  #noiseMax = 0.01
  #randomNo = 20
  #maxMix = ncol(selectGroupData)
  random =  matrix(0, nrow = length(geneIndexes), ncol = randomNo)
  contMatrix = matrix (0, nrow = ncol(selectGroupData), ncol = randomNo)
  
  if (mode == 'constant'){
    noise = rep(noiseMax,randomNo)
    noOfSamples = maxMix
  } else if (mode == 'random'){
    noise = runif(randomNo, 0, noiseMax)
    noOfSamples = floor(runif(1,2,maxMix))
  }
  for (i in 1:randomNo){
    #select how many samples to mix and which samples to mix
    sampleWeights = runif(noOfSamples, 0, 1)
    if (is.na(mixSamp)){
    whichSamples = sample(ncol(selectGroupData), noOfSamples )
    #adding noise based on one of the samples but with shuffled titles
    } else {
      whichSamples = mixSamp
      sampleWeights = runif(length(mixSamp),0,1)
      noOfSamples = (length(mixSamp))
    }
    
    
    addNoise = selectGroupData [sample(nrow(selectGroupData), nrow(selectGroupData)), sample(ncol(selectGroupData), 1)] * noise[i]
    contWrite = sampleWeights/sum(sampleWeights)
    sampleWeights = {sampleWeights-sampleWeights*noise[i]}/sum(sampleWeights)
    contMatrix [whichSamples, i] = contWrite
    sampleMult = lapply(sampleWeights, rep , nrow(selectGroupData))
    sampleMult = unlist(sampleMult)
    dim(sampleMult) = c(nrow(selectGroupData),noOfSamples)
    random[, i] = apply(selectGroupData[, whichSamples] * sampleMult , 1, sum) +addNoise
    if (!any(is.na(noisyGenes))){
      random[noisyGenes, i] = random[noisyGenes, i] + rnorm(length(noisyGenes), geneNoiseMean, geneNoiseSD)
    }
    
    
  }
  outlist = list (random = random, noise = noise, contMat = contMatrix)
  return(outlist)
}

repIndiv = function (aVector, n){
  output = vector(length = length(aVector) * n)
  step = 1
  for (i in aVector){
    output[(step * n - n + 1):(n * step)] = rep(i, n)
    step = step + 1
  }
  return(output)
}

repRow = function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

repCol = function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

trimElement = function (aVector,e){
  return(aVector[aVector!=e])
}

giveGI = function(daGene){
  return(match(daGene, geneData$Gene.Symbol))  
}


# As always open the bloody files. ourNaming field added to the design file ----------------
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


# Get the list and indexes of final selected genes. Also create the selected matrix --------------------
filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/cheatsyClust/Results/min8FinerGroups",include.dirs = FALSE)
fileContents = lapply(paste('C:/Users/Ogan/Dropbox/Rotation 3/cheatsyClust/Results/min8FinerGroups/', filenames, sep = ''), read.table)
geneList = vector(mode = 'character', length = 68000)
step = 1
for (i in fileContents){
  length(i$V1)
  geneList[step:(step + length(i$V1)-1)] = as.character(i$V1)
  step = step+length(i$V1)
}

geneList = trimElement(geneList,'')


#replace it with their gene list comment out when returning to your own sample
#allDataPre = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/Data/mostVariableQuantileNormalized', header = T)
#geneData = allDataPre[,1:3]
#exprData = allDataPre[,4:ncol(allDataPre)]
#geneList2= as.character(read.table('C:/Users/Ogan/Dropbox/Rotation 3/Data/usedGenes')$V1)

#just to compare two two lists. not related to the rest of the code.
#insist(VennDiagram)
#Venn = list(ours = geneList, theirs= geneList2)
#venn.diagram(Venn, fill = c("red", "green"),filename='bok.png')


length(unique(geneList))
geneList = unique(geneList)
geneIndexes = trimNAs(giveGI(geneList)) #note that the data is already trimmed for low expressed ones so indexes arent from original exprData

# also remove the non desirables
selectGeneData = geneData[geneIndexes, ]
#old criteria
#selectExprData = exprData[geneIndexes, (design$age>=14)&(!design$cellType2=='Astroglia')&(!design$originalIndex %in% c(6,7,8))]
selectExprData = exprData[geneIndexes, !is.na(design$ourNaming)]
selectDes = design[!is.na(design$ourNaming),]


#group replicates together by... mean? ------------------------
selectGroupData = matrix(0, nrow = nrow(selectExprData), ncol = length(unique(selectDes$ourNaming)))
groupNames = unique(selectDes$ourNaming)
colnames(selectGroupData) = groupNames
for (i in 1:length(groupNames)){
  groupIndex = selectDes$ourNaming == groupNames[i]
  groupAverage =  tryCatch({apply(selectExprData[, groupIndex ], 1, mean)},
                           error = function(cond){
                             print('fuu')
                             return(selectExprData[, groupIndex ])
                           })
  selectGroupData[, i] = groupAverage 
}  


#Generate random samples. like... very random. mix and add noise -------------------------------------

random.0 = generateSamples(0)
random.01 = generateSamples(0.01)
random.2 = generateSamples(0.2)
random.35 = generateSamples(0.35)
random.5 = generateSamples(0.5)

pred.0 = predCon(random.0$random,selectGroupData)
pred.01 = predCon(random.01$random,selectGroupData)
pred.2 = predCon(random.2$random,selectGroupData)
pred.35 = predCon(random.35$random,selectGroupData)
pred.5 = predCon(random.5$random,selectGroupData)

fr = data.frame(Pearson = c(checkPred(random.0$contMat, pred.0),
                            checkPred(random.01$contMat, pred.01),
                            checkPred(random.2$contMat, pred.2),
                            checkPred(random.35$contMat, pred.35),
                            checkPred(random.5$contMat, pred.5)
), 
RandomNoiseRatio = as.factor(repIndiv(c(0,0.01,0.2,0.35,0.5),length(checkPred(random.0$contMat, pred.0))))
)

insist(ggplot2)

(ggplot(fr, aes(x = RandomNoiseRatio, y = Pearson))+geom_boxplot()+ ggtitle("Simple Noise") + theme(text = element_text(size=20)))


fr = data.frame(Pearson = c(checkPred2(random.0$contMat, pred.0),
                            checkPred2(random.01$contMat, pred.01),
                            checkPred2(random.2$contMat, pred.2),
                            checkPred2(random.35$contMat, pred.35),
                            checkPred2(random.5$contMat, pred.5)
                            ), 
                RandomNoiseRatio = as.factor(repIndiv(c(0,0.01,0.2,0.35,0.5),length(checkPred2(random.0$contMat, pred.0))))
                )
(ggplot(fr, aes(x = RandomNoiseRatio, y = Pearson))+geom_boxplot()+ ggtitle("Simple Noise. Correlation for individual samples")+theme_bw() + theme(text = element_text(size=20)))


# one sample as noise ---------------------
#function is identical to the previous one except it takes one of the samples as a constant noise. it also adds the new known sample matrix
#to the output list


selectedNoise.2 = vector(mode = 'list', length= ncol(selectGroupData))
for (i in 1:ncol(selectGroupData)){
  selectedNoise.2[[i]] = generateSamplesSetNoise(0.2, noisyOne = i)
}


selectedNoise.35 = vector(mode = 'list', length= ncol(selectGroupData))
for (i in 1:ncol(selectGroupData)){
  selectedNoise.35[[i]] = generateSamplesSetNoise(0.35, noisyOne = i)
}

selectedNoise.5 = vector(mode = 'list', length= ncol(selectGroupData))
for (i in 1:ncol(selectGroupData)){
  selectedNoise.5[[i]] = generateSamplesSetNoise(0.5, noisyOne = i)
}




predSele.2 = vector (mode = 'list', length = length(selectedNoise.2))
predSele.35 = vector (mode = 'list', length = length(selectedNoise.2))
predSele.5 = vector (mode = 'list', length = length(selectedNoise.5))

for (i in 1:length(selectedNoise.2)){
  predSele.2[[i]] = predCon (selectedNoise.2[[i]]$random, selectedNoise.2[[i]]$sourceSemps)
  predSele.35[[i]] = predCon (selectedNoise.35[[i]]$random, selectedNoise.35[[i]]$sourceSemps)
  predSele.5[[i]] = predCon (selectedNoise.5[[i]]$random, selectedNoise.5[[i]]$sourceSemps)
}

qualSele.2 = vector(mode = 'list', length = length(predSele.2))
qualSele.5 = vector(mode = 'list', length = length(predSele.2))
qualSele.35 = vector(mode = 'list', length = length(predSele.2))

for (i in 1:length(predSele.2)){
  qualSele.2[[i]] = checkPred(selectedNoise.2[[i]]$contMat, predSele.2[[i]])
  qualSele.5[[i]] = checkPred(selectedNoise.5[[i]]$contMat, predSele.5[[i]])
  qualSele.35[[i]] = checkPred(selectedNoise.35[[i]]$contMat, predSele.35[[i]])
  
}

fr = data.frame(Pearson = c(unlist(qualSele.2), unlist(qualSele.35), unlist(qualSele.5)),
                RandomNoiseRatio = as.factor(repIndiv(c(0.2, 0.35, 0.5), length(unlist(qualSele.2)))),
                SampleRemoved = as.factor(rep(repIndiv(1:length(qualSele.2), length(qualSele.2[[1]])), 3)))

(ggplot(fr, aes(x = RandomNoiseRatio, y = Pearson))+geom_boxplot()+ ggtitle("Selected Samples as Noise") + theme(text = element_text(size=20)))


#One sample and noise -------------------
oneSampleNoise.2 = vector(mode = 'list', length= ncol(selectGroupData))
oneSampleNoise.35 = vector(mode = 'list', length= ncol(selectGroupData))
oneSampleNoise.5 = vector(mode = 'list', length= ncol(selectGroupData))
oneSampleNoise.6 = vector(mode = 'list', length= ncol(selectGroupData))

oneSampleNoise.8 = vector(mode = 'list', length= ncol(selectGroupData))

for (i in 1:ncol(selectGroupData)){
  oneSampleNoise.2[[i]] = generateSamples(.2, mixSamp = i)
  oneSampleNoise.35[[i]] = generateSamples(.35, mixSamp = i)
  oneSampleNoise.5[[i]] = generateSamples(.5, mixSamp = i)
  oneSampleNoise.6[[i]] = generateSamples(.6, mixSamp = i)
  
  oneSampleNoise.8[[i]] = generateSamples(.8, mixSamp = i)
  
}




predOneSampleNoise.2 = vector (mode = 'list', length = length(oneSampleNoise.2))
predOneSampleNoise.35 = vector (mode = 'list', length = length(oneSampleNoise.2))
predOneSampleNoise.5 = vector (mode = 'list', length = length(oneSampleNoise.2))
predOneSampleNoise.6 = vector (mode = 'list', length = length(oneSampleNoise.2))

predOneSampleNoise.8 = vector (mode = 'list', length = length(oneSampleNoise.2))

for (i in 1:length(oneSampleNoise.2)){
  predOneSampleNoise.2[[i]] = predCon (oneSampleNoise.2[[i]]$random, selectGroupData)
  predOneSampleNoise.35[[i]] = predCon (oneSampleNoise.35[[i]]$random,selectGroupData)
  predOneSampleNoise.5[[i]] = predCon (oneSampleNoise.5[[i]]$random, selectGroupData)
  predOneSampleNoise.6[[i]] = predCon (oneSampleNoise.6[[i]]$random, selectGroupData)
  predOneSampleNoise.8[[i]] = predCon (oneSampleNoise.8[[i]]$random, selectGroupData)
  
}



qualOneSampleNoise.2 = vector(mode = 'list', length = length(predSele.2))
qualOneSampleNoise.35 = vector(mode = 'list', length = length(predSele.2))
qualOneSampleNoise.5 = vector(mode = 'list', length = length(predSele.2))
qualOneSampleNoise.6 = vector(mode = 'list', length = length(predSele.2))
qualOneSampleNoise.8 = vector(mode = 'list', length = length(predSele.2))

for (i in 1:length(predOneSampleNoise.2)){
  qualOneSampleNoise.2[[i]] = checkPred(oneSampleNoise.2[[i]]$contMat, predOneSampleNoise.2[[i]])
  qualOneSampleNoise.35[[i]] = checkPred(oneSampleNoise.35[[i]]$contMat, predOneSampleNoise.35[[i]])
  qualOneSampleNoise.5[[i]] = checkPred(oneSampleNoise.5[[i]]$contMat, predOneSampleNoise.5[[i]])
  qualOneSampleNoise.6[[i]] = checkPred(oneSampleNoise.6[[i]]$contMat, predOneSampleNoise.6[[i]])
  qualOneSampleNoise.8[[i]] = checkPred(oneSampleNoise.8[[i]]$contMat, predOneSampleNoise.8[[i]])
  
}


fr = data.frame(Pearson =c(unlist(qualOneSampleNoise.2),unlist(qualOneSampleNoise.35),unlist(qualOneSampleNoise.5),unlist(qualOneSampleNoise.6),unlist(qualOneSampleNoise.8)),
                RandomNoiseRatio = as.factor(repIndiv(c(0.2, 0.35, 0.5, 0.6, 0.8), length(unlist(qualOneSampleNoise.2)))),
                WhichSample = as.factor(rep(repIndiv(1:length(qualOneSampleNoise.2), length(qualOneSampleNoise.2[[1]])), 5)))

(ggplot(fr, aes(x = RandomNoiseRatio, y = Pearson))+geom_boxplot()+ ggtitle("Selected single samples + random noise") + theme(text = element_text(size=20)))


# Noisy GO terms -----------------
#adds additioLnal noise to certain genes in a certain direction


mart = useMart("ensembl", dataset="mmusculus_gene_ensembl")

termsDeGo = getBM(attributes = c("mgi_symbol", "goslim_goa_accession",
                     "goslim_goa_description"), filters = "mgi_symbol",
      values = as.character(selectGeneData$Gene.Symbol), mart = mart)


geneGroupGo = c('ATPase activity', 'protein folding','cell death', 'protein complex assembly', 'aging', 'developmental maturation', "cell cycle")

geneNoise.0_3 = vector(mode = 'list', length = length(geneGroupGo))
geneNoise.0_5 = vector(mode = 'list', length = length(geneGroupGo))
geneNoise.2_3 = vector(mode = 'list', length = length(geneGroupGo))
geneNoise.35_3 = vector(mode = 'list', length = length(geneGroupGo))
geneNoise.5_3 = vector(mode = 'list', length = length(geneGroupGo))

for (i in 1:length(geneGroupGo)){
  genesToNoise = which(selectGeneData$Gene.Symbol %in% termsDeGo[termsDeGo$goslim_goa_description %in%  geneGroupGo[i], ]$mgi_symbol)
  geneNoise.0_3[[i]] = generateSamples(0, noisyGenes = genesToNoise, geneNoiseMean = 3, geneNoiseSD = 1)
  geneNoise.0_5[[i]] = generateSamples(0, noisyGenes = genesToNoise, geneNoiseMean = 5, geneNoiseSD = 1)
  geneNoise.2_3[[i]] = generateSamples(.2, noisyGenes = genesToNoise, geneNoiseMean = 3, geneNoiseSD = 1)
  geneNoise.35_3[[i]] = generateSamples(.35, noisyGenes = genesToNoise, geneNoiseMean = 3, geneNoiseSD = 1)
  geneNoise.5_3[[i]] = generateSamples(.5, noisyGenes = genesToNoise, geneNoiseMean = 3, geneNoiseSD = 1)
}



predGeneNoise.0_3 = vector (mode = 'list', length = length(geneNoise.0_3))
predGeneNoise.0_5 = vector (mode = 'list', length = length(geneNoise.0_5))
predGeneNoise.2_3 = vector (mode = 'list', length = length(geneNoise.0_3))
predGeneNoise.35_3 = vector (mode = 'list', length = length(geneNoise.0_3))
predGeneNoise.5_3 = vector (mode = 'list', length = length(geneNoise.0_3))

for (i in 1:length(geneNoise.0_3)){
  predGeneNoise.0_3[[i]] =  predCon (geneNoise.0_3[[i]]$random, selectGroupData)
  predGeneNoise.0_5[[i]] =  predCon (geneNoise.0_5[[i]]$random, selectGroupData)
  predGeneNoise.2_3[[i]] =  predCon (geneNoise.2_3[[i]]$random, selectGroupData)
  predGeneNoise.35_3[[i]] =  predCon (geneNoise.35_3[[i]]$random, selectGroupData)
  predGeneNoise.5_3[[i]]=  predCon (geneNoise.5_3[[i]]$random, selectGroupData)
}



qualGene.0_3 = vector(mode = 'list', length = length(geneNoise.0_3))
qualGene.0_5 = vector(mode = 'list', length = length(geneNoise.0_3))
qualGene.2_3 = vector(mode = 'list', length = length(geneNoise.0_3))
qualGene.35_3 = vector(mode = 'list', length = length(geneNoise.0_3))
qualGene.5_3 = vector(mode = 'list', length = length(geneNoise.0_3))
for (i in 1:length(geneNoise.0_3)){
  qualGene.0_3[[i]] = checkPred(geneNoise.0_3[[i]]$contMat, predGeneNoise.0_3[[i]])
  qualGene.0_5[[i]] = checkPred(geneNoise.0_5[[i]]$contMat, predGeneNoise.0_5[[i]])
  
  qualGene.2_3[[i]] = checkPred(geneNoise.2_3[[i]]$contMat, predGeneNoise.2_3[[i]])
  qualGene.35_3[[i]] = checkPred(geneNoise.35_3[[i]]$contMat, predGeneNoise.35_3[[i]])
  
  qualGene.5_3[[i]] = checkPred(geneNoise.5_3[[i]]$contMat, predGeneNoise.5_3[[i]])
}

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
fr = data.frame(Pearson = c(unlist(qualGene.0_3), unlist(qualGene.2_3), unlist(qualGene.35_3), unlist(qualGene.5_3)),
                RandomNoiseRatio = as.factor(repIndiv(c( 0, 0.2, 0.35, 0.5), length(unlist(qualGene.0_3)))),
                GeneAddMean = as.factor(c(rep(3, length(unlist(qualGene.0_3))) , rep(3, 3*length(unlist(qualGene.0_3))) )),
                GOTerms = rep(repIndiv(geneGroupGo, length(qualGene.0_3[[1]])), 4)
)
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
(ggplot(fr, aes(x = RandomNoiseRatio, y = Pearson))+geom_boxplot(aes())+geom_jitter(aes(color = GOTerms), size = 3, alpha = 0.8)+ ggtitle("Selected Gene Groups Noise and regular Noise (std = 1, mean = mostly 3)") + scale_colour_manual(values=cbbPalette)
+ theme_bw() + theme(text = element_text(size=20)))



# Generate completely shuffled nonsense samples

#shuffleAll =  selectGroupData [sample(nrow(selectGroupData), nrow(selectGroupData)),]





#######################################







#close(matlab)




source('/C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')
require(ggplot2)
require(reshape)
generateSamples = function(noiseMax = 0, randomNo = 30, maxMix = 5,mode = 'constant', noisyGenes = NA, geneNoiseMean = 0, geneNoiseSD = 1, mixSamp = NA, props = NA){
    #noiseMax = 0.01
    #randomNo = 20
    #maxMix = ncol(exprData)
    selectGroupData = 2^selectGroupData
    random =  matrix(0, nrow = nrow(selectGroupData), ncol = randomNo)
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
        if (any(is.na(mixSamp))){
            whichSamples = sample(ncol(selectGroupData), noOfSamples )
            #adding noise based on one of the samples but with shuffled titles
        } else {
            whichSamples = mixSamp
            sampleWeights = runif(length(mixSamp),0,1)
            noOfSamples = (length(mixSamp))
        }
        
        if (!(any(is.na(is.na(props))))){
            sampleWeights = props
        }
        
        addNoise = (selectGroupData[sample(nrow(selectGroupData), nrow(selectGroupData)), sample(ncol(selectGroupData), 1)]) * noise[i]
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
    outlist = list (random = log(random)/log(2), noise = noise, contMat = contMatrix)
    return(outlist)
}

daFolder = list.files("/C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween(1 July)(A moment later)",include.dirs = T)



# As always open the bloody files. ourNaming field added to the design file ----------------
design = read.table('/C:/Users/Ogan/Dropbox/Rotation 3/Data/normalizedDesign.csv',header=T,sep='\t')
design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)
allDataPre = read.csv('/C:/Users/Ogan/Dropbox/Rotation 3/Data/mostVariableQuantileNormalized', header = T)
geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData = exprData
design = design[match(colnames(exprData),design$sampleName,),]
rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<5)

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]
discludeSamples = which(is.na(design$ourNaming))
exprData = exprData[, -discludeSamples]
design = design[-discludeSamples, ]

guide = c('cellKind', 'ourNamingGeneral', 'ourNaming', 'justGaba', 'justGabaPV', 'someNaming')
step = 1
for (folder in daFolder){
    
    
    #fold genes ----
    filenames = list.files("/C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween(2 July)/"+folder,include.dirs = FALSE)
    fileContents = lapply("/C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween(2 July)/" + folder+'/'+ filenames, read.table)
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
    
    # get group means ----
    selectGroupData = matrix(0, nrow = nrow(exprData), ncol = length(trimNAs(unique(design[,guide[step]]))))
    groupNames = names(geneList)
    
    colnames(selectGroupData) = groupNames
    for (i in 1:length(groupNames)){
        groupIndex = design[,guide[step]] == groupNames[i]
        groupAverage =  tryCatch({apply(exprData[, groupIndex ], 1, median)},
                                 error = function(cond){
                                     print('fuu')
                                     return(exprData[, groupIndex ])
                                 })
        selectGroupData[, i] = groupAverage 
    }  
    
    
    
    
    for ( i in 1:ncol(selectGroupData)){
        corCheck = which(geneData$Gene.Symbol %in% foldGenes[[which(names(foldGenes) == unique(groups)[i])]])
        temp = matrix(nrow = length(corCheck), ncol = length(corCheck))
        
        
        sampNames = repIndiv(seq(0.9,0.1,-0.1),(sum(lower.tri(temp))) )
        
        corrs = vector( length = length(seq(0.9,0.1, -0.1)) * (sum(lower.tri(temp))))
        for (j in seq(0.9,0.1,-0.1)){
            m = j
            samplez = matrix(0, ncol = 30, nrow = nrow(selectGroupData))
            for ( z in 1:30){
                m = rnorm(1, j, sd = 0.05)
                if (m > 1){
                    m = 1
                }
                if (m < 0){
                    m = 0
                }
                randomSemp = generateSamples(noiseMax = 0.1, randomNo = 1, mixSamp = c(i, sample({1:ncol(selectGroupData)}[-i],4)) , props = c(m, rep((1-m)/4,4)))
                tryCatch({samplez[,z] = randomSemp$random}, warning = function(w){print(m)}) 
            }
            corr = cor(t(samplez[corCheck, ]))
            write.table( apply(samplez[corCheck, ],1 ,var), file = colnames(selectGroupData)[i]+'var' + j,  row.names = F, col.names = F)
            write.table(samplez[corCheck, ], file = colnames(selectGroupData)[i] + j,  row.names = F, col.names = F)
            corrs = corrs %c% unlist(corr[lower.tri(corr)])
        }
        
        fr = data.frame (corrs = abs(corrs), proportions = sampNames)
        write.table(fr, file = colnames(selectGroupData)[i], row.names = F, col.names = F)
        
        leFit = lm(corrs ~ proportions, data = fr)
        #median(fr[fr[,2] %in% 0.1,1])    
        cof = coef(leFit)
        
        prediction = predict(leFit,data.frame(proportions = seq(0.9,0.1, -0.1)) , interval = 'confidence')
        prediction = cbind(prediction, seq(0.9,0.1, -0.1))
        predFr = as.data.frame( prediction)
        
        predFrame = data.frame(lwr = repIndiv(prediction[,2], (sum(lower.tri(temp)))), upr  = repIndiv(prediction[,3], (sum(lower.tri(temp)))))
        fr = cbind(fr, predFrame)
        
        (ggplot(fr , aes_string(x= colnames(fr)[2], y= colnames(fr)[1]))
         + geom_violin(alpha=0.3, aes(group = fr[,2]))
         #+ geom_boxplot(aes(group = fr[,2]))
         + stat_summary(fun.y=median.quartile,color = 'red', geom='point')
         + ggtitle(colnames(selectGroupData)[i])
         + geom_abline(intercept =cof[1], slope = cof[2])
         + geom_ribbon(data = predFr,alpha = 0.4,aes(x =V4 , ymax= upr, ymin= lwr, y = NULL))
        )
        
        ggsave(filename =folder+'/coexpGen/' + colnames(selectGroupData)[i]+' ' + unique(groups)[j] + '.png')
        
    }
    
    
}

loc=getwd()
setwd('..')
parent = getwd()
setwd(loc)

source(paste0(parent,'/ogbox.r'))
require(ggplot2)
library(reshape2)

#load and arrange sample expression data ------
mouseRegionData = read.csv(paste0(parent,'/mouseRegions/mouseRegions'))
load(paste0(parent,'/Data/GSE12649_GSE5388_expression.RData'))
bpCntScz = aned_high_GSE12649
bpCnt = aned_high_GSE5388
bpCntSczDes = Samples_GSE12649
bpCntDes = Samples_GSE5388
rm(aned_high_GSE12649)
rm(aned_high_GSE5388)
rm(Samples_GSE12649)
rm(Samples_GSE5388)

# load files ----------
design = read.table(paste0(parent,'/Data/normalizedDesign.csv'),header=T,sep='\t')
design$anatomical = gsub('Corpus Striatum', 'Striatum', design$anatomical)
design$anatomical = gsub('Motor Cortex', 'Motor cortex', design$anatomical)
design$anatomical = gsub('Layer 5A Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 5B Cortex', 'Cortex', design$anatomical)
design$anatomical = gsub('Layer 6 Cortex', 'Cortex', design$anatomical)

design$cellType[grepl('^Neurons', design$cellType)] = 'Mixed Neurons'
design$cellType=gsub(', P.*?$','',design$cellType)
allDataPre = read.csv(paste0(parent,'/Data/mostVariableQuantileNormalized'), header = T)
geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]
exprData = exprData
design = design[match(colnames(exprData),design$sampleName,),]
rowmax = apply(exprData, 1, max)
discludeGenes = which(rowmax<5)

exprData = exprData[-discludeGenes,]
geneData = geneData[-discludeGenes,]



# dealing with gene selection --------
filenames = list.files(paste0(parent,"/AllGeneralComparison/Rest/InBetween/someNaming"),include.dirs = FALSE)
fileContents = lapply(paste0(parent,'/AllGeneralComparison/Rest/InBetween/someNaming/', filenames), read.table)
geneList = vector(mode = 'list', length = length(fileContents))
makePos = vector(mode = 'list', length = length(fileContents))
names(geneList) = filenames
for (i in 1:length(fileContents)){
    geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>log(10,base=2)])
    makePos[[i]] =  as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>log(10,base=2)])
}


# microglia exception ----
micro = read.table('C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/MicrogliaFromPaper')
geneList$Microglia = as.character(micro$V1[as.numeric(as.character(micro$V2))>log(10,base=2)])

puristList = vector(mode = 'list', length = length(geneList))
puristPos = vector(mode = 'list', length = length(geneList))
for (i in 1:length(geneList)){
    puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
    puristPos[[i]] = trimElement(makePos[[i]], unlist(makePos[-i]))
}


makePos = puristPos
names(puristList) = names(geneList)
commonGround = puristList
#commonGround= geneList
# 
# filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/silhouetteSelection/results/Groups 18 June",include.dirs = FALSE)
# fileContents = lapply(paste('C:/Users/Ogan/Dropbox/Rotation 3/silhouetteSelection/results/Groups 18 June/', filenames, sep = ''), read.table)
# 
# geneList = vector(mode = 'list', length = length(fileContents))
# names(geneList) = filenames
# options(warn=-1)
# for (i in 1:length(fileContents)){
#     geneList[[i]] = trimNAs(as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>0.5]))
# }
# options(warn=1)
# siloGenes = geneList
# 
# 
# commonGround = vector(mode = 'list', length = length(foldGenes))
# for (i in 1:length(foldGenes)){
#     commonGround[[i]] = intersect(siloGenes[[i]], foldGenes[[i]])
# }
# 
# names(commonGround) = names(siloGenes)
# 
# commonGround = foldGenes

# humanChip Adaptation --------
orthoInfo = read.csv(paste0(parent,'/Data/HG-U133A.na33.ortholog.csv'))
orthoInfo = orthoInfo[orthoInfo$Ortholog.Array == 'MOE430A',1:5 ]
humanGroundScz = vector(mode = 'list', length = length(commonGround))
humanGround = vector(mode = 'list', length = length(commonGround))
posGroundScz = vector(mode = 'list', length = length(commonGround))
posGround =  vector(mode = 'list', length = length(commonGround))
names(humanGround) = names(commonGround)
names(humanGroundScz) = names(commonGround)
for (i in 1:length(commonGround)){
    probeNames = geneData$Probe[geneData$Gene.Symbol %in% commonGround[[i]]]
    posProbes =  geneData$Probe[geneData$Gene.Symbol %in% makePos[[i]]]
    humanGroundScz[[i]] = as.character(bpCntScz[bpCntScz$Probe %in% orthoInfo$Probe.Set.ID[tolower(orthoInfo$Ortholog.Probe.Set) %in% probeNames], 'Gene Symbol'])
    humanGround[[i]] = as.character(bpCnt[bpCnt$Probe %in% orthoInfo$Probe.Set.ID[tolower(orthoInfo$Ortholog.Probe.Set) %in% probeNames], 'Gene Symbol'])
    
    posGroundScz[[i]] =  as.character(bpCntScz[bpCntScz$Probe %in% orthoInfo$Probe.Set.ID[tolower(orthoInfo$Ortholog.Probe.Set) %in% posProbes], 'Gene Symbol'])
    posGround[[i]] = as.character(bpCnt[bpCnt$Probe %in% orthoInfo$Probe.Set.ID[tolower(orthoInfo$Ortholog.Probe.Set) %in% posProbes], 'Gene Symbol'])
}



# marker change -----
mouseRegionExpr = mouseRegionData[,4:ncol(mouseRegionData)]
mouseRegionGenes = mouseRegionData[,1:3]


bpCntSczExpr = bpCntScz[,4:ncol(bpCntScz)]
bpCntSczGenes = bpCntScz[,1:3]


colnames(bpCntSczGenes) = colnames(geneData)

bpCntExpr = bpCnt[,4:ncol(bpCnt)]
bpCntGenes = bpCnt[,1:3]


colnames(bpCntGenes) = colnames(geneData)

windowSize = 4

rownames(bpCntSczExpr) = bpCntSczGenes$Gene.Symbol
rownames(bpCntExpr) = bpCntGenes$Gene.Symbol

usedStuff = vector(mode = 'list', length = length(commonGround))
names(usedStuff) = names(commonGround)
for (i in 1:length(commonGround)){

    trainer = exprData[geneData$Gene.Symbol %in% commonGround[[i]],]
    groups = design
    
    #bpCntSczR = t(scale(t(bpCntSczExpr[bpCntSczGenes$Gene.Symbol %in% humanGroundScz[[i]], ])))
    #bpCntR = t(scale(t(bpCntExpr[which(bpCntGenes$Gene.Symbol %in% humanGround[[i]]), ])))
    
    #bpCntSczR= apply(bpCntSczR,2, mean)
    #bpCntR = apply(bpCntR,2, mean)
    bpCntSczR = bpCntSczExpr[bpCntSczGenes$Gene.Symbol %in% humanGroundScz[[i]], ]
    
    #probing... get variance of genes
    varScz = apply(bpCntSczR,1,var)
    meanScz = apply(bpCntSczR,1,mean)
    
    bpCntR = bpCntExpr[which(bpCntGenes$Gene.Symbol %in% humanGround[[i]]), ]
    varCnt = apply(bpCntR, 1, var)
    meanCnt = apply(bpCntR,1,mean)
    
    pcaCntScz = prcomp(t(bpCntSczR), scale = T)
    
    pcaCntScz$rotation = pcaCntScz$rotation * ((sum(pcaCntScz$rotation[(rownames(pcaCntScz$rotation) %in% posGround[[i]]),1])<0)*(-2)+1)
    
    #report negative ones
    sczNeg = rownames(pcaCntScz$rotation)[pcaCntScz$rotation[,1]<(0)]
    sczVeryNeg = rownames(pcaCntScz$rotation)[pcaCntScz$rotation[,1]<(-0.1)]
    print(i)
    print(pcaCntScz$rotation[pcaCntScz$rotation[,1]<0,1])
    print('####')
    
    pcaCntScz$x = t(as.matrix(t(scale(t(bpCntSczR))))) %*% as.matrix(pcaCntScz$rotation)
    
    pcaCnt = prcomp(t(bpCntR), scale = T)
    
    pcaCnt$rotation = pcaCnt$rotation * ((sum(pcaCnt$rotation[(rownames(pcaCnt$rotation) %in% posGround[[i]]),1])<0)*(-2)+1)
    
    #report negative ones
    sczCnt = rownames(pcaCnt$rotation)[pcaCnt$rotation[,1]<(0)]
    print(pcaCnt$rotation[pcaCnt$rotation[,1]<0,1])
    
    #report common negatives
    print(intersect(sczNeg,sczCnt))
    
    
    pcaCnt$x = t(as.matrix(t(scale(t(bpCntR))))) %*% as.matrix(pcaCnt$rotation)
    
    
    bpCntSczR = t(pcaCntScz$x[,1])
    bpCntR = t(pcaCnt$x[,1])
    
    #bpCntSczR = bpCntSczExpr[bpCntSczGenes$Gene.Symbol %in% humanGroundScz[[i]], ]
    #bpCntR = bpCntExpr[bpCntGenes$Gene.Symbol %in% humanGround[[i]], ]
    
    #cntSczFrame = data.frame(zScore = unlist(as.data.frame(bpCntSczR)), sample = repIndiv(bpCntSczDes['disease_state',], nrow(bpCntSczR)))
    cntSczFrame = data.frame(PC1 = unlist(as.data.frame(bpCntSczR)), sample = bpCntSczDes['disease_state',])
    
    
    
    #collect info about the markers
    toMerge = list(varScz = (as.numeric(varScz)), 
                   varCnt = (as.numeric(varCnt)),
                   sczRot = (as.numeric(pcaCntScz$rotation[,1])),
                   CntRot = (as.numeric(pcaCnt$rotation[,1])),
                   meanScz = (as.numeric(meanScz)),
                   meanCnt = (as.numeric(meanCnt)))
        
        toMerge1 = list(varScz, 
                       varCnt)
    names(toMerge1[[1]]) = names(varScz)
    names(toMerge1[[2]]) = names(varCnt)
    
    toMerge2 = list(pcaCntScz$rotation[,1],
                    pcaCnt$rotation[,1])
    names(toMerge2[[1]]) = rownames(pcaCntScz$rotation)
    names(toMerge2[[2]]) = rownames(pcaCnt$rotation)
    
    toMerge3 = list(meanScz,
                    meanCnt)
    
    names(toMerge3[[1]]) = names(meanScz)
    names(toMerge3[[2]]) = names(meanCnt)
    
    usedStuff[[i]] = cbind(do.call(merge, c(toMerge1,by=0, all=TRUE)), 
                           do.call(merge, c(toMerge2,by=0, all=TRUE)),
                           do.call(merge, c(toMerge3,by=0, all=TRUE)))
    
    usedStuff[[i]]= usedStuff[[i]][,-c(4,7)]
    colnames(usedStuff[[i]]) = c('gene','varScz', 'varCnt', 'sczRot', 'cntRot', 'meanScz', 'meanCnt')
    
    contBP=  tryCatch({
        wilcox.test(cntSczFrame[cntSczFrame$sample=='BP',1], cntSczFrame[cntSczFrame$sample=='Cont',1])
    },error = function(e){
        return('Er')
    })
    contScz = tryCatch({
     wilcox.test(cntSczFrame[cntSczFrame$sample=='SCZ',1], cntSczFrame[cntSczFrame$sample=='Cont',1])
    },error = function(e){
       return('Er')
    })   
    tryCatch({
        (p = ggplot(cntSczFrame, aes(x =sample, y =PC1  ))
         + geom_violin( color="#C4C4C4", fill="#C4C4C4")
         + geom_boxplot(width=0.1,fill = 'lightblue')
         + ggtitle('GSE12649 ' + names(commonGround)[i])
         + annotate('text' , x = 1.5, y =1, label = 'p = ' +  round(contBP$p.value,digits = 5), size = 8)
         + annotate('text' , x = 2.5, y =1, label = 'p = ' +  round(contScz$p.value,digits = 5), size = 8)
         + scale_y_continuous(limits=c(-windowSize, windowSize))
         + theme_bw()
         + theme(axis.text.x  = element_text(vjust=windowSize/2, size=20), axis.title.y = element_text(vjust=0.5, size=20),axis.title.x = element_text(vjust=0.5, size=0) , title = element_text(vjust=0.5, size=20))
        )

    ggsave(filename = 'GSE12649 ' + names(commonGround)[i] + '.png')
    }, error = function(e){
        
    })
    
    #cntFrame = data.frame(zScore = unlist(as.data.frame(bpCntR)), sample = gsub('_t','',repIndiv(bpCntDes['disease_state',], nrow(bpCntR))))
    cntFrame = data.frame(PC1 = unlist(as.data.frame(bpCntR)), sample = gsub('_t','',bpCntDes['disease_state',]))
    
    
    contBP = tryCatch({
    wilcox.test(cntFrame[cntFrame$sample=='BP',1], cntFrame[cntFrame$sample=='Cont',1])
    },error = function(e){
        return( 'Er')
    })
    tryCatch({
    (p = ggplot(cntFrame, aes(x =sample, y =PC1  ))
     + geom_violin( color="#C4C4C4", fill="#C4C4C4")
     + ggtitle('GSE5388 ' + names(commonGround)[i])
     + geom_boxplot(width=0.1,fill = 'lightblue')
     + annotate('text' , x = 1.5, y =1, label = 'p = ' +  round(contBP$p.value,digits = 5), size = 8)
     + scale_y_continuous(limits=c(-windowSize, windowSize))
     + theme_bw()
     + theme(axis.text.x  = element_text(vjust=windowSize/2, size=20), axis.title.y = element_text(vjust=0.5, size=20),axis.title.x = element_text(vjust=0.5, size=0) , title = element_text(vjust=0.5, size=20))
    )
    
    ggsave(filename = 'GSE5388 ' + names(commonGround)[i] + '.png')
    }, error = function(e){
    
    })
    
    
    mouseR = t(scale(t(mouseRegionExpr[mouseRegionGenes$Probe %in% (geneData$Probe[geneData$Gene.Symbol %in%  commonGround[[i]] ]), ])))
    mouseFrame = data.frame(PC1 = unlist(as.data.frame(mouseR)), sample = repIndiv(c(rep('cereb',3),rep('cort',3)), nrow(mouseR)))
    cerebCort = tryCatch({
      wilcox.test(mouseFrame[mouseFrame$sample=='cereb',1], mouseFrame[mouseFrame$sample=='cort',1])
    },error = function(e){
        return('Er')
    })
    tryCatch({
    (p = ggplot(mouseFrame, aes(x =sample, y =PC1  ))
     + geom_violin( color="#C4C4C4", fill="#C4C4C4")
     + geom_boxplot(width=0.1,fill = 'lightblue')
     + ggtitle('mouseRegion ' + names(commonGround)[i])
     + annotate('text' , x = 1.5, y =1, label = 'p = ' +  round(cerebCort$p.value,digits = 5), size = 8)
     + scale_y_continuous(limits=c(-windowSize, windowSize))
     + theme(axis.text.x  = element_text(vjust=windowSize/2, size=20), axis.title.y = element_text(vjust=0.5, size=20),axis.title.x = element_text(vjust=0.5, size=0) , title = element_text(vjust=0.5, size=20))
     + theme_bw()
    )
    ggsave(filename = 'mouseRegion ' + names(commonGround)[i] + '.png')
    }, error = function(e){
        
    })
 
    
}


i=6
usedStuff[[i]][order(usedStuff[[i]]$sczRot,usedStuff[[i]]$cntRot ), ]




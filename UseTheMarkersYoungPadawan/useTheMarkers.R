source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')
require(ggplot2)
library(reshape2)

#load and arrange sample expression data ------
mouseRegionData = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/mouseRegions/mouseRegions')
load('C:/Users/Ogan/Dropbox/Rotation 3/Data/GSE12649_GSE5388_expression.RData')
bpCntScz = aned_high_GSE12649
bpCnt = aned_high_GSE5388
bpCntSczDes = Samples_GSE12649
bpCntDes = Samples_GSE5388
rm(aned_high_GSE12649)
rm(aned_high_GSE5388)
rm(Samples_GSE12649)
rm(Samples_GSE5388)

# load files ----------
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



# dealing with gene selection --------
filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/justGabaPV",include.dirs = FALSE)
fileContents = lapply(paste('C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/justGabaPV/', filenames, sep = ''), read.table)
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
commonGround = puristList
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
orthoInfo = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/Data/HG-U133A.na33.ortholog.csv')
orthoInfo = orthoInfo[orthoInfo$Ortholog.Array == 'MOE430A',1:5 ]
humanGroundScz = vector(mode = 'list', length = length(commonGround))
humanGround = vector(mode = 'list', length = length(commonGround))
names(humanGround) = names(commonGround)
names(humanGroundScz) = names(commonGround)
for (i in 1:length(commonGround)){
    probeNames = geneData$Probe[geneData$Gene.Symbol %in% commonGround[[i]]]
    humanGroundScz[[i]] = as.character(bpCntScz[bpCntScz$Probe %in% orthoInfo$Probe.Set.ID[tolower(orthoInfo$Ortholog.Probe.Set) %in% probeNames], 'Gene Symbol'])
    humanGround[[i]] = as.character(bpCnt[bpCnt$Probe %in% orthoInfo$Probe.Set.ID[tolower(orthoInfo$Ortholog.Probe.Set) %in% probeNames], 'Gene Symbol'])

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


for (i in 1:length(commonGround)){

    trainer = exprData[geneData$Gene.Symbol %in% commonGround[[i]],]
    groups = design
    
    bpCntSczR = t(scale(t(bpCntSczExpr[bpCntSczGenes$Gene.Symbol %in% humanGroundScz[[i]], ])))
    bpCntR = t(scale(t(bpCntExpr[bpCntGenes$Gene.Symbol %in% humanGround[[i]], ])))
    #bpCntSczR = bpCntSczExpr[bpCntSczGenes$Gene.Symbol %in% humanGroundScz[[i]], ]
    #bpCntR = bpCntExpr[bpCntGenes$Gene.Symbol %in% humanGround[[i]], ]
    
    cntSczFrame = data.frame(zScore = unlist(as.data.frame(bpCntSczR)), sample = repIndiv(bpCntSczDes['disease_state',], nrow(bpCntSczR)))
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
        (p = ggplot(cntSczFrame, aes(x =sample, y =zScore  ))
         + geom_violin(alpha=0.3)
         #+ geom_boxplot()
         + ggtitle('GSE12649 ' + names(commonGround)[i])
         + stat_summary(fun.y=median.quartile,shape='_',size=10,color = 'red', geom='point')
         + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
         + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
         + annotate('text' , x = 1.5, y =2, label = 'p = ' +  round(contBP$p.value,digits = 5))
         + annotate('text' , x = 2.5, y =2, label = 'p = ' +  round(contScz$p.value,digits = 5))
         + scale_y_continuous(limits=c(-4, 4))
         + theme(text = element_text(size=20))
        )

    ggsave(filename = 'GSE12649 ' + names(commonGround)[i] + '.png')
    }, error = function(e){
        
    })
    
    cntFrame = data.frame(zScore = unlist(as.data.frame(bpCntR)), sample = gsub('_t','',repIndiv(bpCntDes['disease_state',], nrow(bpCntR))))
    contBP = tryCatch({
    wilcox.test(cntFrame[cntFrame$sample=='BP',1], cntFrame[cntFrame$sample=='Cont',1])
    },error = function(e){
        return( 'Er')
    })
    tryCatch({
    (p = ggplot(cntFrame, aes(x =sample, y =zScore  ))
     + geom_violin(alpha=0.3)
     #+ geom_boxplot()
     + ggtitle('GSE5388 ' + names(commonGround)[i])
     + stat_summary(fun.y=median.quartile,shape='_',size=10,color = 'red', geom='point')
     + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
     + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
     + annotate('text' , x = 1.5, y =2, label = 'p = ' +  round(contBP$p.value,digits = 5))
     + theme(text = element_text(size=20))
     + scale_y_continuous(limits=c(-4, 4))
    )
    
    ggsave(filename = 'GSE5388 ' + names(commonGround)[i] + '.png')
    }, error = function(e){
    
    })
    
    
    mouseR = t(scale(t(mouseRegionExpr[mouseRegionGenes$Probe %in% (geneData$Probe[geneData$Gene.Symbol %in%  commonGround[[i]] ]), ])))
    mouseFrame = data.frame(zScore = unlist(as.data.frame(mouseR)), sample = repIndiv(c(rep('cereb',3),rep('cort',3)), nrow(mouseR)))
    cerebCort = tryCatch({
      wilcox.test(mouseFrame[mouseFrame$sample=='cereb',1], mouseFrame[mouseFrame$sample=='cort',1])
    },error = function(e){
        return('Er')
    })
    tryCatch({
    (p = ggplot(mouseFrame, aes(x =sample, y =zScore  ))
     + geom_violin(alpha=0.3)
     #+ geom_boxplot()
     + ggtitle('mouseRegion ' + names(commonGround)[i])
     + annotate('text' , x = 1.5, y =2, label = 'p = ' +  round(cerebCort$p.value,digits = 5))
     + stat_summary(fun.y=median.quartile,shape='_',size=10,color = 'red', geom='point')
     + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
     + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
     + theme(text = element_text(size=20))
     + scale_y_continuous(limits=c(-4, 4))
    )
    ggsave(filename = 'mouseRegion ' + names(commonGround)[i] + '.png')
    }, error = function(e){
        
    })
 
    
}
    
    
    
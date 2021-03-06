#standardize everything in lists. live happily ever after


source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')
insist(ggplot2)
insist(caret)
require(rgl)
insist(e1071)
require(reshape)
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
discludeSamples = which(is.na(design$ourNaming))
exprData = exprData[, -discludeSamples]
design = design[-discludeSamples, ]


filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/justGaba",include.dirs = FALSE)
fileContents = lapply(paste('C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/justGaba/', filenames, sep = ''), read.table)
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

exprData = exprData[geneData$Gene.Symbol %in% unlist(puristList),]
geneData = geneData[geneData$Gene.Symbol %in% unlist(puristList),]

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

# fix naming
colnames(bpCnt)[2] = 'Gene.Symbol'
colnames(bpCntScz)[2] = 'Gene.Symbol'



# initial PCA ===================
#expression values are stored here
testExpr = list (GSE12649 = bpCntScz[4:ncol(bpCntScz)],GSE5388 = bpCnt[4:ncol(bpCnt)], mouseRegion = mouseRegionData[4:ncol(mouseRegionData)] )

#gene data is stored here
testGene = list(GSE12649 = bpCntScz[1:3],GSE5388 = bpCnt[1:3], mouseRegion = mouseRegionData[1:3] )

#
orthoInfo = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/Data/HG-U133A.na33.ortholog.csv')
orthoInfo = orthoInfo[orthoInfo$Ortholog.Array == 'MOE430A',1:5 ]



exprDatas = vector(mode = 'list', length = length(testGene))
geneDatas = vector(mode = 'list', length = length(testGene))

groups = design$ourNamingGeneral
groups = design$justGaba
groups2 = design$ourNaming

testGroups = list (
    GSE12649 = bpCntSczDes['disease_state',],
    GSE5388 = gsub('_t','',bpCntDes['disease_state', ]),
    mouseRegion = c('cere','cere','cere', 'cort', 'cort', 'cort')
    )

for (i in 1:length(testExpr)){
    #find genes with matching orthologues in the humans. do not apply to mouse samples.
    if (!names(testExpr)[i] %in% 'mouseRegion'){
        usePairs = orthoInfo[{tolower(orthoInfo$Ortholog.Probe.Set) %in% tolower(geneData$Probe)}&
                             {tolower(orthoInfo$Probe.Set.ID) %in% tolower(testGene[[i]]$Probe)}
                             ,c(1,3)]
        usePairs = usePairs[!(duplicated(usePairs[,1]) | duplicated(usePairs[,2])),]
        
        exprDatas[[i]] = exprData[match(tolower(usePairs$Ortholog.Probe.Set), tolower(geneData$Probe)), ] 
        geneDatas[[i]] = geneData[match(tolower(usePairs$Ortholog.Probe.Set), tolower(geneData$Probe)), ]
        
        testExpr[[i]] = testExpr[[i]][match(usePairs$Probe.Set.ID, testGene[[i]]$Probe), ] 
        testGene[[i]] = testGene[[i]][match(usePairs$Probe.Set.ID, testGene[[i]]$Probe), ] 
    } else {
        exprDatas[[i]] = exprData[ - which( !geneData$Probe %in% testGene[[i]]$Probe ), ]
        geneDatas[[i]] = geneData[ - which( !geneData$Probe %in% testGene[[i]]$Probe ), ]
        
        testExpr[[i]] = testExpr[[i]][ - which( !testGene[[i]]$Probe %in% geneDatas[[i]]$Probe ), ]
        testGene[[i]] = testGene[[i]][ - which( !testGene[[i]]$Probe %in% geneDatas[[i]]$Probe ), ]
        
        testExpr[[i]] = testExpr[[i]][match(geneDatas[[i]]$Probe, testGene[[i]]$Probe), ] 
        testGene[[i]] = testGene[[i]][match(geneDatas[[i]]$Probe, testGene[[i]]$Probe), ] 
    }
    
    allTogether =scale(t(cbind(exprDatas[[i]], testExpr[[i]])))
    
    
    #use first 4 pcas to train the classifier
    pca = prcomp(allTogether[1:ncol(exprDatas[[i]]), ], scale = T)
    
    testPC = as.matrix(allTogether[(ncol(exprDatas[[i]])+1):nrow(allTogether) , ]) %*% as.matrix(pca$rotation)
    classifier = naiveBayes(pca$x[,1:3], groups)
    prediction = predict(classifier, testPC[,1:3], 'raw',threshold = 0.001)
    groupMerge = cbind(prediction, testGroups[[i]])
    
    
    
    
    for (j in 5:length(unique(groups))){
        print(j)
        
        if (length(unique(testGroups[[i]]))== 3){
            fr = as.data.frame(groupMerge[ , c(j,ncol(groupMerge))])
            fr[,1] = (as.numeric(as.character(fr[,1])))
            p1 = wilcox.test(fr[,1][fr$V2 == unique(fr$V2)[1]], fr[,1][fr$V2 == unique(fr$V2)[2] ])
            p2 = wilcox.test(fr[,1][fr$V2 == unique(fr$V2)[2]], fr[,1][fr$V2 == unique(fr$V2)[3] ])
            fr [,1] = log(fr[,1])
          
            (ggplot(fr , aes_string(x= colnames(fr)[2], y= colnames(fr)[1]))
             + geom_violin(alpha=0.3)
             + geom_boxplot()
             + ggtitle(names(testExpr)[i]+' ' + colnames(fr)[1])
             + annotate('text' , x = 1.5, y = max(fr[,1]), label = 'p = ' +  round(p1$p.value,digits = 5))
             + annotate('text' , x = 2.5, y = max(fr[,1]), label = 'p = ' +  round(p2$p.value,digits = 5))
            )
            ggsave(filename = names(testExpr)[i]+' ' + colnames(fr)[1] + '.png')
        } else if (length(unique(testGroups[[i]]))== 2){
            
            fr = as.data.frame(groupMerge[ , c(j,ncol(groupMerge))])
            fr[,1] = (as.numeric(as.character(fr[,1])))
            p1 = wilcox.test(fr[,1][fr$V2 == unique(fr$V2)[1]], fr[,1][fr$V2 == unique(fr$V2)[2] ])
            fr[,1] = log(fr[,1])
            
           
            (ggplot(fr , aes_string(x= colnames(fr)[2], y= colnames(fr)[1]))
             + geom_violin(alpha=0.3)
             + geom_boxplot()
             + ggtitle(names(testExpr)[i]+' ' + colnames(fr)[1])
             + annotate('text' , x = 1.5, y = max(fr[,1]), label = 'p = ' +  round(p1$p.value,digits = 5))
            )
            ggsave(filename = names(testExpr)[i]+' ' + colnames(fr)[1] + '.png')
        }
            
        
    }
    
}


colors = toColor(groups,rainbow(length(unique(groups))))

plot3d(pca$x[,1], pca$x[,2], pca$x[,3], col = colors, size = 9 )

usePairs = orthoInfo[{tolower(orthoInfo$Ortholog.Probe.Set) %in% tolower(geneData$Probe)}&
{tolower(orthoInfo$Probe.Set.ID) %in% tolower(testGene[[i]]$Probe)}
,c(1,3)]

usePairs[(duplicated(usePairs[,1]) | duplicated(usePairs[,2])),]



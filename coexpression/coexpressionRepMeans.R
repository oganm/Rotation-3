source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')

daFolder = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween",include.dirs = T)
require(ggplot2)
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

mouseRegionData = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/mouseRegions/mouseRegions')
mouseRegionData = cbind(mouseRegionData, log((2^mouseRegionData[,4:6]+2^mouseRegionData[,7:9])/2, base= 2))

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






testGroups = list (
    GSE12649 = bpCntSczDes['disease_state',],
    GSE5388 = gsub('_t','',bpCntDes['disease_state', ]),
    mouseRegion = c('cere','cere','cere', 'cort', 'cort', 'cort', 'mix', 'mix', 'mix')
)



for (folder in daFolder){
    exprDatas = vector(mode = 'list', length = length(testGene))
    geneDatas = vector(mode = 'list', length = length(testGene))
    testExpr = list (GSE12649 = bpCntScz[4:ncol(bpCntScz)],GSE5388 = bpCnt[4:ncol(bpCnt)], mouseRegion = mouseRegionData[4:ncol(mouseRegionData)] )
    #gene data is stored here
    testGene = list(GSE12649 = bpCntScz[1:3],GSE5388 = bpCnt[1:3], mouseRegion = mouseRegionData[1:3] )
    dir.create('sigs/')
    dir.create(folder + '/')
    dir.create('sigs/'+folder)
    
    filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/"+folder,include.dirs = FALSE)
    fileContents = lapply(paste('C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/'+folder+'/', filenames, sep = ''), read.table)
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
    
    for (i in 1:length(testGene)){
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
            # exprDatas[[i]] = exprData[ -which( !geneData$Probe %in% testGene[[i]]$Probe ), ]
            # geneDatas[[i]] = geneData[ -which( !geneData$Probe %in% testGene[[i]]$Probe ), ]
            exprDatas[[i]] = exprData
            geneDatas[[i]] = geneData
            
            testExpr[[i]] = testExpr[[i]][ - which( !testGene[[i]]$Probe %in% geneDatas[[i]]$Probe ), ]
            testGene[[i]] = testGene[[i]][ - which( !testGene[[i]]$Probe %in% geneDatas[[i]]$Probe ), ]
            
            testExpr[[i]] = testExpr[[i]][match(geneDatas[[i]]$Probe, testGene[[i]]$Probe), ] 
            testGene[[i]] = testGene[[i]][match(geneDatas[[i]]$Probe, testGene[[i]]$Probe), ] 
        }
        
        dir.create(folder + '/neg')
        dir.create('sigs/'+ folder + '/neg')
        
        dir.create(folder + '/pos')
        dir.create('sigs/' + folder + '/pos')
        
        for (j in 1:length(filenames)){
            corCheck = testExpr[[i]][geneDatas[[i]]$Gene.Symbol %in% foldGenes[[which(names(foldGenes)==filenames[j])]],]
            rownames(corCheck) = geneDatas[[i]]$Gene.Symbol[geneDatas[[i]]$Gene.Symbol %in% foldGenes[[which(names(foldGenes)==filenames[j])]]]
            corCheck = corCheck[apply(corCheck,1,mean)>8,]
            temp = matrix(nrow = nrow(corCheck), ncol = nrow(corCheck))
            groupNames = unique(testGroups[[i]])
            correlations = vector(length = length(groupNames)*(sum(lower.tri(temp))))
            if (length(correlations)==0){
                next
            }
            sampNames = repIndiv(groupNames, (sum(lower.tri(temp))))
    
            for (z in 1:length(groupNames)){
                corelasyon = cor(t(corCheck[, testGroups[[i]] %in% groupNames[z]]), method ='spearman')
                apply(corelasyon<0,2, any)
        
                correlations[{(sum(lower.tri(temp)))*z-(sum(lower.tri(temp)))+1}:{(sum(lower.tri(temp)))*z}] = corelasyon[lower.tri(corelasyon, diag = FALSE)]
            }
    
            daFrame = data.frame((correlations), sampNames)
            posDaFrame = daFrame[daFrame$X.correlations.>0,]
            negDaFrame = daFrame[daFrame$X.correlations.<0,]
            daFrame$X.correlations. = abs( daFrame$X.correlations.)
    
            if (length(unique(testGroups[[i]]))== 3){
                tryCatch({p1 = wilcox.test(daFrame[,1][daFrame[,2] == groupNames[1]], daFrame[,1][daFrame[,2] == groupNames[2] ])},error= function(e){p1.pvalue=1})
                tryCatch({p2 = wilcox.test(daFrame[,1][daFrame[,2] == groupNames[2]], daFrame[,1][daFrame[,2] == groupNames[3] ])},error= function(e){p2.pvalue=1})
                
                p1$p.value= (p1$p.value * length(filenames))
                p2$p.value= (p2$p.value * length(filenames))
                
                
                tryCatch({(ggplot(daFrame , aes_string(x= colnames(daFrame)[2], y= colnames(daFrame)[1]))
                + geom_violin(alpha=0.3)
                #+ geom_boxplot()
                + ggtitle(names(testExpr)[i]+' ' + filenames[j])
                + annotate('text' , x = 1.5, y = max(daFrame[,1]), label = 'p = ' +  round(p1$p.value,digits = 5))
                + annotate('text' , x = 2.5, y = max(daFrame[,1]), label = 'p = ' +  round(p2$p.value,digits = 5))
                + stat_summary(fun.y=median.quartile,shape='_',size=10,color = 'red', geom='point')
                + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
                + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
                + theme(text = element_text(size=20))
                )
                aster = ''
                if (p1$p.value<0.05 | p2$p.value<0.05){
                    aster ='(sig)'
                    ggsave(filename ='sigs/'+ folder+'/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                    
                }
                ggsave(filename =folder+'/'+ names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                print(names(testExpr)[i]+' ' + filenames[j] + '.png'+' 3')}, error = function(e){})
                ###### nedDaFrame
                tryCatch({p1 = wilcox.test(negDaFrame[,1][negDaFrame[,2] == groupNames[1]], negDaFrame[,1][negDaFrame[,2] == groupNames[2] ])},error= function(e){p1.pvalue=1})
                tryCatch({p2 = wilcox.test(negDaFrame[,1][negDaFrame[,2] == groupNames[2]], negDaFrame[,1][negDaFrame[,2] == groupNames[3] ])},error= function(e){p2.pvalue=1})
                
                p1$p.value= (p1$p.value * length(filenames))
                p2$p.value= (p2$p.value * length(filenames))
                tryCatch({(ggplot(negDaFrame , aes_string(x= colnames(negDaFrame)[2], y= colnames(negDaFrame)[1]))
                 + geom_violin(alpha=0.3)
                 #+ geom_boxplot()
                 + ggtitle(names(testExpr)[i]+' ' + filenames[j])
                 + annotate('text' , x = 1.5, y = max(negDaFrame[,1]), label = 'p = ' +  round(p1$p.value,digits = 5))
                 + annotate('text' , x = 2.5, y = max(negDaFrame[,1]), label = 'p = ' +  round(p2$p.value,digits = 5))
                 + stat_summary(fun.y=median.quartile,shape='_',size=10,color = 'red', geom='point')
                 + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
                 + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
                 + theme(text = element_text(size=20))
                 )
                aster = ''
                if (p1$p.value<0.05 | p2$p.value<0.05){
                    aster ='(sig)'
                    ggsave(filename ='sigs/'+ folder + '/neg/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                    
                }
                ggsave(filename =folder+'/neg/'+ names(testExpr)[i]+' ' + filenames[j] + aster+ '.png')
                print(names(testExpr)[i]+' ' + filenames[j] + '.png'+' 3')}, error = function(e){})
                
                #### posDaFrame
                tryCatch({p1 = wilcox.test(posDaFrame[,1][posDaFrame[,2] == groupNames[1]], posDaFrame[,1][posDaFrame[,2] == groupNames[2] ])},error= function(e){p1.pvalue=1})
                tryCatch({p2 = wilcox.test(posDaFrame[,1][posDaFrame[,2] == groupNames[2]], posDaFrame[,1][posDaFrame[,2] == groupNames[3] ])},error= function(e){p2.pvalue=1})
                
                p1$p.value= (p1$p.value * length(filenames))
                p2$p.value= (p2$p.value * length(filenames))
                tryCatch({(ggplot(posDaFrame , aes_string(x= colnames(posDaFrame)[2], y= colnames(posDaFrame)[1]))
                 + geom_violin(alpha=0.3)
                 #+ geom_boxplot()
                 + ggtitle(names(testExpr)[i]+' ' + filenames[j])
                 + annotate('text' , x = 1.5, y = max(posDaFrame[,1]), label = 'p = ' +  round(p1$p.value,digits = 5))
                 + annotate('text' , x = 2.5, y = max(posDaFrame[,1]), label = 'p = ' +  round(p2$p.value,digits = 5))
                 + stat_summary(fun.y=median.quartile,shape='_',size=10,color = 'red', geom='point')
                 + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
                 + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
                 +theme(text = element_text(size=20))
                )
                aster = ''
                if (p1$p.value<0.05 | p2$p.value<0.05){
                    aster ='(sig)'
                    ggsave(filename ='sigs/' + folder + '/pos/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                    
                }
                ggsave(filename =folder+'/pos/'+ names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                print(names(testExpr)[i]+' ' + filenames[j] + '.png'+' 3')}, error = function(e){})
                
                
                
            }
            if (length(unique(testGroups[[i]])) == 2){
                tryCatch({p1 = wilcox.test(daFrame[,1][daFrame[,2] == groupNames[1]], daFrame[,1][daFrame[,2] == groupNames[2] ])},error= function(e){p1.pvalue=1})
                
                p1$p.value= (p1$p.value * length(filenames))
                tryCatch({(ggplot(daFrame , aes_string(x= colnames(daFrame)[2], y= colnames(daFrame)[1]))
                + geom_violin(alpha=0.3)
                #+ geom_boxplot()
                + ggtitle(names(testExpr)[i]+' ' + filenames[j])
                + annotate('text' , x = 1.5, y = max(daFrame[,1]), label = 'p = ' +  round(p1$p.value,digits = 5))
                + stat_summary(fun.y=median.quartile,shape='_',size=10,color = 'red', geom='point')
                + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
                + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
                + theme(text = element_text(size=20))
                )
                aster = ''
                if (p1$p.value<0.05){
                    aster ='(sig)'
                    ggsave(filename ='sigs/'+ folder +'/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                    
                }
                ggsave(filename = folder + '/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                print(names(testExpr)[i]+' ' + filenames[j] + '.png' + ' 2')}, error = function(e){})
                
                ### negDaFrame 
                tryCatch({p1 = wilcox.test(negDaFrame[,1][daFrame[,2] == groupNames[1]], negDaFrame[,1][negDaFrame[,2] == groupNames[2] ])},error= function(e){p1.pvalue=1})
                p1$p.value= (p1$p.value * length(filenames))
                
                tryCatch({(ggplot(negDaFrame , aes_string(x= colnames(negDaFrame)[2], y= colnames(negDaFrame)[1]))
                 + geom_violin(alpha=0.3)
                 #+ geom_boxplot()
                 + ggtitle(names(testExpr)[i]+' ' + filenames[j])
                 + annotate('text' , x = 1.5, y = max(negDaFrame[,1]), label = 'p = ' +  round(p1$p.value,digits = 5))
                 + stat_summary(fun.y=median.quartile,shape='_',color = 'red', geom='point')
                 + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
                 + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
                 + theme(text = element_text(size=20))
                )
                aster = ''
                if (p1$p.value<0.05){
                    aster ='(sig)'
                    ggsave(filename ='sigs/'+ folder + '/neg/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                    
                }
                ggsave(filename = folder + '/neg/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                print(names(testExpr)[i]+' ' + filenames[j] + '.png' + ' 2')}, error = function(e){})
                
                #posDaFrame
                tryCatch({p1 = wilcox.test(posDaFrame[,1][posDaFrame[,2] == groupNames[1]], posDaFrame[,1][posDaFrame[,2] == groupNames[2] ])},error= function(e){p1.pvalue=1})
                p1$p.value= (p1$p.value * length(filenames))
                
                tryCatch({(ggplot(posDaFrame , aes_string(x= colnames(posDaFrame)[2], y= colnames(posDaFrame)[1]))
                 + geom_violin(alpha=0.3)
                 #+ geom_boxplot()
                 + ggtitle(names(testExpr)[i]+' ' + filenames[j])
                 + annotate('text' , x = 1.5, y = max(posDaFrame[,1]), label = 'p = ' +  round(p1$p.value,digits = 5))
                 + stat_summary(fun.y=median.quartile,shape='_',size=10,color = 'red', geom='point')
                 + stat_summary(fun.y=threeQuartile,shape='_',size=10,color = 'orange', geom='point')
                 + stat_summary(fun.y=median,shape = '_', size = 10, color = 'red',geom= 'point')
                 + theme(text = element_text(size=20))
                )
                aster = ''
                if (p1$p.value<0.05){
                    aster ='(sig)'
                    ggsave(filename = 'sigs/' + folder + '/pos/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                    
                }
                ggsave(filename = folder + '/pos/' + names(testExpr)[i]+' ' + filenames[j] + aster + '.png')
                print(names(testExpr)[i]+' ' + filenames[j] + '.png' + ' 2')}, error = function(e){})
            }
    
        }
    
    }

}
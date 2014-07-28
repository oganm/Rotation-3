source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')
require(cluster)
require(gplots)


daFolder = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween(2 July)",include.dirs = T)

source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')
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


daFolder = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween(2 July)",include.dirs = T)


siloMat = matrix(nrow =length(daFolder) , ncol = 2)
rownames(siloMat) = daFolder
stepFol = 1
for (folder in daFolder){
    filenames = list.files("C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween(2 July)/"+folder,include.dirs = FALSE)
    fileContents = lapply(paste('C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween(2 July)/'+folder+'/', filenames, sep = ''), read.table)
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
    stepi = 1    
    for (i in list(puristList, geneList)){
         clustering = vector(length = nrow(design))
         for (j in 1:length(filenames)){
             clustering[design[, folder]==filenames[j]] = j
         }
         clustering =  trimElement(clustering,0)
         data = t(exprData[geneData$Gene.Symbol %in% unlist(i) ,which(!is.na(design[,folder]))])
         
         cluster = list(clustering = clustering, data = data)    
         silo = silhouette(cluster,dist(data))
         siloMat[stepFol, stepi] = mean(silo[,3])
         
         
         png(folder +' '+c('purist','nonPurist')[stepi] + '.jpg', width =600, height = 600)
         colors = toColor(clustering,rainbow(length(unique(clustering))))
         heatmap.2(cor(t(data)), dendrogram = 'column', trace = 'none', ColSideColors= colors )
         legend("bottomleft", legend= filenames,
                fill = colors[match(1:max(clustering), clustering)],cex=0.8)
         dev.off()
         stepi = stepi + 1
    }
    
    stepFol = stepFol+1
}





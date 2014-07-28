source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')
source('specificity_index.R')

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

nameGroups = list(someNaming = design$someNaming, 
                  ourNamingGeneral = design$ourNamingGeneral,
                  ourNaming = design$ourNaming,
                  cellKind = design$cellKind,
                  justGaba = design$justGaba,
                  justGabaPV = design$justGabaPV)

fileNames = names(nameGroups)

stepi = 1
for (i in nameGroups[5:6]){
    dir.create(fileNames[stepi])
    groupNames = trimNAs(unique(i))
    realGroups = vector(mode = 'list', length = length(groupNames))
    names(realGroups) = groupNames
    for (j in 1:length(groupNames)){
        realGroups[[j]] = which(i == groupNames[j])           
    }
    groupAverages = list()
    
    
    for (j in realGroups){
        groupAverage = tryCatch({apply(exprData[,j], 1, median)},
                                error = function(cond){
                                    print('fuu')
                                    return(exprData[,j])
                                })
        groupAverages = c(groupAverages, list(groupAverage))
    }
    
    names(groupAverages)= groupNames
    groupAverages = (as.data.frame(groupAverages))
    
    datComb = specificity_index(2^groupAverages)
    rownames(datComb) = 1:nrow(datComb)
    print(stepi)
    for (k in 1:length(realGroups)){
        temp = datComb[which(datComb[,2*k]<1e-4),(2*k-2+1):(2*k)]
        if (is.null(dim(temp))){
            dim(temp) = c(1,2)
            rownames(temp) = rownames(datComb)[which(datComb[,2*k]<1e-4)]
        } else{
            temp = temp[order(temp[,1]), ]   
        }
        
        write.table(geneData[as.numeric(rownames(temp)),'Gene.Symbol' ], file = fileNames[stepi] + '/' + names(realGroups[k]),quote=F,row.names=F,col.names=F)
        
        
    }
    
    stepi= stepi + 1
        
}

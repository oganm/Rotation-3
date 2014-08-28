whichFile = 'Data/mostVariableExp'
outFile = 'Data/mostVariableQuantileNormalized'

quantileNorm = function(whichFile,outFile){
    require(preprocessCore)
    
    allDataPre = read.csv(whichFile, header = T)
    
    geneData = allDataPre[,1:3]
    exprData = allDataPre[,4:ncol(allDataPre)]
    newExprData = normalize.quantiles(as.matrix(exprData))
    boxplot(newExprData)
    newExprData = as.data.frame(newExprData)
    colnames(newExprData) = colnames(exprData)
    newAllData = cbind(geneData, newExprData)
    write.csv(newAllData, file = outFile, row.names=FALSE)
}

quantileNorm(whichFile,outFile)

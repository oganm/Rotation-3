getParent = function(step = 1){
  wd = getwd()
  setwd('..')
  parent = getwd()
  setwd(wd)
  return(parent)
}

dpaste = function (...){
  paste(..., sep='')
}

parent = getParent()
parent = dpaste(parent,'/') 

allDataPre = read.csv(dpaste(parent,'Data/allNormalized'), header = T)

geneData = allDataPre[,1:3]
exprData = allDataPre[,4:ncol(allDataPre)]

newExprData = vector(mode = 'double',
                     length = ncol(exprData) * length(unique(geneData$Gene.Symbol))
                    )

dim(newExprData) = c(length(unique(geneData$Gene.Symbol)), ncol(exprData))
colnames(newExprData) = colnames(exprData)

prog = txtProgressBar(min = 1, max = length(unique(geneData$Gene.Symbol)),style=3)

for (i in 1:length(unique(geneData$Gene.Symbol))){
  indexes = which(geneData$Gene.Symbol %in% unique(geneData$Gene.Symbol)[i])
  groupData = exprData[indexes,]
  chosen = which.max(apply(groupData,1,var))
  newExprData[i,] = as.double(groupData[chosen,])
  setTxtProgressBar(prog, i)
}
close(prog)

newExprData = as.data.frame(newExprData)

newGeneData = geneData[match(unique(geneData$Gene.Symbol),geneData$Gene.Symbol),]


newAllData = cbind(newGeneData, newExprData)
rownames(newAllData) = NULL
write.csv(newAllData, file = 'mostVariableExp', row.names=FALSE)


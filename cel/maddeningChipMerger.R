source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')

wd = getwd()
setwd('..')
parent = getwd()
setwd(wd)

findInList = function(object, aList){
  indexes = vector()
  for (i in 1:length(aList)){
    if (object %in% aList[[i]]){
      indexes = c(indexes, i)
    }
  }
  return(indexes)
}


require(affy)
insist(compare)
affydata <- ReadAffy()

rma

setwd('MOE430A')
affydataOld = ReadAffy()
setwd('..')

setpNListOld = probeNames(affydataOld)
pNList = probeNames(affydata)


subsetList = pNList[pNList %in% pNListOld]


subsetPm = pm(affydata, unique(subsetList))
subsetPmOldOrdered = pm(affydataOld, unique(subsetList))
subsetPmOld = pm(affydataOld)

allPm = cbind(subsetPm, subsetPmOldOrdered)

rownames(allPm) = rownames(subsetPmOld)





subset = NULL
verbose = TRUE
destructive = TRUE
normalize = TRUE
background = TRUE
bgversion = 2
ngenes = length(geneNames(affydataOld))
allpNList = split(0:(length(pNListOld) - 1), pNListOld)


#rownames(allPm) = 1:nrow(allPm)
exprs <- .Call("rma_c_complete", allPm, 
               allpNList, ngenes, normalize, background, bgversion, 
               verbose, PACKAGE = "affy")

phenoD = combine(phenoData(affydata), phenoData(affydataOld))
annot =  annotation(affydataOld)
protocolD = combine(protocolData(affydata), protocolData(affydataOld))
experimentD = experimentData(affydataOld)


newNormalized = new("ExpressionSet", phenoData = phenoD, annotation = annot, 
    protocolData = protocolD, experimentData = experimentD, 
    exprs = exprs)

#rm(affydata)
#rm(affydataOld)


newNormalized -> nvals

require(moe430a.db)
insist(gdata)

ned <- exprs(nvals)
nsamp <- sampleNames(nvals)
nprobes <- featureNames(nvals)

x = moe430aGENENAME
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

x <- moe430aSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
sdf <- data.frame(probe=names(vals), gene=vals)


sadf <- merge(sdf,adf, by="probe", all.x=TRUE, sort=FALSE)

colnames(sadf) <- c("Probe","Gene Symbol","Annotation")

aned <- merge(sadf,ned, by.x="Probe", by.y="row.names", all.x=TRUE, sort=FALSE)
header = gsub('.cel', '', gsub('.CEL','', colnames(aned)[4:ncol(aned)]))
colnames(aned) = c(colnames(aned)[1:3], header)
write.csv(aned, "allNormalized", row.names=FALSE)
boxplot(aned[,4:ncol(aned)])
df = read.xls("Design.xls")
gsms = regmatches(df[, 1], gregexpr("GSM\\d\\d\\d\\d\\d(\\d|)", df[, 1],perl=T))

colnames(aned) %in% gsms[[65]]


indexes = vector()
for (i in 1:length(header)){
   indexes = c(indexes, findInList(header[i], gsms))
}


newDesign = data.frame(sampleName = header, originalIndex = indexes, df[indexes,])
colnames(newDesign) = c('sampleName','originalIndex','officialName','cellType','anatomical','age','method','isolation','amplification','RNAamount','platform','reference')
newDesign$originalIndex = as.numeric(newDesign$originalIndex)

newDesign$age = gsub('~', '', newDesign$age)
newDesign$age = gsub('P', '', newDesign$age)
newDesign$age = gsub('7-8', '7.5', newDesign$age)
newDesign$age[grepl('(precise age not given)',newDesign$age)] = 60
newDesign$age = as.numeric(newDesign$age)
newDesign = newDesign[order(as.numeric(rownames(newDesign))),]

write.csv(newDesign, "normalizedDesignNonOrdered", row.names=FALSE)
aned['Gene Symbol',]


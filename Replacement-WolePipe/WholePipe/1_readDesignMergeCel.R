# Reads CEL files from a given directory according to a design file
#essential parts of the file are
#1. the first collumn must include the GSMs
#2. platform must be written to the Platform collumn

source('https://raw.githubusercontent.com/oganm/toSource/master/ogbox.r')
source('https://raw.githubusercontent.com/oganm/toSource/master/mergeChips.R')
require(affy)
require(compare)

desFile='Data2/Design.csv'
# just selects S1s of dopaminergics   (((?<=:)|(?<=[,]))A\d.*?S1_M430A)
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
celDir ='cel'
outFolder='Data2'

design = read.table(desFile,quote='',header=T,sep='\t')
gsms = regmatches(design[, 1], gregexpr(celRegex, design[, 1],perl=T))

platforms = unique(design$Platform)

affies = lapply( rep("AffyBatch", len(platforms)), new )
names(affies) = platforms

#create a vector of affy objects based on how many chips are there just works with 2 right now
#implementation of a third generation is straightforward. in fact i could have done it instead
#of writing this but oh well...
for (i in 1:len(affies)){
    celsInFolder = list.files(
        paste0(celDir,'/',platforms[i]))
    celsNoExtension = gsub('[.](C|c)(E|e)(L|l)','',celsInFolder)
    gsms = regmatches(design[design$Platform == platforms[i], 1], gregexpr(celRegex, design[design$Platform == platforms[i], 1],perl=T))
    gsms = unlist(gsms)
    
    relevant = celsInFolder[celsNoExtension %in% gsms]
    
    affies[i] = ReadAffy(filenames = paste0(celDir,'/',platforms[i],'/',relevant))
}

newNormalized = mergeChips(affies[[1]],affies[[2]])


nvals = newNormalized
require(mouse430a2.db)

ned <- exprs(nvals)
nsamp <- sampleNames(nvals)
nprobes <- featureNames(nvals)

x = mouse430a2GENENAME
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

x <- mouse430a2SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
sdf <- data.frame(probe=names(vals), gene=vals)


sadf <- merge(sdf,adf, by="probe", all.x=TRUE, sort=FALSE)

colnames(sadf) <- c("Probe","Gene Symbol","Annotation")

aned <- merge(sadf,ned, by.x="Probe", by.y="row.names", all.x=TRUE, sort=FALSE)
header = gsub('.cel', '', gsub('.CEL','', colnames(aned)[4:ncol(aned)]))
colnames(aned) = c(colnames(aned)[1:3], header)
write.csv(aned, paste0(outFolder,"/allNormalized"), row.names=FALSE)
boxplot(aned[,4:ncol(aned)])
gsms = regmatches(design[, 1], gregexpr(celRegex, design[, 1],perl=T))

#gsms = regmatches(df[, 1], gregexpr("(GSM\\d\\d\\d\\d\\d(\\d|))|(PC\\d....)|(Y+.*?((?=(,))|\\d+))|(((?<=:)|(?<=,))A\\d.*?30A)|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))", df[, 1],perl=T))


indexes = vector()
for (i in 1:length(header)){
    indexes = c(indexes, findInList(header[i], gsms))
}
header[!header %in% unlist(gsms)]
newDesign = data.frame(sampleName = header, originalIndex = indexes, design[indexes,])
colnames(newDesign) = c('sampleName','originalIndex',names(design))
#newDesign$originalIndex = as.numeric(newDesign$originalIndex)

#newDesign$age = gsub('~', '', newDesign$age)
#newDesign$age = gsub('P', '', newDesign$age)
#newDesign$age = gsub('7-8', '7.5', newDesign$age)
#newDesign$age[grepl('(precise age not given)',newDesign$age)] = 60
#newDesign$age = as.numeric(newDesign$age)
newDesign = newDesign[order(as.numeric(rownames(newDesign))),]

write.table(newDesign, paste0(outFolder,"/normalizedDesign"), row.names=FALSE,sep = '\t', quote=F)

rm(aned)
rm(ned)
rm(sadf)
rm(adf)
rm(affies)
rm(nvals)
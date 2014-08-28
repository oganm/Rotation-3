wd = getwd()
setwd('..')
parent = getwd()
setwd(wd)

require(gdata)

files = list.files()

df = read.xls( "Design.xls")

gsms = regmatches(df[, 1], gregexpr("GSM\\d\\d\\d\\d\\d(\\d|)", df[, 1],perl=T))


gsms = regmatches(df[, 1], gregexpr("(GSM\\d\\d\\d\\d\\d(\\d|))|(PC\\d....)|(Y+.*?((?=(,))|\\d+))|(((?<=:)|(?<=,))A\\d.*?30A)|(v2_(?![G,H,r]).*?((?=(,))|($)))", df[, 1],perl=T))




source("http://www.bioconductor.org/biocLite.R")

# biocLite installs or updates Bioconductor and CRAN packages, ensuring that packages 
# from the appropriate version of Bioconductor are installed, and that all packages
# remain up to date. 

biocLite("affy") 
biocLite("affyPLM")
#biocLite("mogene10stv1cdf")
#biocLite("org.Mm.eg.db")
#biocLite("mogene10sttranscriptcluster.db")
biocLite("mouse4302cdf")
biocLite("mouse4302.db")
biocLite("moe430a.db")
#biocLite("hgu133plus2cdf")
#biocLite("hgu133plus2.db")
#biocLite("hgu133acdf")
#biocLite("hgu133a.db")

library(affy)
library(affyPLM)
#library(mogene10stv1cdf)
#library(org.Mm.eg.db)
#library(mogene10sttranscriptcluster.db)
library(mouse4302cdf)
library (mouse4302.db)
require(moe430a.db)
#library(hgu133plus2cdf)
#library (hgu133plus2.db)
#library(hgu133plus2cdf)
#library (hgu133a.db)
#library(hgu133acdf)

# 1 Read in probe level data

# read the data. In order to read all the CEL files in the wd, use (), otherwise - 
#cel_files<-as.vector(cel_files[which(samples[1,]=="GPL96")])
#affydata <- ReadAffy(filenames = as.character(cel_files[which(samples[1,]=="GPL96")]))
affydata <- ReadAffy()

setwd('MOE430A')

affydataOld = ReadAffy()

setwd('..')

# raw expression data
ed <- exprs(affydata)
edOld = exprs(affydataOld)


samp <- sampleNames(affydata)
sampOld = sampleNames(affydataOld)
probes <- featureNames(affydata)
probesOld = featureNames(affydataOld)

newSubset = ed[probes %in% ,]


# 2 Normalizing Data   
#
# The Affy package has implementations of a number of normalization methods
# for single-channel arrays. This includes (among others):
#   - mas5() - Affymetrix's own normalization program
#   - rma() - 'Robust Multi-Chip' average
#   - gcrma() - A bias-corrected RMA
# GCRMA is good but takes significantly longer than RMA, so RMA is the
# most commonly used
nvals <- rma(affydata,affydataOld)

# normalised expression data
ned <- exprs(nvals)

nsamp <- sampleNames(nvals)
nprobes <- featureNames(nvals)

#Add Gene symbol and annotation

#match annotation
x <- mouse4302GENENAME
#x = moe430aGENENAME
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)


#match Gene Symbol

#x <- moe430aSYMBOL
x = mouse4302SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
sdf <- data.frame(probe=names(vals), gene=vals)


sadf <- merge(sdf,adf, by="probe", all.x=TRUE, sort=FALSE)

colnames(sadf) <- c("Probe","Gene Symbol","Annotation")

#Add all to the expression data

aned <- merge(sadf,ned, by.x="Probe", by.y="row.names", all.x=TRUE, sort=FALSE)

write.csv(aned, "oldChipSoloNorm", row.names=FALSE)


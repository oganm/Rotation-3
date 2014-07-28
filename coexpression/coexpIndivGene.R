# match found genes in mice to humans
# extends the marker network by coexpression
# 
# looks at coexpression changes of individual marker gene pairs
# only uses control samples for
require('dynamicTreeCut')
require(EBcoexpress)
library(topGO)
require(biomaRt)
require(ReactomePA)
outFold = 'ourNamingGeneral'
geneListLoc = 'C:/Users/Ogan/Dropbox/Rotation 3/AllGeneralComparison/Rest/InBetween/ourNamingGeneral'
source('C:/Users/Ogan/Dropbox/Rotation 3/ogbox.r')

# load design and expr ---------------------------------
design = read.table('C:/Users/Ogan/Dropbox/Rotation 3/Data/normalizedDesign.csv',header=T,sep='\t')
orthoInfo = read.csv('C:/Users/Ogan/Dropbox/Rotation 3/Data/HG-U133A.na33.ortholog.csv')
orthoInfo = orthoInfo[orthoInfo$Ortholog.Array == 'MOE430A',1:5 ]

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

# load human bipolar data ---------------------------------------
diagnosis = read.delim("C:/Users/Ogan/Dropbox/Rotation 3/Data/Stanley_Bipolar_Controls.csv",sep = "\t")

load('C:/Users/Ogan/Dropbox/Rotation 3/Data/GSE12649_GSE5388_expression.RData')

GSE12649 = aned_high_GSE12649
GSE5388 = aned_high_GSE5388
GSE12649Genes = GSE12649[,1:3]
GSE5388Genes = GSE5388[,1:3]
GSE12649Expr = GSE12649[,4:ncol(GSE12649)]
GSE5388Expr = GSE5388[,4:ncol(GSE5388)]
    


GSE12649Des = Samples_GSE12649
GSE5388Des = Samples_GSE5388
rm(aned_high_GSE12649)
rm(aned_high_GSE5388)
rm(Samples_GSE12649)
rm(Samples_GSE5388)


removeGSE12649 = which(GSE12649Des['characteristics',] %in% diagnosis$GSE12649[diagnosis$Axis.I.Primary.Diagnosis == 'BP II '])
removeGSE5388 = which(GSE5388Des['characteristics',] %in% diagnosis$GSE5388[diagnosis$Axis.I.Primary.Diagnosis == 'BP II '])

GSE12649Expr = GSE12649Expr[, -removeGSE12649]
GSE5388Expr = GSE5388Expr[, -removeGSE5388]

GSE12649Des = GSE12649Des[, -removeGSE12649]
GSE5388Des = GSE5388Des[, - removeGSE5388]

testGroups = list (
    GSE12649 = GSE12649Des['disease_state',],
    GSE5388 = gsub('_t','',GSE5388Des['disease_state', ])
    )


testExpr = list (
    GSE12649 = GSE12649Expr,
    GSE5388 = GSE5388Expr
    )
colnames(GSE12649Genes) = colnames(geneData)
colnames(GSE5388Genes) = colnames(geneData)
testGene = list(
    GSE12649 = GSE12649Genes,
    GSE5388 = GSE5388Genes)


# load a priori markers.----------------------------
filenames = list.files(geneListLoc,include.dirs = FALSE)
fileContents = lapply(geneListLoc + '/' + filenames, read.table)
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


# contract and expand marker gene network by using controls from both samples----------------------------
testExprC = vector (mode = 'list', length = length(testExpr))
testGeneC = vector (mode = 'list', length = length(testExpr))
exprDatas = vector (mode = 'list', length = length(testExpr))
geneDatas = vector (mode = 'list', length = length(testExpr))
mart <- useMart("ensembl",dataset = 'hsapiens_gene_ensembl')
attributes <- c("entrezgene", "hgnc_symbol")

for (i in 1: length(testExpr)){
    usePairs = orthoInfo[{tolower(orthoInfo$Ortholog.Probe.Set) %in% tolower(geneData$Probe)}&
                         {tolower(orthoInfo$Probe.Set.ID) %in% tolower(testGenes[[i]]$Probe)}
                        ,c(1,3)]

    usePairs = usePairs[!(duplicated(usePairs[,1]) | duplicated(usePairs[,2])),]
    exprDatas[[i]] = exprData[match(tolower(usePairs$Ortholog.Probe.Set), tolower(geneData$Probe)), ] 
    geneDatas[[i]] = geneData[match(tolower(usePairs$Ortholog.Probe.Set), tolower(geneData$Probe)), ]
    
    testExprC[[i]] = testExpr[[i]][match(usePairs$Probe.Set.ID, testGene[[i]]$Probe), ] 
    testGeneC[[i]] = testGene[[i]][match(usePairs$Probe.Set.ID, testGene[[i]]$Probe), ]
    
    allIDs = getBM(attributes=attributes, filters = "hgnc_symbol", values = testGeneC[[i]]$Gene.Symbol, mart = mart, uniqueRows=T)
    
    for (j in 1:length(filenames)){
        #### attempt to expand selection--------------
        #just takes the data from control to avoid differential expression artefacts
        #         corCheck = testExprC[[i]][geneDatas[[i]]$Gene.Symbol %in% foldGenes[[which(names(foldGenes)==filenames[j])]],testGroups[[i]]=='Cont']
        #         rownames(corCheck) = testGeneC[[i]]$Gene.Symbol[geneDatas[[i]]$Gene.Symbol %in% foldGenes[[which(names(foldGenes)==filenames[j])]]]
        # 
        #      
        #         net = blockwiseModules(t(corCheck), power = 6,
        #                                TOMType = "unsigned", minModuleSize = 10,
        #                                reassignThreshold = 0, mergeCutHeight = 0.25,
        #                                numericLabels = TRUE, pamRespectsDendro = FALSE,
        #                                verbose = 3)
        #      
        #         chosen = mode(net$colors)
        #         smallSet = rownames(corCheck)[net$colors == chosen]
        #                  
        #         netFull = blockwiseModules(t(testExpr[[i]][, testGroups[[i]]=='Cont']), power = 6,
        #                                TOMType = "unsigned", minModuleSize = 30,
        #                                reassignThreshold = 0, mergeCutHeight = 0.25,
        #                                numericLabels = TRUE, pamRespectsDendro = FALSE,
        #                                verbose = 3)
        #         
        #         which(testGene[[i]]$Gene.Symbol %in% smallSet)
        #         
        #         
        # 
        #         # did not took the absolute value while clustering as we are not looking for negatively corrolated ones
        #         corr = cor(t(corCheck)^6)
        #         distMat = dist(corr)
        #         hc = hclust(distMat)
        #         clusters = cutreeDynamic(hc , minClusterSize = 0, distM = as.matrix(distMat), minSplitHeight = 0.5)
        #         # some test code
        #         #mean(apply(corCheck[clusters ==1,], 2,var))
        #         #mean(abs(corr[clusters == ,]))
        #         #sum(clusters==1)
        #         
        #         #take the largest cluster
        #         chosens = rownames(corCheck)[clusters == mode(clusters)]
        #         
        #         
        #         #now to expand the set
        #         corr = cor(t(testExpr[[i]])^6)
        #         distMat = dist(corr)
        #         
        ######
        #change this part if you manage to get the previous park work
        corCheck = testExprC[[i]][geneDatas[[i]]$Gene.Symbol %in% foldGenes[[which(names(foldGenes)==filenames[j])]], testGroups[[i]]!='SCZ']
        rownames(corCheck) = testGeneC[[i]]$Gene.Symbol[geneDatas[[i]]$Gene.Symbol %in% foldGenes[[which(names(foldGenes)==filenames[j])]] ]
        
       
        
        condition = as.numeric(as.factor(testGroups[[i]][testGroups[[i]]!='SCZ']))
        
        median(abs(cor(t(corCheck[, condition == 1]))))
        
        pattern =ebPatterns(c("1,1", "1,2"))
        # BWMC is used. apparently better than pearson
        D = makeMyD(corCheck, condition, useBWMC = T) 
        #sum(abs(D[,1]))
        #sum(abs(D[,2]))
        initHP <- initializeHP(D, condition, seed = 1)
        oneStep = ebCoexpressOneStep(D, condition, pattern, initHP)
        tail(oneStep[[2]][order(oneStep[[2]][,2]),])
        genePairs = rownames(D[oneStep[[2]][,'DC1']>0.95,])
        geneSingles = unique(unlist(strsplit(genePairs, split = '~')))
        # get probe names for GO
        probes = testGeneC[[i]]$Probe[testGeneC[[i]]$Gene.Symbol %in% geneSingles]
        background =  testGeneC[[i]]$Probe[testGeneC[[i]]$Gene.Symbol %in% rownames(corCheck)]
        allBack = testGeneC[[i]]$Probe
        
        # ID conversation for enrichment analysis
        IDConv = getBM(attributes=attributes, filters = "hgnc_symbol", values = rownames(corCheck), mart = mart, uniqueRows=T)
        
        #look at enrichment with other cell type specific ones as control
        pathEnrich = enrichPathway(IDConv$entrezgene[IDConv$hgnc_symbol %in% geneSingles], universe = trimNAs(IDConv$entrezgene) ,pvalueCutoff = 0.05, readable = T)
        if (nrow(summary(pathEnrich)) != 0){
            png('enrichment/' +'DC_' + names(testExpr)[i] + '_' + filenames[j] + '.png', width = 800)
            print(barplot(pathEnrich, showCategory = 10))
            dev.off()
            write.table(summary(pathEnrich), file = 'enrichment/' + 'DC_' + names(testExpr)[i] + '_' + filenames[j] )
        }else {
            print('fuuu1')
        }
        #look at enrichment without background of cell type specific stuff
        pathEnrich = enrichPathway(IDConv$entrezgene[IDConv$hgnc_symbol %in% geneSingles], universe = trimNAs(allIDs$entrezgene))
        if (nrow(summary(pathEnrich)) != 0){
            png('enrichment/' +'DC_'+'NoBack_'+ names(testExpr)[i] + '_' + filenames[j] + '.png', width = 800)
            print(barplot(pathEnrich, showCategory = 10))
            dev.off()
            write.table(summary(pathEnrich), file = 'enrichment/' + 'DC_' +'NoBack_'+ names(testExpr)[i] + '_' + filenames[j] )
        }else {
            print('fuuu2')
        }
        
        #look at erichment of cell type specific genes
        pathEnrich = enrichPathway(trimNAs(IDConv$entrezgene), universe = trimNAs(allIDs$entrezgene))
        if (nrow(summary(pathEnrich)) != 0){
            print('whaa')
            png('enrichment/' +'CellType_'+ names(testExpr)[i] + '_' + filenames[j] + '.png', width = 800)
            print(barplot(pathEnrich, showCategory = 10))
            print('whoo')
            dev.off()
            write.table(summary(pathEnrich), file = 'enrichment/' + 'CellType_' + names(testExpr)[i] + '_' + filenames[j] )
        }else {
            print('fuuu3')
        }
        
        # go enrichment with cell type with background
#         geneList = rep(1, length(background))
#         geneList[background %in% probes]=0
#         names(geneList) = background
#         
#         geneList = rep(1, length(allBack))
#         geneList[allBack %in% probes]=0
#         names(geneList) = allBack
#         
#         geneList = rep(1, length(allBack))
#         geneList[allBack %in% background]=0
#         names(geneList) = allBack
#         
#         sampleGOdata <- new("topGOdata",
#                             description = "Simple session", ontology = "BP",
#                             allGenes = geneList, geneSel = topDiffGenes,
#                             nodeSize = 10,
#                             annot = annFUN.db, affyLib = affyLib)
        
        
    }
    
    
}

# unfinished
    
    
#     
 ebPatterns(c("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1",
               "1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2
2 2 2"), TRUE)
    


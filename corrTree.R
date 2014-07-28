require(ape)
require(gdata)
#transcriptomeSources = read.xls('journal.pone.0016493.s010.xls')
transcriptomeSources = read.xls('journal.pone.0016493.s010.xls',perl='C:/Perl64/bin/perl.exe')
allData = t(read.csv('Data/allData', header = F))
names = read.delim('Data/colNames', header = F)
usedData = t(read.csv('Data/usedData', header = F))
transcriptomeSources$No=1:64
fullNames=do.call(paste,transcriptomeSources[c('No','Description', 'Anatomical.Region')])


correlation = cor(usedData, method = "pearson")
#rownames(correlation)=fullNames
#colnames(correlation)=fullNames
cluster=hclust(dist(t(allData)))
cluster=hclust(dist(correlation))

colnames(correlation)=NULL
rownames(correlation)=NULL

plot(cluster, hang = -1)


#colors to be used. must exceed to no of longest list of labels
useColors = c('green',
              'black',
              'blue',
              'orange',
              'brown',
              'red',
              'darkgreen',
              'yellow',
              'darkgray',
              'cyan',
              'darkseagreen',
              'deeppink',
              'aquamarine4',
              'darkkhaki',
              'darkolivegreen',
              'red4',
              'tan2',
              'thistle',
              'slateblue4',
              'thistle',
              'powderblue',
              'lightslategray',
              'maroon',
              'thistle4'
              )
###################################################################################################
reference = as.character(transcriptomeSources$Reference)
refColor = rep('white', 64)

refNames = unique(reference)

for (i in 1:length(refNames)){
  refColor[grepl(refNames[i], reference)] = useColors[i]
}


legend("bottomleft", legend= refNames,
       fill = c(useColors))


################################################################################################
method = as.character(transcriptomeSources$Method)
methodColor = rep('white', 64)
methodNames = unique(method)

for (i in 1:length(methodNames)){
  methodColor[grepl(methodNames[i], method)] = useColors[i]
}

legend("bottomleft", legend= methodNames,
       fill = c(useColors))

###############################################################################################

anatomical = as.character(transcriptomeSources$Anatomical.Region)
anatomical = gsub('Corpus Striatum', 'Striatum', anatomical)
anatomical[grepl('Somatosensory',anatomical)] = 'Somatosensory cortex'

anatomicalColor = rep('white', 64)
anatomicalNames = unique(anatomical)



for (i in 1:length(anatomicalNames)){
  anatomicalColor[grepl(anatomicalNames[i], anatomical)] = useColors[i]
}

legend("bottomleft", legend= anatomicalNames,
       fill = c(useColors))


#################################################################################################

age = as.character(transcriptomeSources$Age.of.Mouse..postnatal.day.)
age = gsub('~', '', age)
age[grepl('P60',age)] = 'Adult'
age[grepl('(precise age not given)',age)] = 'Adult'

ageColor = rep('white', 64)
ageNames = unique(age)



for (i in 1:length(ageNames)){
  ageColor[grepl(ageNames[i], age)] = useColors[i]
}

legend("bottomleft", legend= ageNames,
       fill = c(useColors))

#################################################################################################

platform = as.character(transcriptomeSources$Microarray.Platform)


platformColor = rep('white', 64)
platformNames = unique(platform)



for (i in 1:length(platformNames)){
  platformColor[grepl(platformNames[i], platform)] = useColors[i]
}

legend("bottomleft", legend= platformNames,
       fill = c(useColors))

###################################################################################################

isolation = as.character(transcriptomeSources$RNA.isolation.method)


isoColor = rep('white', 64)
isoNames = unique(isolation)



for (i in 1:length(isoNames)){
  isoColor[grepl(isoNames[i], isolation)] = useColors[i]
}

legend("bottomleft", legend= isoNames,
       fill = c(useColors))
###################################################################################################

description = as.character(transcriptomeSources$Description)
description[grepl('^Neurons', description)] = 'Mixed Neurons'
description[grepl('Medium Spiny Neurons', description)] = 'Medium Spiny Neurons'
description[grepl('^Cholinergic', description)] = 'Cholinergic'
description[grepl('Oligodendrocyte', description)] = 'Oligodendrocytes'
description[grepl('Dopaminergic', description)] = 'Dopaminergic'
description[grepl('Motor', description)] = 'Motor+Cholinergic'

description=gsub(', P.*?$','',description)

desColor = rep('white', 64)
desNames = unique(description)



for (i in 1:length(desNames)){
  desColor[grepl(desNames[i], description)] = useColors[i]
}

legend("bottom", legend= desNames,
       fill = c(useColors))


allColors=cbind(reference=refColor, method=methodColor, anatomical=anatomicalColor, platform=platformColor, description=desColor)

palette = colorRampPalette(c("orange", "white"))(n = 1000)
a = heatmap.3(allData, trace = "none", Rowv = F, Colv = T, 
          col = palette, ColSideColors = as.matrix(allColors))
plot.new()
title('description')
legend("center", legend= desNames,
       fill = c(useColors))
#plot.new()
#title('isolation')
#legend("bottomleft", legend= isoNames,
#       fill = c(useColors))
plot.new()
title('platform')
legend("center", legend= platformNames,
       fill = c(useColors))
#plot.new()
#title('age')
#legend("bottomleft", legend= ageNames,
#       fill = c(useColors))
plot.new()
title('anatomical')
legend('center', legend= anatomicalNames,
       fill = c(useColors))
plot.new()
title('method')
legend("center", legend= methodNames,
       fill = c(useColors))
plot.new()
title('reference')
legend("center", legend= refNames,
       fill = c(useColors))


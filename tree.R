require(ape)

mouse=read.table('namesMIR')
mouse=as.character(t(mouse))
mouse=tolower(mouse)
  
mytree <- read.tree(file.choose())

mytree$tip.label=tolower(gsub("id.*?_",'',mytree$tip.label))
#temp is the list of capped from the other program

capped=gsub('-','_',temp)


plot.phylo(mytree,no.margin=T,tip.color=colors,direction='downwards')

relaxedCapped=tolower(gsub('-','_',rownames(evenMoreInterestin)))


which(mytree$tip.label %in% mouse)
which(mytree$tip.label %in% capped)

colors=rep('black',length(mytree$tip.label))

colors[which(mytree$tip.label %in% mouse)]='green'
colors[which(mytree$tip.label %in% tolower(capped))]='red'
colors[which(mytree$tip.label %in% relaxedCapped)]='blue'

plot.phylo(mytree,no.margin=T,tip.color=colors,direction='downwards')

gsub("id.*?_",'',mytree$tip.label)
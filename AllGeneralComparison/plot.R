require(plot3D)
require(scatterplot3d)
library(rgl)

#excecute humanBipol.R to get useBipolarExpr1 and 2


open3d()
pca = prcomp(t(useBipolarExpr1), scale = T)

a = pca$x

colooor = vector(length = length(rownames(a)))
colooor[grep('SCZ', rownames(a))] = 'blue'
colooor[grep('Cont', rownames(a))] = 'green'
colooor[grep('BP', rownames(a))] = 'red'

plot3d(x=  pca$x[,1], y =  pca$x[,2], z=  pca$x[,3],size= 8, col = colooor)


insist(ggplot2)


pca1 = prcomp(t(useBipolarExpr1), scale = T)
pca2 = prcomp(t(useBipolarExpr2), scale = T)
pca3 = prcomp(t(useSelectExprData1), scale = T)

for (i in 1:ncol(pca2$rotation)){
  if (sign(pca1$rotation['14866', i]) != sign(pca2$rotation['14866', i])){
    pca2$rotation[, i] = -pca2$rotation[, i]
    pca2$x[, i] = - pca2$x[ ,i]
  }
}

for (i in 1:ncol(pca2$rotation)){
  if (sign(pca1$rotation['14866', i]) != sign(pca3$rotation['12802', i])){
    pca3$rotation[, i] = -pca3$rotation[, i]
    pca3$x[, i] = - pca3$x[ ,i]
  }
}


colooor = vector(length = length(rownames(pca1$x)))
colooor[grep('SCZ', rownames(pca1$x))] = 'SCZ'
colooor[grep('Cont', rownames(pca1$x))] = 'Cont'
colooor[grep('BP', rownames(pca1$x))] = 'BP'
pca1$x = cbind(pca1$x, colooor)
a = as.data.frame(pca1$x)
a$PC1 =as.numeric(as.character(a$PC1))
a$PC2 =as.numeric(as.character(a$PC2))

cbPalette = c('red', 'black', 'green')
(ggplot(as.data.frame(a), aes(x = PC1, y = PC2, col = colooor)) + geom_point(size=3)+geom_vline(xintercept = 0)+ geom_hline(yintercept = 0)+scale_colour_manual(values=cbPalette))


colooor = vector(length = length(rownames(pca2$x)))
colooor[grep('SCZ', rownames(pca2$x))] = 'SCZ'
colooor[grep('Cont', rownames(pca2$x))] = 'Cont'
colooor[grep('BP', rownames(pca2$x))] = 'BP'
pca2$x = cbind(pca2$x, colooor)


b = as.data.frame(pca2$x)
b$PC1 =as.numeric(as.character(b$PC1))
b$PC2 =as.numeric(as.character(b$PC2))
(ggplot(as.data.frame(b), aes(x = PC1, y = PC2, col = colooor)) + geom_point(size=3)+geom_vline(xintercept = 0)+ geom_hline(yintercept = 0)+scale_colour_manual(values=cbPalette))


c = as.data.frame(pca3$x)
c = cbind(c, selectDes$ourNaming)
c$PC1 =as.numeric(as.character(c$PC1))
c$PC2 =as.numeric(as.character(c$PC2))
(ggplot(c, aes(x = PC1, y = PC2, col = selectDes$ourNaming)) + geom_point(size=3)+geom_vline(xintercept = 0)+ geom_hline(yintercept = 0))#+scale_colour_manual(values=cbPalette))



bipolar1 = useBipolarExpr1[,1:28]
control1 = useBipolarExpr1[,29:(29+30)]
sch1 =useBipolarExpr1[,(29+30+1):(29+30+35)]






library(ggplot2)
vps <- baseViewports()
pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot

coloDendro = function (data, sampleColors, maskNames = F, w = 20, h = 7){
  options(device = "windows")
  dev.new(width = w, height = h)
  #sampleColors takes in a data frame with titles in rows and samples to collumns
  
  if (maskNames == T){
    rownames(data)= NULL
  }
  hc = hclust(dist(data))
  
  
  par(mfrow=c(2,1))
  par(mar=c(0, 3.6/10*w, 0, 0))
  plot(as.dendrogram(hc), xlab="",sub="",axes = F , yaxs = "i")
  
  
  imageMatr = matrix(0, nrow = nrow(sampleColors), ncol = ncol(sampleColors))
  
  allColors = vector()

  for (i in 1:nrow(sampleColors)){
    allColors =c (allColors,unique(sampleColors[i,]))
  }
  
  allColors = unique(allColors)
  for (i in 1:nrow(sampleColors)){
    for (j in 1:ncol(sampleColors)){
     imageMatr[i,j] = which(allColors == sampleColors[i,j])
    }
  }
  
  
  par(mar=c(0, 5.05/10*w,0 , 1.45/10*w))
  
  image(t(imageMatr)[hc$order,nrow(imageMatr):1], col = allColors, axes = F)
    
  text(par("usr")[1], seq(1,0,-1/(nrow(sampleColors)-1)), labels = rownames(sampleColors), srt = 0, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
  
}


(p1)
p2 = ggplot

require(stats)

hc=hclust(dist(a))
as.dendrogram(hc)
plot(hc)

par(mfrow=c(2,1))
par(mar=c(0, 0, 0, 0))
plot(hc, xlab="",sub="", hang=-1)
image(matrix(runif(10*10,1,10),nrow=10),col=rainbow(10), axes = FALSE)

text(par("usr")[1], seq(0,1,1/4), labels = c('ggdhs','fhshfdhfd','fhfdhfdh','fhfdhfdhfd'), srt = 0, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
axis(2)

m <- matrix(rnorm(100), ncol=10) 
image(m) 
text(0:9/9, 0:9/9, 0:9,srt = 0) 

plot(as.phylo(hc), cex = 0.9, label.offset = 1)

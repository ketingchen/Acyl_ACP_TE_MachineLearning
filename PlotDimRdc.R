require(RColorBrewer)
require(car)
#relabel=function(data=isu.traits, traits="SL", cluster=kmeans.isu$cluster){
#  medians=aggregate(data[,traits]~cluster,FUN=median); colnames(medians)[2]="median"
#  cl.max=medians$cluster[which.max(medians$median)]
#  cl.min=medians$cluster[which.min(medians$median)]
#  cl.new=rep(2, length(cluster))
#  cl.new[cluster==cl.max]=3
#  cl.new[cluster==cl.min]=1
#  return(cl.new)  
#}

color.category=function(category, col=NULL){
  ncolors=nlevels(factor(category))
  if(is.null(col)){
    if(ncolors==2) colors=c("black","red") else colors=brewer.pal(ncolors, name="Dark2")
    return(colors[factor(category)])
  } else return(col[factor(category)])
}


##### scores: dataframes or matrices for principal components, with PCs as columns and sample id as rows
##### dims: the index of principal components to be plotted
##### colorBy: a vector with the same number of rows as in scores, specifyin the colors assigned for every sample
##### col: the colors to be used, and it can be generated automatically when leave this argument as default (NULL)
##### explained.variance: the explained proportion of variance (<1) by the model. For PLS, it should be R2Y
##### legend: a vector assigning labels for different colors and generating the corresponding legend
##### Ellipse: 95% confidence ellipse labeling the data points distribution with 95% probability
##### EllipseBy: a vector assigning groups for different ellipses
##### textBy: a vector with the same number of rows as in scores, specifying the texts that are used as symbols for every sample
##### legend.location: specify the location of legend either by:
#####                  "topright","center","topleft","bottomright","bottonleft"
#####                  or by coordinates at x and y axis

plotDimRdc=function(scores, dims=c(1:2), colorBy, col=NULL, explained.variance=NULL, legend=NULL, legend.location="topright", main=NULL, Ellipse=F, EllipseBy=NULL, textBy=NULL, line=3, line.axis=-0.5, symbol.size=1){
  scores=as.data.frame(scores)[,dims]
  rng=abs(t(apply(scores, 2, range)))*1.2*sign(t(apply(scores, 2, range)))
  colorScheme=color.category(colorBy, col=col)  
  
  colnames(scores)=paste("PC",1:ncol(scores),sep="")
  
  plot(scores, type="n",xlab="", xlim=rng[1,], ylim=rng[2,], ylab="", axes=F);box()
  axis(1, label=F, tck=-0.02)
  axis(1, line=line.axis, tick=F)
  axis(2, label=F, tck=-0.02)
  axis(2, line=line.axis, tick=F) 
  
  if(!is.null(explained.variance)){
    mtext(side=1, text=paste("PC", dims[1]," (", round(explained.variance[1]*100),"%)", sep=""), cex=1, line=line)
    mtext(side=2, text=paste("PC", dims[2]," (", round(explained.variance[2]*100),"%)", sep=""), cex=1, line=line)
  } 
    
  if(is.null(textBy)) points(scores, col=colorScheme,pch=19, cex=symbol.size) else text(scores, col=colorScheme, labels=textBy,cex=symbol.size)
  
  if(!is.null(legend)) legend(legend.location, fill=unique(colorScheme), legend=legend, bty="n") 
  
  if(!is.null(main)) mtext(main,cex=1)
  
  if(Ellipse==T)
  {
    sources=unique(EllipseBy)
    for(i in sources)
    {
    with(scores, dataEllipse(PC1[EllipseBy==i], PC2[EllipseBy==i], levels=.95, plot.points=FALSE, center.pch=F, col="black", lwd=0.5))
    }
  }
}
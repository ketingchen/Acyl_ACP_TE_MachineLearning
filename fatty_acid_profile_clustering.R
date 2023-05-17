#################################################################
################## Clustering based on FA profiles ##############
############################################## ##################
setwd("/Users/kchen/Box/SILKY/Keting/Projects/Fuyuan")
require(dynamicTreeCut)
data <- read.csv("pairwise comparison of TE sequences encoded with labels.csv", row.names=1)
substrate <- read.csv("TE substrate specificities.csv", row.names=1)
fa <- read.csv("TE fatty acid profiles.csv", row.names=1)

d <- as.matrix(dist(scale(fa)))
hc <- hclust(d=dist(scale(fa)), method="ward.D"); den.h <- as.dendrogram(hc)
cls.h <- cutreeDynamic(hc, method="hybrid",distM=d, cutHeight=50, minClusterSize=10)
names(cls.h) <- hc$labels

### The results for cutreeDynamc is not stable, as the cluster label often changes but clustering membership does not change
### here the cluster labels are forced to be consistent for each indepdent run
cls.h[cls.h==cls.h["rTE28"]] <- 33
cls.h[cls.h==cls.h["TEGm162"]] <- 22
cls.h[cls.h==cls.h["CpFatB1"]] <- 11
cls.h[cls.h==33] <- 3
cls.h[cls.h==22] <- 2
cls.h[cls.h==11] <- 1

### the scripts below plot the cluster dendrogram, it is redundant with the next block, and can be ignored
require(dendextend)
require(RColorBrewer)
#cls.kmeans <- paste0("Cls", cls.kmeans); names(cls.kmeans) <- rownames(fa)
colorCodes <- brewer.pal(3, "Dark2")
#colorCodes <- c(Cls1="black", Cls2="blue")
 names(colorCodes) <- c("3", "2", "1")
labels_colors(den.h) <- colorCodes[cls.h][hc$order]

par(mar=c(3, 5, 0.5, 5),mfrow=c(1,1), cex=0.7)
plot(den.h, horiz=T)


##### principal component analysis for fatty acid profiles, and dendrogram ####
##### the remaining scripts in this file generate fig 4 in the manuscript
fa <- read.csv("TE fatty acid profiles.csv", row.names=1)
fa <- data.frame(cluster=cls.h, fa)

##### load two customized scripts for PCA and PCA plotting
source("/Users/kchen/Box/SILKY/Keting/Projects/Essential scripts/PlotDimRdc.R")
source("/Users/kchen/Box/Silky/Keting/Projects/Essential scripts/PLS.R")
pca = pca.nipals(fa[,-1], scale=T)
scores = data.frame(pca$Scores)

## PCA plot
pdf("PCA_cluster_fa_v4.pdf", width=7.09, height=7.28, pointsize=7)
m=matrix(c(1, 1, 1, 1, 2, 3, 4, 5), 4, 2)
layout(m, heights=rep.int(0.8, nrow(m)))
par(cex=1, ps=7, mar=c(2, 1, 2, 12), oma=c(1, 1, 1, 1))
plot(den.h, horiz=T)
text(25, 0, labels="Euclidean distance")
arrows(x0=-17, y0=1, x1=-17, y1=25, length=0, xpd=T)
arrows(x0=-17, y0=26, x1=-17, y1=47, length=0, xpd=T)
arrows(x0=-17, y0=48, x1=-17, y1=57, length=0, xpd=T)
mtext(expression(bold("Cluster C")), side=4, at=12.5, adj=0, line=4.4, col=colorCodes[3])
mtext(expression(bold("Cluster B")), side=4, at=35.5, adj=0, line=4.4, col=colorCodes[2])
mtext(expression(bold("Cluster A")), side=4, at=50.5, adj=0, line=4.4, col=colorCodes[1])
mtext(expression(bold("a")), at=40, line=0, cex=1)
#legend("topleft", paste0("Cluster ",sort(unique(fa$cluster))),fill=rev(colorCodes))

par(mar=c(1.5, 4, 1, 1))
# plotDimRdc(scores, colorBy=fa$cluster, col=colorCodes, explained.variance=c(0.43, 0.16), legend=NULL, legend.location="topleft",textBy=rownames(fa), Ellipse=T, EllipseBy=fa$cluster, line=2, symbol.size=1)
plotDimRdc(scores, colorBy=fa$cluster, col=colorCodes, explained.variance=c(0.43, 0.16), legend=NULL, legend.location="topleft", Ellipse=F, EllipseBy=fa$cluster, line=2, line.axis=-0.3, symbol.size=1)

sources=sort(unique(fa$cluster), decreasing=T)
for(s in sources){
	with(scores, dataEllipse(PC1[fa$cluster==s], PC2[fa$cluster==s], levels=.95, plot.points=FALSE, center.pch=F, col=colorCodes[s], lwd=0.5))
}
text(-3.5, -1, labels=expression(bold("Cluster C")), col=colorCodes[3])
text(-1, 3.7, labels=expression(bold("Cluster B")), col=colorCodes[2])
text(-0.2, -3, labels=expression(bold("Cluster A")), col=colorCodes[1])
legend("topleft", legend=expression(bold("b")), bty="n")
#mtext(expression(bold("b")), at=-6, line=0, cex=1)


## bar chart
avg <- as.matrix(aggregate(.~cluster, data=fa, FUN=mean))
se <- as.matrix(aggregate(.~cluster, data=fa, function(x)sd(x)/sqrt(length(x))))
a <- gsub("\\.", ":", colnames(avg)[-1])
a <- gsub("C", "", a)
colnames(avg)[-1] <- colnames(se)[-1] <- a


par(mar=c(0, 4, 2, 1))
x1=barplot(avg[1,-1], col=colorCodes[1], ylim=c(0,67), las=2, axisname=F, axes=F);box()
errbar(x1, y=avg[1,-1], bar=se[1,-1])
axis(1, labels=F, at=x1, tck=-0.02)
axis(2, labels=F, tck=-0.02)
axis(2, at=c(0, 20, 40, 60), tick=F, line=-0.3)
mtext("Concentration (% mol)", side=2, line=2)
#mtext("Concentration", side=2, line=2)
#mtext(expression(bold("c")), at=-2.6, line=0, cex=1)
legend("topleft", legend=expression(bold("c")), bty="n")

par(mar=c(0, 4, 0, 1))
x2=barplot(avg[2,-1], col=colorCodes[2], ylim=c(0,67), las=2, axisname=F, axes=F);box()
errbar(x2, y=avg[2,-1], bar=se[2,-1])
axis(1, labels=F, at=x2, tck=-0.02)
axis(2, labels=F, tck=-0.02)
axis(2, at=c(0, 20, 40, 60), tick=F, line=-0.3)
mtext("Concentration (% mol)", side=2, line=2)
#mtext("Concentration", side=2, line=2)
#mtext(expression(bold("d")), at=-2.6, line=0, cex=1)
legend("topleft", legend=expression(bold("d")), bty="n")

par(mar=c(4, 4, 0, 1))
x3=barplot(avg[3,-1], col=colorCodes[3], ylim=c(0,67), las=2, axes=F);box()
errbar(x3, y=avg[3,-1], bar=se[3,-1])
axis(1, labels=F, at=x3, tck=-0.02)
axis(2, labels=F, tck=-0.02)
axis(2, at=c(0, 20, 40, 60), tick=F, line=-0.3)
mtext("Concentration (% mol)", side=2, line=2)
#mtext("% mol FA", side=2, line=1)
#mtext("FA concentration", side=2, line=2)
#mtext(expression(bold("e")), at=-2.6, line=0, cex=1)
legend("topleft", legend=expression(bold("e")), bty="n")
title(xlab="Fatty acid")

dev.off()

fuyuan.label <- data$label
names(fuyuan.label) <- rownames(data)
pir1 <- names(fuyuan.label)[fuyuan.label==1]
enz1 <- lapply(pir1, function(p)strsplit(p, "->")[[1]])
enz1 <- lapply(enz1, function(e)cls.h[e])
enz.diff1 <- do.call(c, lapply(enz1, function(e)abs(diff(e))))

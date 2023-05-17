require(ranger)
require(dynamicTreeCut)
setwd("~/Box/SILKY/Keting/Projects/Fuyuan")
data <- read.csv("pairwise comparison of TE sequences encoded with labels.csv", row.names=1)
#substrate <- read.csv("TE substrate specificities.csv", row.names=1)
fa <- read.csv("TE fatty acid profiles.csv", row.names=1)

d <- as.matrix(dist(scale(fa)))
hc <- hclust(d=dist(scale(fa)), method="ward.D"); den.h <- as.dendrogram(hc)
cls.h <- cutreeDynamic(hc, method="hybrid",distM=d, cutHeight=50, minClusterSize=10)
names(cls.h) <- hc$labels

cls.h[cls.h==cls.h["rTE28"]] <- 33
cls.h[cls.h==cls.h["rTE40"]] <- 22
cls.h[cls.h==cls.h["rTE48"]] <- 11
cls.h[cls.h==33] <- 3
cls.h[cls.h==22] <- 2
cls.h[cls.h==11] <- 1

### redefine the response for each instance  
### each instance is a pairwise comparison of two TEs
### the features are the sequence variation between two TEs at each amino acid position (0, same amino acid; 1, different amino acids)
### the response is defined based on the clustering analysis of the fatty acid profile (1, two TEs belonging to the same cluster; 0, different clusters)

label = rep(1, nrow(data)); names(label)=rownames(data)
for(i in 1:length(label)){
	nms <- strsplit(names(label)[i], "->")[[1]]
	c <- cls.h[nms]
	if(c[1]!=c[2]) label[i]=0
}


### random forest with new data label
data1 <- data; data1$label = label

##### prepare the output matrix
##### the random forest classifier are implemented 10 times to account for the randomness in the classifier construction
##### imp.corr.matrix records the importance score of all features for each classifier
##### pvalue.matrix records the corresponding pvalues
##### pred.err compares the predictive performance using the old datasets (with old response labels) and the new datasets (with the newly defined response labels)
##### the results of pred.err were not used in the manuscripts; it is for other in-house training purposes, and can be ignore for this study
#imp.matrix <- matrix(0, 10, 351, dimnames=list(paste0("RF",1:10), colnames(data)[1:351]))
imp.corr.matrix <- matrix(0, 10, 351, dimnames=list(paste0("RF",1:10), colnames(data)[1:351]))
pvalue.matrix <- matrix(0, 10, 351, dimnames=list(paste0("RF",1:10), colnames(data)[1:351]))
pred.err <- matrix(0, 2, 10, dimnames=list(c("Original", "Transformed"), paste0("Run",1:10)))

pb <- txtProgressBar(min=1, max=10, style=3)
for(i in 1:10){
	#### here two classifiers are established, one using the old response labels (rf0), and the other one uses the new response labels (rf)
	#### only the results from rf are used in the publication
	rf0 <- ranger(label~., data=data, importance="impurity_corrected", classification=T)
	rf <- ranger(label~., data=data1, importance="impurity_corrected", classification=T)
	imp.corr.matrix[i,] <- importance(rf)
	pvalue.matrix[i,] <- importance_pvalues(rf)[,2]
	pred.err[1,i] <- rf0$prediction.error
	pred.err[2,i] <- rf$prediction.error
	setTxtProgressBar(pb, i)
}

##### set the zero pvalue to a minimum value for pvalue correction
pvalue.matrix[pvalue.matrix==0] <- 1e-167
pvalue.corr <- apply(pvalue.matrix, 1, function(x)p.adjust(x, method="BH"))
sig.p <- apply(pvalue.corr, 1, function(x)sum(x<0.001))
pvalue.corr <- cbind(pvalue.corr, sig.p)
##### for each feature, select the maximum pvalue out of ten classifiers as the final pvalue
fin.p <- apply(pvalue.corr, 1, function(x){p=x[1:10]; max(p)})
pvalue.corr <- cbind(pvalue.corr, fin.p)
imp.corr.matrix <- t(imp.corr.matrix)
write.table(pvalue.corr,"pvalues_10Runs_RF.txt")
write.table(imp.corr.matrix, "importance_score_10Runs_RF.txt")

##### for each feature, the final importance is the average among importance scores derived from 10 classifiers
#####					the final pvalue is the maximum pvalue out of ten classifiers 
imp <- apply(imp.corr.matrix, 1, mean)
pvalues <- fin.p

##### an importance rank was made based on the final importance score
pos <- data.frame(imp, pvalues)
pos <- pos[order(pos$imp, decreasing=T),]

rank <- rownames(pos)

rownames(imp.corr.matrix) <- rownames(pvalue.corr) <- gsub("X","", rownames(imp.corr.matrix))

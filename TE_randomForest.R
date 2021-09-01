require(ranger)
require(dynamicTreeCut)

data <- read.csv("pairwise comparison of TE sequences encoded with labels.csv", row.names=1)
substrate <- read.csv("TE substrate specificities.csv", row.names=1)
fa <- read.csv("TE fatty acid profiles.csv", row.names=1)

### Hierachical clustering of TEs
d <- as.matrix(dist(scale(fa)))
hc <- hclust(d=dist(scale(fa)), method="ward.D"); den.h <- as.dendrogram(hc)
cls.h <- cutreeDynamic(hc, method="hybrid",distM=d, cutHeight=50, minClusterSize=10)
names(cls.h) <- hc$label
print(table(cls.h))

c1 <- cls.h["rTE28"]
c2 <- cls.h["TEGm162"]
c3 <- cls.h["CpFatB1"]
print(c(c1, c2, c3))

cls.h[cls.h==c1] <- 33
cls.h[cls.h==c2] <- 22
cls.h[cls.h==c3] <- 11
cls.h[cls.h==33] <- 3
cls.h[cls.h==22] <- 2
cls.h[cls.h==11] <- 1

### new data label: binary
label = rep(1, nrow(data)); names(label)=rownames(data)
for(i in 1:length(label)){
	nms <- strsplit(names(label)[i], "->")[[1]]
	c <- cls.h[nms]
	if(c[1]!=c[2]) label[i]=0
}


### random forest with new data label
data1 <- data; data1$label = label

imp.matrix <- matrix(0, 10, 351, dimnames=list(paste0("RF",1:10), colnames(data)[1:351]))
imp.corr.matrix <- matrix(0, 10, 351, dimnames=list(paste0("RF",1:10), colnames(data)[1:351]))
pvalue.matrix <- matrix(0, 10, 351, dimnames=list(paste0("RF",1:10), colnames(data)[1:351]))
pred.err <- matrix(0, 2, 10, dimnames=list(c("Original", "Transformed"), paste0("Run",1:10)))

pb <- txtProgressBar(min=1, max=10, style=3)
for(i in 1:10){
	rf0 <- ranger(label~., data=data, importance="impurity_corrected", classification=T)
	rf <- ranger(label~., data=data1, importance="impurity_corrected", classification=T)
	imp.corr.matrix[i,] <- importance(rf)
	pvalue.matrix[i,] <- importance_pvalues(rf)[,2]
	pred.err[1,i] <- rf0$prediction.error
	pred.err[2,i] <- rf$prediction.error
	setTxtProgressBar(pb, i)
}

pvalue.matrix[pvalue.matrix==0] <- 1e-167
pvalue.corr <- apply(pvalue.matrix, 1, function(x)p.adjust(x, method="BH"))
sig.p <- apply(pvalue.corr, 1, function(x)sum(x<0.001))
pvalue.corr <- cbind(pvalue.corr, sig.p)
fin.p <- apply(pvalue.corr, 1, function(x){p=x[1:10]; max(p)})
pvalue.corr <- cbind(pvalue.corr, fin.p)
imp.corr.matrix <- t(imp.corr.matrix)

imp <- apply(imp.corr.matrix, 1, mean)
pvalues <- fin.p

pos <- data.frame(imp, pvalues)
pos <- pos[order(pos$imp, decreasing=T),]

rank <- rownames(pos)

rownames(imp.corr.matrix) <- rownames(pvalue.corr) <- gsub("X","", rownames(imp.corr.matrix))


################## Feature selection: incremental feature selection approach ##################
################## Necessary inputs: Substrate specificities
################## Cross-validation strategy: leave-one-out
data <- data1
substrate <- read.csv("TE substrate specificities.csv", row.names=1)

require(doSNOW)
cross=function(data, verbose=1){
	TEs <- rownames(substrate)
	MCC <- TPR <- TNR <- rep(0, length(TEs))
	names(MCC) <- names(TPR) <-names(TNR) <-TEs
    
	i <- 1
	pv <- round(length(TEs)/10)
	for(te in TEs){
		
		idx <- which(grepl(te, rownames(data)))
		train <- data[-idx, ]
		test <- data[idx, ]
	    
	    	xm <- colnames(data)
	    	xm <- xm[xm!="label"]
   	 	xnew <- test[, xm]
    		ynew <- test[, "label"]
    
		mod <- ranger(label ~., data=train, classification=T)
		yhat <- predict(mod, data=xnew)$predictions
		names(ynew) <- names(yhat) <- rownames(test)
	
		p0 <- names(ynew)[ynew == 1]
		n0 <- names(ynew)[ynew == 0]
		p1 <- "none"; if(sum(yhat==1)>0) p1 <- names(yhat)[yhat == 1]
		n1 <- "none"; if(sum(yhat==0)>0) n1 <- names(yhat)[yhat == 0]
		
		tp <- sum(p1 %in% p0)
		fp <- sum(!p1 %in% p0)
		tn <- sum(n1 %in% n0)
	   	 fn <- sum(!n1 %in% n0)
		
		TPR[te] <- tp / length(p0)
		TNR[te] <- tn / length(n0)
		mcc1 <- tp*tn - fp*fn
		mcc2 <- (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)
		MCC[te] <- mcc1 / sqrt(mcc2)

		if(verbose==1){
			if(i ==1) cat("start: ") else{
			if(i %% pv == 0) cat("->", i/pv*10,"%")
			if(i == length(TEs)) cat("->100%|", "\n")
			}
		}		
		
		i <- i+1
	}
	return(rbind(TPR, TNR, MCC))
}

ifs <- function(data, rank, start=20, nThreads=3){
	data0 <- data
		
	cl <- makeCluster(nThreads)
	registerDoSNOW(cl)
	pb <- txtProgressBar(max=length(rank), style=3)
	progress <- function(k) setTxtProgressBar(pb, k)
	opts <- list(progress=progress)
	
	result <- foreach(k=start:length(rank), .combine=cbind, .options.snow=opts, .packages="ranger",.export="cross") %dopar%{
		substrate <- read.csv("TE substrate specificities.csv", row.names=1)
		data1 <- data0[,rank[1:k]]
		data1 <- cbind(data1, label=data0[,"label"])
		eval <- cross(data=data1, verbose=0)
		tpr <- mean(eval[1,])
		tnr <- mean(eval[2,])
		mcc <- mean(eval[3,])
		return(c(tpr=tpr, tnr=tnr, mcc=mcc))
	}
    close(pb)
    stopCluster(cl)
    colnames(result) <- paste0("#Features:", start:length(rank))
	return(result)
}


eval <- list()

for(i in 1:20){
	cat("IFS:Run",i,"\n")
	eval_i <- ifs(data, rank, 2, nThreads=16)
	eval[[i]] <- eval_i
}

readme <- "This rda file contains the results of 20 IFS: a list of eval"
save(eval, readme, file="IFS_10runs_hclust.rda")

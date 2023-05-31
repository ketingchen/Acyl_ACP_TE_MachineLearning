library(ranger)
library(doSNOW)
data <- read.csv("pairwise_comparison_of_TE_sequences_encoded_with_labels.csv", row.names=1)
substrate <- read.csv("TE_substrate_specificities.csv", row.names=1)

rank <- rownames(read.table("importance_rank.txt"))

################## Feature selection: incremental feature selection approach ##################
################## Necessary inputs: Substrate specificities
################## Cross-validation strategy: ten-fold cross validation

##### change it to 10-fold cross validation
cross <- function(data, verbose=1){
        #TEs <- rownames(substrate)
        MCC <- TPR <- TNR <- rep(0, 10)
        names(MCC) <- names(TPR) <-names(TNR) <- paste("fold", 1:10, sep="_")

        #i <- 1
        #pv <- round(length(TEs)/10)
        ##### create the dataset for ten-fold cross validation
        x <- array(1:nrow(data))
        sp <- split(sample(x), 1:10)
        for(i in 1:length(sp)){
                ### crosss-validation
                ######### prepare training and testing dataset
                idx <- sp[[i]]
                train <- data[-idx, ]
                test <- data[idx, ]

           	xm <- colnames(data)
            	xm <- xm[xm!="label"]
                xnew <- test[, xm]
        	ynew <- test[, "label"]

        	######### model construction and testing 
                mod <- ranger(label ~., data=train, classification=T)
                yhat <- predict(mod, data=xnew)$predictions
                names(ynew) <- names(yhat) <- rownames(test)

                ######### evaluation of the classification performance
                ######### original data (true labels)
                p0 <- names(ynew)[ynew == 1]
                n0 <- names(ynew)[ynew == 0]
                ######### predicated labels by the RF classifier
                p1 <- "none"; if(sum(yhat==1)>0) p1 <- names(yhat)[yhat == 1]
                n1 <- "none"; if(sum(yhat==0)>0) n1 <- names(yhat)[yhat == 0]

                ######### calculation of true & false positives, and true & false negatives
                tp <- sum(p1 %in% p0)
                fp <- sum(!p1 %in% p0)
                tn <- sum(n1 %in% n0)
            	fn <- sum(!n1 %in% n0)

                ######### calculation of recall (TPR), precision (TNR), and MCC
                TPR[i] <- tp / length(p0)
                TNR[i] <- tn / length(n0)
                mcc1 <- tp*tn - fp*fn
                mcc2 <- (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)
                MCC[i] <- mcc1 / sqrt(mcc2)
		
		if(verbose==1){
                        if(i ==1) cat("start: 10%") else {
                        if(i <10 ) cat("->", i*10,"%")
                        if(i == 10) cat("->100%|", "\n")
                        }
                }
        }
        return(rbind(TPR, TNR, MCC))
}

###### Incremental feature selection strategy to identify the classifier with the minimum number of features but having an optimal performance
ifs <- function(data, rank, start=20, nThreads=3){
        data0 <- data

        cl <- makeCluster(nThreads)
        registerDoSNOW(cl)
        pb <- txtProgressBar(max=length(rank), style=3)
        progress <- function(k) setTxtProgressBar(pb, k)
        opts <- list(progress=progress)

        result <- foreach(k=start:length(rank), .combine=cbind, .options.snow=opts, .packages="ranger",.export="cross") %dopar%{
                substrate <- read.csv("TE_substrate_specificities.csv", row.names=1)
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
##### 20 runs of IFS
for(i in 1:20){
        cat("IFS:Run",i,"\n")
        eval_i <- ifs(data, rank, 2, nThreads=16)
        eval[[i]] <- eval_i
}

saveRDS(eval, "IFS_evaluation.rds")



setwd("~/Box/SILKY/Keting/Projects/Fuyuan")
eval <- readRDS("IFS_evaluation.rds")

eval0 <- matrix(0, nrow(eval[[1]]), ncol(eval[[1]]))
tpr <- tnr <- mcc <- NULL
for(i in 1:length(eval)){
	eval0 <- eval0 + eval[[i]]
	tpr <- rbind(tpr, eval[[i]]["tpr", ])
	tnr <- rbind(tnr, eval[[i]]["tnr", ])
	mcc <- rbind(mcc, eval[[i]]["mcc", ])
}
eval0 <- eval0/length(eval)

tpr.sd <- apply(tpr, 2, sd)
tnr.sd <- apply(tnr, 2, sd)
mcc.sd <- apply(mcc, 2, sd)

eval <- eval0

#########  z-transformation of correlation coefficients #########
z.trans = function(r){
	0.5*(log(1+r) - log(1-r))
}
logit.trans = function(p){
	log(p/(1-p))
}

tpr.logit = logit.trans(tpr)
tnr.logit = logit.trans(tnr)
mcc.z = z.trans(mcc)


pval <- NULL
for(i in 1:(ncol(mcc.z)-1)){
	a <- mcc.z[,i]
	b <- mcc.z[,i+1]
	mean.a <- mean(a)
	mean.b <- mean(b)
	delta <- mean.a - mean.b
	test <- t.test(a, b, paired=T)
	conf <- paste0("(",round(test$conf.int[1], 3), ",", round(test$conf.int[2],3),")")
	p <- test$p.value
	comparison <- paste0(i+1, "vs", i+2)
	##### use the non-transformed data for the table
	mean.a <- mean(mcc[,i])
	mean.b <- mean(mcc[,i+1])
	delta <- mean.a - mean.b
	pval <- rbind(pval, data.frame(comparison, mean.a, mean.b, delta, p))
}
pval$padj <- p.adjust(pval$p, method="BH")
source("/Users/kchen/Box/SILKY/Keting/Projects/Essential scripts/asterisk.R")
pval$ast <- asterisk(pval$padj)
pval.mcc <- pval
write.csv(pval.mcc, "Revision/IFS_model_comparison_mcc.csv", row.names=F)

################## Plotting based on importance score and p-values ##################
source("/Users/kchen/Box/SILKY/Keting/Projects/Essential scripts/errbar.R")

pdf("Revision/Fig6.pdf", width=3.46, height=5.11, pointsize=7)
par(mfrow=c(3,1), mar=c(4, 4, 1.5, 1), cex=1)
##########################################################################################
#### a: residue position with different levels of importance for substrate specificity
imp <- read.table("importance_score_10Runs_RF.txt")
imp <- apply(imp, 1, mean)
pval <- read.table("pvalues_10Runs_RF.txt")
pval <- pval$fin.p

col <- rep("black", length(imp)) #### non-signficant positions
col[pval<0.01] <- "dodgerblue" #### signficiant positions
col[imp>=sort(imp, decreasing=T)[26]] <- "tomato" #### most important positions as selected by IFS apporach

plot(x=1:351, y=imp, ylim=c(-2, 21), xlim=c(0, 350), axes=F, type="h", col=col, xlab="Residue position", ylab="Importance score"); box()
points(x=which(col=="tomato"), y=imp[which(col=="tomato")], pch=16, cex=0.7, col="tomato")
axis(2)
axis(1)

legend("topleft", legend=expression(bold("a")), xjust=0, bty="n")
legend(x=20, y=20, fill=c("black", "dodgerblue", "tomato"), 
		legend=c("Non-significant", "Signficant", "Essential"),
		horiz=T,
		bty="n")
##########################################################################################
#### b: model performance by IFS strategy
plot(2:351, eval[1,], type="l", ylim=c(0.5, 1.05), xlim=c(0, 350), xlab="Residue ranking", ylab="Predictive performance",col="magenta", lwd=1, axes=F);box()
axis(1)
axis(2, at=seq(0.4, 1, 0.2))
errbar(2:351, eval[1,], col="hotpink", bar=tpr.sd)
points(2:351, eval[2,], type="l", col="coral", lwd=1)
errbar(2:351, eval[2,], col="indianred2", bar=tnr.sd)
points(2:351, eval[3,], type="l", col="forestgreen", lwd=1)
errbar(2:351, eval[3,], col="green4", bar=mcc.sd)

legend("topleft", legend=expression(bold("b")), xjust=0, bty="n")
legend(x=20, y=1.05, fill=c("coral", "hotpink", "forestgreen"), legend=c("Recall","Specificity", "MCC"), horiz=T, bty="n")
##########################################################################################
#### c: model performance by IFS strategy (a zoom in view that only include models with no more than 30 features)
plot(2:30, eval[1,2:30], type="l", ylim=c(0.4, 1.05), xlab="Residue ranking", ylab="Predictive performance",col="magenta", lwd=1, axes=F); box()
axis(1)
axis(2, at=seq(0.4, 1, 0.2))
errbar(2:30, eval[1,2:30], col="hotpink", bar=tpr.sd[2:30])
points(2:30, eval[2,2:30], type="l", col="coral", lwd=1)
errbar(2:30, eval[2,2:30], col="indianred2", bar=tnr.sd[2:30])
points(2:30, eval[3,2:30], type="l", col="forestgreen", lwd=1)
errbar(2:30, eval[3,2:30], col="green4", bar=mcc.sd[2:30])
points(27, eval[1, 27], pch=19, col="magenta", cex=1)
points(27, eval[2, 27], pch=19, col="coral", cex=1)
points(27, eval[3, 27], pch=19, col="forestgreen", cex=1)

legend("topleft", legend=expression(bold("c")), xjust=0, bty="n")

dev.off()
##########################################################################################
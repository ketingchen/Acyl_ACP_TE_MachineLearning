pca.nipals = function(X, thresh=1e-6, lv = 2, continue = T, scale=F)
{
  X = scale(as.matrix(X), center=T,scale=scale)
  scores = loadings = NULL
  variance = rep(0, lv)
  TotalSS = sum(X^2)
  # For each component
    for(i in 1:lv)
    {
        t = X[,sample(1:ncol(X), 1)] #pick a random starting point for t
        continue = T
        #For each iteration to search the optimal component
        while(continue)
        {
          t_u = t/sqrt(sum(t^2)) #scale t to be of unit size
          p_t = crossprod(t_u, X) #p = t(T) %*% X/(t(T)%*%T)
          p_t = p_t/sqrt(sum(p_t^2)) #scale p to be of unit size
          t.old = t
          t = tcrossprod(X, p_t) #T = XP/(t(P)%*%P)
          if(crossprod(t.old-t) <= thresh) continue = F
        }
        X = X - t%*%p_t
        scores =cbind(scores, t)
        loadings = rbind(loadings, p_t)
        variance[i] = 1 - sum(X^2)/TotalSS
        if(i>1) variance[i] = variance[i]-sum(variance[1:i-1])
    }
  colnames(scores)=paste("PC",1:lv,sep="")
  colnames(loadings)=colnames(X); rownames(loadings)=paste("PC",1:lv, sep="")
  names(variance)=paste("PC",1:lv,sep="")
    
  return(list("Scores"=scores, "Loadings"=loadings, "Explained"=variance))
}

##### PCA by SVD #####
pca.svd=function(X, lv=2, scale=T){
  if(nrow(X)==1) return("Sample size has to be >1")
  if(lv>ncol(X)) lv=ncol(X)
  #scale the data
  Xs=scale(X, center=T, scale=scale)
  # no need to generate the left singular vectors at this moment, sigmas and right vectors (loadings) are more important
  Xs.svd=svd(Xs, nu=0, nv=lv)
  
  sigma2=(Xs.svd$d/sqrt(nrow(X)-1))^2
  pct=sigma2[1:lv]/sum(sigma2)
  cum=rep(0, length(pct));for(i in 1:length(pct)) cum[i]=ifelse(i==1, pct[i], cum[i-1]+pct[i])
  
  explained=rbind(pct, cum)
  
  loadings=Xs.svd$v[,1:lv, drop=F]
  
  # Estimate scores based on Xs and loadings
  scores=Xs %*% loadings
  
  colnames(scores)=colnames(loadings)=colnames(explained)=paste0("PC",1:lv)
  
  return(invisible(list(Scores=scores, Loadings=loadings, Explained=explained)))
}


##### NIPALS algorithm to build a PLS model
pls.nipals = function(X, Y, thresh=1e-06, lv=2, continue=T, scale=F,center=T)
{
  #center X and Y
  X0 = scale(X, center=center, scale=scale)
  Y0 = scale(Y, center=center, scale=scale)
  
  X0[is.na(X0)]=0; X0.1=X0
  Y0[is.na(Y0)]=0; Y0.1=Y0
  
  SSX0 = sum(X0^2); SSY0 = sum(Y0^2)
  R2X = RSS = R2Y = NULL
  #Initialize returned values
  TT = UU = PP = QQ = WW = BB = NULL
  ############ Each component ############
  for (i in 1:lv)
  {
    #Initialize u
    u = Y0[,sample(1:ncol(Y0), 1),drop=F]
  
    continue = T
    t_old = matrix(0, nrow(X0), 1)
  ############ Each iteration within a component: first determine w_t, t, q_t, u
      while(continue)
      {
        #X block
        w = crossprod(X0, u)/drop(crossprod(u))
        w = w/sqrt(drop(crossprod(w))) #Normalize
        t = X0 %*% w
  
        #Y block
        q = crossprod(Y0, t)/drop(crossprod(t))
        q = q/sqrt(drop(crossprod(q))) #Normalize
        u = Y0 %*% q
    
        #check convergence
        if(crossprod(t_old-t)>thresh) t_old = t
        else continue = F
      }
  ############ After convergence, determine the returned values for the component
    p = crossprod(X0, t)/drop(crossprod(t))
    
    TT = cbind(TT, t)
    PP = cbind(PP, p)
    UU = cbind(UU, u)
    QQ = cbind(QQ, q)
    WW = cbind(WW, w)
      
    #Calculate the coefficient for inner relation between u and t
    b = drop(crossprod(u, t)/drop(crossprod(t)))
    BB = c(BB, b)  
    #Calculate the residuals
    X0 = X0 - t%*%t(p)
    Y0 = Y0 - b*t%*%t(q)
    R2X = c(R2X, 1 - sum(X0^2)/SSX0)
    R2Y = c(R2Y, 1 - sum(Y0^2)/SSY0)
    RSS = c(RSS, sum(Y0^2))
  }
  RR = WW %*% solve(t(PP) %*% WW)
  R2X = diff(c(0, R2X))
  R2Y = diff(c(0, R2Y))
  
  
  colnames(TT)=colnames(UU)=colnames(PP)=colnames(QQ)=colnames(WW)=colnames(RR)=names(R2X)=names(R2Y)=paste("LV",1:lv,sep="")
  return(list("X0"=X0.1, "Y0"=Y0.1, "LatentNum"=lv,"ScoresX"=TT, "ScoresY"=UU, "LoadingsX"=PP, "LoadingsY"=QQ,"InnerXY"=BB,"Rotation" = RR, "Weights"=WW, "R2X" = R2X, "R2Y" = R2Y, "RSS"=RSS))
}

##### SIMPLSalgorithm to build a PLS model
pls.sim = function(X, Y, lv=1, scale=F, center=T)
{
  X0 = scale(X, center=center, scale=scale)
  Y0 = scale(Y, center=center, scale=scale)
  
  X0[is.na(X0)]=0
  Y0[is.na(Y0)]=0
  # Calculate the total sum of squares for X and Y
    SSX0 = sum(X0^2)
    SSY0 = sum(Y0^2)
  
  n = dim(X0)[1]
  RR=TT=PP=QQ=UU=VV=RSS=NULL
  S= t(X0)%*%Y0
  for(i in 1:lv)
  {
    qq = svd(S)$v[,1] #qq is the eigenvector of S
    rr = S%*%qq  #Rotation vector for X
    tt = X0%*%rr #X block factor scores
    tt = scale(tt, center=T, scale=F)
    tnorm = as.numeric(sqrt(t(tt)%*%tt))
    tt = tt/tnorm  # Normalize tt
    rr = rr/tnorm  # Normalize rr
    pp = crossprod(X0, tt) # X blocking factor loadings
    qq = crossprod(Y0, tt) # Y blocking factor loadings
    uu = Y0 %*% qq # Y blocking factor scores
    vv = pp #Initiate orthogonal loadings
    if (i > 1)
    {
      vv = vv - VV %*%crossprod(VV, pp)
      uu = uu - TT %*%crossprod(TT, uu)
    }
    vv = vv/sqrt(drop(crossprod(vv))) ##vv is the orthogonal loadings of X, it is the same output from NIPALS algorithm
    S = S - vv %*% crossprod(vv, S)
    RR = cbind(RR, rr)
    TT = cbind(TT, tt)
    PP = cbind(PP, pp)
    QQ = cbind(QQ, qq)
    UU = cbind(UU, uu)
    VV = cbind(VV, vv)
    
    Y1 = X0 %*% RR %*% t(QQ)
    RSS = c(RSS, sum((Y0 - Y1)^2))
  }
  
  
   R2X = diag(crossprod(PP))/SSX0
   R2Y = diag(crossprod(QQ))/SSY0
 
  colnames(RR) = colnames(TT) = colnames(PP) = colnames(QQ) = colnames(UU) = colnames(VV) = names(R2X) = names(R2Y) = paste(rep("LV", lv), 1:lv, sep="")
  return(list("X0"=X0, "Y0"=Y0, "LatentNum"=lv, "Weights"=RR, "ScoresX"=TT, "LoadingsX"=PP, "ScoresY"=UU, "LoadingsY"=QQ, "OrthogonalX"=VV, "R2X"=R2X, "R2Y"=R2Y, "RSS" = RSS))
}


# Calculate the relationship between X and Y for each component
RtoQ = function(pls, lv=1)
{
	RR = pls$Rotation
	QQ = pls$LoadingsY
	sapply(QQ[,lv], function(x) RR[,lv]/x)
}

 
# calcuate vip for each variable, only for NIPALS at this stage
# pls: pls object created by pls.nipals
vips = function(pls)
{
  vip = NULL
  p = dim(pls$Weights)[1]
  R2Y = pls$R2Y
  w2 = pls$Weights^2
  if(dim(w2)[2]>1)
  	vips = apply(w2 %*% diag(R2Y), 1, sum)/sum(R2Y)
  else
    {vips = c(w2); names(vips) = rownames(w2)} 
  sqrt(p)*sqrt(vips)
}

#### Permutation tests to assess VIP ######
permute = function(X, Y, lv=2, scale=F, method=c("nipals", "sim"), run=1000,report=100,adjust=F)
{
  num = dim(X)[1]
  Y1 = Y
  if(method[1] == "nipals")
  	pls0 = pls.nipals(X, Y, lv=lv, scale=scale)
  else
    pls0 = pls.sim(X, Y, lv=lv, scale=scale)	
  
  ss0 = vips(pls0)*sum(pls0$R2Y)
    
  pval = rep(1, dim(X)[2])
  names(pval) = colnames(X) 
  pb = txtProgressBar(min=0, max=run, style=3)
  ##Create the random population of weight for each sequence by permutation tests
  for(i in 1:run)
       {
          X1 = X[sample(1:num),]
          if(method[1] == "nipals")
          	pls1 = pls.nipals(X1, Y1, lv=lv, scale=scale)
          else
            pls1 = pls.sim(X1, Y1, lv=lv, scale=scale)
            
          ss1 = vips(pls1)*sum(pls1$R2Y)
          
          pval = pval+as.numeric(ss1 >= ss0)
          setTxtProgressBar(pb, i)
  }
  cat("\n")
  pval=pval/(run+1) 
  if(adjust) pval = p.adjust(pval, method ="BH")
  pval
}


### Calculate the portions of response variance explained by selected transcripts
### pls: pls object created by pls.nipals
explain.SSy = function(pls, seqs)
{
  p = dim(pls$Weights)[1]
  lv = pls$LatentNum
  w = vip(pls)/sqrt(p)
  if (length(seqs)==1) w[seqs]^2*pls$R2Y[lv]
  else sum(w[seqs]^2)*pls$R2Y[lv]
}

#################### cross validation for pls object####################
# pls: an object returned by pls.nipals or pls.sim
# k: number of folds in cross-validation
# progress: a logical value of whether to show the progress bar
# fig: a logical value of whether to plot the cumulative Q2
Q2Y.cv = function(pls, k=10, thresh=1e-6, progress=T, fig=F)
{
  X=pls$X0; Y=pls$Y0
  n = dim(X)[1]
  vmax = pls$LatentNum
  
  PRESS.inside = Q2 = matrix(0, nrow = vmax, ncol = dim(Y)[2])
  press.mat = lapply(1:vmax, function(x){matrix(NA, nrow = n, ncol = dim(Y)[2])})
  RSS.indiv = lapply(1:(vmax+1), function(x){matrix(NA, nrow = n, ncol = dim(Y)[2])})
  RSS.indiv[[1]]=Y
  RSS = matrix(0, nrow = vmax+1, ncol = dim(Y)[2])
  RSS[1,] = colSums(Y^2)
  folds = createfolds(n, k)
  
  if(progress){pb = txtProgressBar(style = 3); nBar=1}
  for (v in 1:vmax)
  {   
    tt = pls$ScoresX[,v]
    u=pls$ScoresY[,v]
    p=crossprod(X, tt)/drop(crossprod(tt))
    q=crossprod(Y, tt)/drop(crossprod(tt))
    
    RSS.indiv[[v+1]]=Y - tt%*%t(q)
    RSS[v+1,]= colSums((Y - tt%*%t(q))^2)
    
    for (i in 1:k)
    {
      if(progress)
      {setTxtProgressBar(pb, nBar/(vmax*k)); nBar=nBar+1}
      
      omit=folds[i,]
      X.train = X[-omit, , drop = FALSE]
      Y.train = Y[-omit, , drop = FALSE]
      X.test = X[omit, , drop = FALSE]
      Y.test = Y[omit, , drop = FALSE]
      u.cv = u[-omit]
      w.old.cv=0
      repeat{
        w.cv = crossprod(X.train, u.cv)/drop(crossprod(u.cv))
        w.cv = w.cv/sqrt(drop(crossprod(w.cv)))
        t.cv = X.train %*% w.cv
        q.cv = crossprod(Y.train, t.cv)/drop(crossprod(t.cv))
        q.cv = q.cv/sqrt(drop(crossprod(q.cv)))
        u.cv = Y.train %*% q.cv
        if(crossprod(w.cv-w.old.cv)<thresh) break
        w.old.cv = w.cv
      }
      t.cv = X.train %*% w.cv
      q.cv = crossprod(Y.train, t.cv)/drop(crossprod(t.cv))
      Y.hat.cv = X.test %*% w.cv %*% t(q.cv)
      press.mat[[v]][omit,]=Y.test-Y.hat.cv
    } # end k-fold cross-validation for component c
    
    PRESS.inside[v,]=apply(press.mat[[v]], 2, crossprod)
    Q2[v,] = 1 - PRESS.inside[v,]/RSS[v,]
    
    #Data deflation: prepare for next component
    X = X - tt %*% t(p)
    Y = Y - tt %*% t(q)
  } # end loop on c of 1:cmax
  if(progress) cat("\n")
  Q2.lv = matrix(1 - rowSums(PRESS.inside)/rowSums(RSS[-(vmax+1),,drop=F]), nrow=1, ncol=vmax,
                 dimnames = list("Q2.lv", 1:vmax))
  rownames(Q2) = paste("LV", 1:vmax, sep="")
  colnames(Q2) = colnames(Y)
  
  Q2.lv = t(Q2.lv)
  Q2.cum = NULL
  for(i in 1:vmax)
  {
    Q2.cum = c(Q2.cum, 1-prod(1-Q2.lv[1:i]))
  }
  names(Q2.cum)=1:vmax
  if(fig) plot(Q2.cum, type="b", xlab="Number of latent variables")
  
  
  cbind(Q2.lv, Q2.cum)
}

########################  Permutation on Q2Y to avoid overfitting ######################## 
Q2Y.perm = function(pls0, run=1000, method = c("NIPALS", "sim"))
{
  X0 = pls0$X0; Y0 = pls0$Y0; lv = pls0$LatentNum
  n = dim(X0)[1]
  Q2Y0 = Q2Y.cv(pls0, progress=F, fig=F)[lv,"Q2.cum"]
  pvalue = 1
  pb = txtProgressBar(min=0, max=run, style = 3)
  
  for (i in 1:run)
  {
    X1 = X0[sample(1:n),,drop=F]; Y1 = Y0[sample(1:n),,drop=F]
    if(method[1]=="NIPALS")
      pls1 = pls.nipals(X=X1, Y=Y1, lv=lv, scale=F, center=F)
    else
      pls1 = pls.sim(X=X1, Y=Y1, lv=lv, scale=F, center=F)
    Q2Y1 = Q2Y.cv(pls1, progress=F, fig=F)[lv,"Q2.cum"]

    pvalue = pvalue + as.numeric(Q2Y1>=Q2Y0)
    setTxtProgressBar(pb, i)
  }
  cat("\n")
  pvalue / (run+1)
}

######  This is a function that compute mean square error of prediction under assigned number of components
###### prev: pls model with c-1 latent variables
###### X: explanatory variables
###### Y: response variables
###### folds: randomization scheme to divide X and Y into K-folds. It is returned from createfolds
###### method: nipals or sim

R2Y.perm = function(X, Y, scale=F, lv=1, run=1000, method=c("nipals", "sim"))
{
	p= 1
	n=dim(X)[1]
	
	if(method[1] == "nipals")
		  R2Y0 = sum(pls.nipals(X=X, Y=Y, scale=scale, lv=lv)$R2Y)
	else
	   	R2Y0 = sum(pls.sim(X=X, Y=Y, scale=scale, lv=lv)$R2Y)
	
	pb = txtProgressBar(min=0, max=run, style = 3)
	
   for(i in 1:run)
   {

   	  X1 = X[sample(1:n),,drop=F]
	    Y1 = Y[sample(1:n),,drop=F]
	  
	  if(method[1] == "nipals")
		  R2Y1 = sum(pls.nipals(X=X1, Y=Y1, scale=scale, lv=lv)$R2Y)
	  else
	   	R2Y1 = sum(pls.sim(X=X1, Y=Y1, scale=scale, lv=lv)$R2Y)

  	  p = p + as.numeric(R2Y1>=R2Y0)
  	  setTxtProgressBar(pb, i)
   }
   p/(1+run)
}


######## This is a function that randomized data into different folds for cross validation purpose ######
###### The returned value will be a k-row matrix ######
###### Numbers in the matrix correspondings to a sample row in a nxp or nxm dataset
###### Each row represent a fold which contains a random set of samples #####
###### Note that the returned contained zero values, which should be removed in the cross-validation script
createfolds = function(n,k)
{
  f = floor(n/k)
  folds = matrix(rep(0, k*ceiling(n/k)), k, ceiling(n/k))
  total = 1:n
  for(i in 1:k)
  {
    if(i < k)
    {
      folds[i,1:f]= sample(total, f)
      total = setdiff(total, folds[i,])
    }
    else 
    {
     #### Till this step, each fold should contain floor(n/k) samples
     #### The last fold, fold K should contain at most ceiling(n/k) samples
     #### However, in some cases, more samples were left for fold K
     #### Therefore, this step allows the extra samples to be randomly assigned to other folds
      d = length(total)-ceiling(n/k)
      if(d > 0) 
      {
        extra = sample(total, d)
        folds[sample(1:(k-1), d), f+1] = extra
        folds[i, ] = setdiff(total, extra)
      }
      else
        folds[i, ] = total
    }
  }
  return(folds)
}



######### This function is to determine a cut-off value for gene selection #######
######## nPC: the number of latent variables ######## 
######## permute----- half: only permute X (permute = 0) or Y (permute = 1) ######## 
########              full: permute both X and Y (permute = 2) ######## 
######## alpha: level of signficance. This script uses alpha to find out the cut off values
########        It means (1- alpha/2) * 100% of genes in the random population is below the positive cut off value
########        And (1 - alpha/2) *100% of genes in the random population is above the negative cut off value
######## side: one-side or two-sided test
########       0: right side, (1-alpha)% percentile
########       1: left side, alpha% percentile
########       2: two-sided, (1-alpha/2)% and (alpha/2)% percentile
cutoff = function(X, Y, lv=2, scale = F, permute=2, alpha=0.05, side=2)
{
  X0 = X
  Y0 = Y
  num = dim(X0)[1]
  
  ###### Permutation ######
  if(permute == 0)
  {
    X0 = X0[sample(1:num, num),]
    }
  else if(permute == 1)
  {
    Y0 = Y0[sample(1:num, num),]
  }
  else
  {
    X0 = X0[sample(1:num, num),]
    Y0 = Y0[sample(1:num, num),]
  }
  
  ###### Find the weight matrix for random sample######
  rand = pls.sim(X0, Y0, lv=lv, scale = scale)$Weights
  ###### Find the cut-off value based on given alpha, note that cutoff values are two-sided
  if(side==0)
    {cut = apply(rand, 2, function(x){quantile(x, 1-alpha)})
     names(cut) = paste(rep("PC", lv), 1:lv, sep = "")}
  else if (side == 1)
    {cut = apply(rand, 2, function(x){quantile(x, alpha)})
     names(cut) = paste(rep("PC", lv), 1:lv, sep = "")}
  else
    {cut =c(apply(rand, 2, function(x){quantile(x, 1-alpha/2)}), apply(rand, 2, function(x){quantile(x,alpha/2)}))
  names(cut) = c(paste(paste(rep("PC", lv), 1:lv, sep=""),"+"),paste(paste(rep("PC", lv), 1:lv, sep=""),"-"))}
  return(cut)
}


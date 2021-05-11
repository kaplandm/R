# Feedback: kaplandm@missouri.edu
# Author: David M. Kaplan
# Nonparametric confidence intervals (etc.) for conditional quantile-related objects.
# 
# quantile.inf.np() to calculate confidence intervals
# plot.quantile.inf.np() to plot results (separate graphs for discrete conditioning variables)
# plot2.quantile.inf.np() to plot results all on same graph
# other functions herein are for internal use only
# 
# References (chronological):
# 'Fractional order statistic approximation for nonparametric conditional quantile inference' by Matt Goldman and David M. Kaplan, 2017, https://doi.org/10.1016/j.jeconom.2016.09.015
# 'Non-parametric inference on conditional quantile differences and linear combinations, using L-statistics' by Matt Goldman and David M. Kaplan, 2018, https://doi.org/10.1111/ectj.12095
# 
# See also https://github.com/kaplandm/R/blob/main/quantile_inf.R
# for related *unconditional* inference.
# 
# Empirical examples: https://github.com/kaplandm/R/blob/main/quantile_inf_np_examples.R
# Simulations: https://github.com/kaplandm/R/blob/main/quantile_inf_np_sims.R


###################
#                 #
# Load packages   #
#                 #
###################
prob.loaded <- (exists("quantile.inf") && exists("quantile.inf.CIuh"))
success <-
  tryCatch({source('https://raw.githubusercontent.com/kaplandm/R/main/quantile_inf.R');TRUE},
           error=function(w) FALSE)
if (!success) {
  success <-
    tryCatch({source("quantile_inf.R");TRUE},
             error=function(w) FALSE)
  if (!success) {
    if (prob.loaded) {
      warning("Couldn't load quantile_inf.R, but it seems like you already did.")
    } else {
      stop("Failed to source() quantile_inf.R from web or local file.  You may download and source() it yourself, or at least make sure it's in your getwd().  Currently available at https://github.com/kaplandm/R")
    }
  }
}
#
tryCatch(library(ks), error=function(w) {cat("Must install package ks\n");stop(w)}) #kde and kdde
tryCatch(library(quantreg), error=function(w) {cat("Must install package quantreg\n");stop(w)}) #for rq
#
#NOW OPTIONAL:
# library(splines)
# library(fields) #qsreg
# library(np) #kernel conditional density estimation
#
# library(stats) #auto-loads


###################
#                 #
# Constants       #
#                 #
###################
HCOR <- 1 #coefficient for finite-sample correction; see below
EPS.MULT <- 1e8


###################
#                 #
# Bandwidth fns   #
#                 #
###################
# See appendix of paper
hhat.2s.fn <- function(p,n,z,Fpp,Fp,f.X,fp.X,DENOM.MIN=0) {
  ifelse(p==0.5,
         n^(-1/3) *(abs(f.X*Fpp+2*fp.X*Fp)+DENOM.MIN)^(-1/3),
         n^(-1/3) *
         ((sign(-(f.X*Fpp+2*fp.X*Fp))*(1-2*p)
           +sqrt((1-2*p)^2+4))  /
         (2*DENOM.MIN+2*abs(f.X*Fpp+2*fp.X*Fp)))^(1/3)
  ) } #2nd version is valid for p==0.5, but slower to compute
#
hhat.1s.fn <- function(p,n,z,Fpp,Fp,f.X,fp.X,DENOM.MIN=0) {
  n^(-3/7) *
  (z / (DENOM.MIN+3*sqrt(p*(1-p)*f.X)
  *ifelse(f.X*Fpp+2*fp.X*Fp>0,
          f.X*Fpp+2*fp.X*Fp,
          (f.X*Fpp+2*fp.X*Fp)*(-5/2))
  ))^(2/7)  }


###################
#                 #
# Validate args   #
#                 #
###################
# Check for valid input data
quantile.inf.np.validate.inputs <- function(Y,X,x0s,p,ALPHA,JOINT.FLAG,ONESIDED,KERNEL.TYPE,METHOD.TYPE,PSI,hs.pt,hs.joint) {
  #Check for NULL values.
  for (varname in c('Y','X','x0s','p','ALPHA','JOINT.FLAG','ONESIDED','KERNEL.TYPE','METHOD.TYPE')) {
    if (is.null(get(varname))) stop(sprintf("%s cannot be NULL.",varname))
  }
  #Y,p,PSI,hs.pt,hs.joint: variable type
  for (varname in c('Y','p','PSI','hs.pt','hs.joint')) {
    if (!is.null(get(varname))) {
      if (!(is.vector(get(varname)) || is.list(get(varname)))) stop(sprintf("%1$s must be a vector (and not a list).  If %1$s is a one-dimensional array or matrix, you can use as.vector(%1$s) in order to cast %1$s as a vector.",varname))
      if (!is.numeric(get(varname))) stop(sprintf("%1$s must be numeric; i.e., is.numeric(%1$s) must be TRUE.",varname))
      if (any(is.nan(get(varname))) || any(is.na(get(varname)))) stop(sprintf("%1$s cannot have NA or NaN values; i.e., any(is.nan(%1$s)) and any(is.na(%1$s)) must be FALSE.",varname))
    }
  }
	n <- length(Y)
	#X type, size
	if (!is.numeric(X)) stop("X must be numeric; i.e., is.numeric(X) must be TRUE.  (For non-numeric discrete components of X, please convert them to integers so that X is numeric.)")
	X.discrete <- NULL;  X.discrete.n <- 1
	if ((is.vector(X) && !is.list(X)) || (is.matrix(X) && (dim(X)[2]==1 || dim(X)[1]==1))) {
	  if (length(X)!=n) stop(sprintf("X must be the same length as Y, but length(Y)=%d and length(X)=%d.",length(Y),length(X)))
	  if (METHOD.TYPE=='qte') stop("X needs to contain at least one binary treatment indicator in addition to the one continuous regressor.")
	} else if (is.matrix(X)) {
	  tmp <- dim(X)
	  if (tmp[1]!=n) {
	    if (tmp[2]!=n) {stop(sprintf("dim(X)[1] should equal length(Y), the number of observations, but length(Y)=%d and dim(X)[1]=%d and dim(X)[2]=%d.",length(Y),dim(X)[1],dim(X)[2]))} else {X <- t(X)}
      }
      if (dim(X)[2]>dim(X)[1]) stop(sprintf("X contains %d regressors, which is more than the number of observations, %d.",dim(X)[2],dim(X)[1]))
      if (METHOD.TYPE=='qte') {
        tmp <- sort(unique(X[,2])) #treatment indicator
        if (!(length(tmp)==2 && tmp[1]==0 && tmp[2]==1)) stop(sprintf("For method qte, the second column of X must be a binary treatment indicator taking only values zero and one."))
      }
    } else stop("X must be a vector (and not a list) or a matrix.")
    #x0s type, size, sorting
    if (!is.vector(x0s) || is.list(x0s)) stop(sprintf("x0s must be a vector (and not a list).  If x0s is a one-dimensional array or matrix, you can use as.vector(x0s) in order to cast x0s as a vector."))
    if (is.unsorted(x0s)) { stop("x0s must be sorted in increasing order.") }
    #PSI, multiple p
    if (length(p)>1) {
      if (METHOD.TYPE %in% c('lincom','qte')) { 
        if (length(p)!=length(PSI)) { stop(sprintf("Arguments p and PSI must have same length, but length(p)=%d and length(PSI)=%d",length(p),length(PSI))) }
        if (!is.null(PSI) && PSI[1]!=1) { stop('Argument PSI should be normalized to have PSI[1]=1') }
      } else if (METHOD.TYPE=='joint' && !is.null(PSI) && !is.na(PSI)) { 
        warning("Ignoring argument PSI since METHOD.TYPE=='joint'");  PSI <- NULL
      }
      if (METHOD.TYPE=='single') { stop(sprintf("For METHOD.TYPE=='single', argument p must be a scalar, but length(p)=%d",length(p))) }
    } else if (!is.null(PSI) && !is.na(PSI)) { 
      warning("Ignoring argument PSI since length(p)=1") 
      if (METHOD.TYPE %in% c('lincom','joint')) { METHOD.TYPE <- 'single' }
    }
  
  #hs: same length as x0s
	if (!is.null(hs.pt)) {
	  if (is.list(hs.pt) && METHOD.TYPE=='qte') {
	    if (length(hs.pt$t)!=length(x0s) || length(hs.pt$c)!=length(x0s)) {
	      stop("Length of hs.pt$t must equal length(x0s)")
	    }
	  } else if (is.vector(hs.pt) && METHOD.TYPE!='qte') {
	    if (length(hs.pt)!=length(x0s)) { stop("Length of hs.pt must equal that of x0s.") }
	  } else {
	    stop(sprintf("hs.pt needs to be either 1) NULL (recommended), 2) a list if METHOD.TYPE is 'qte', 3) a vector if METHOD.TYPE is not 'qte'; see comments above function code for quantile.inf.np(). class(hs.pt)=%s",class(hs.pt)))
	  }
	}
	if (!is.null(hs.joint)) {
	  if (is.list(hs.joint) && METHOD.TYPE=='qte') {
	    if (length(hs.joint$t)!=length(x0s) || length(hs.joint$c)!=length(x0s)) {
	      stop("Length of hs.joint$t must equal length(x0s)")
	    }
	  } else if (is.vector(hs.joint) && METHOD.TYPE!='qte') {
	    if (length(hs.joint)!=length(x0s)) { stop("Length of hs.joint must equal that of x0s.") }
	  } else {
	    stop("hs.joint needs to be either 1) NULL (recommended), 2) a list if METHOD.TYPE is 'qte', 3) a vector if METHOD.TYPE is not 'qte'; see comments above function code for quantile.inf.np().")
	  }
	}
}


###################
#                 #
# Obj estimation  #
#                 #
###################
# LOCAL: use local (to x0s) nearest-neighbor-type estimators
# currently only supports one-dimensional X (this is *after* conditioning on any discrete X)
quantile.inf.np.est.objs <- function(Y,X,x0s,p,KERNEL.TYPE='gaussian',LOCAL=NA,kNN=FALSE) {
  #bandwidth: (normal rule-of-thumb-type)
  kde.h <- hns(X,deriv.order=0)
  f.X.hats <- kde(X,eval.points=x0s,h=kde.h)$estimate #or another estimator
  fp.X.hats <- kdde(X, h=hns(X,deriv.order=1), deriv.order=1, eval.points=x0s)$estimate
  #
  LOCALeff <- LOCAL
  if (is.na(LOCALeff)) LOCALeff <- (length(Y)>=1e3)
  if (!LOCALeff) {
    if (!(require(splines) && require(fields) && require(np))) {
      LOCALeff <- TRUE
    }
  }
  if (LOCALeff) {
    if (kNN) {
      n <- length(Y); C <- 2/5
      Nn <- min(n,max(10,round(C*n^(4/5)))) #optimal local poly estimation (incl. derivatives) for smoothness sQ=2 as assumed in paper (well, 2.00...1). happens to also be optimal for KDE (though not KDDE)
      sortinds <- order(X,na.last=NA) #don't include NA, but keep original inds
      sortXinds <- cbind(X[sortinds],sortinds)
      rm(sortinds)
      #     wheres <- apply(X=matrix(X,1),MARGIN=2,FUN=function(x)sum(x<x0s))
      minX <- sortXinds[1,1];  maxX <- sortXinds[n,1]
      tmpX <- vector("list",length(x0s));  cntr <- 1
      if (min(x0s)<minX) { stop("Error: min(x0s)<min(X)") }
      if (max(x0s)>maxX) { stop("Error: max(x0s)>max(X)") }
      for (i in 1:n) {
        if (sortXinds[i+Nn,1]>=x0s[cntr]) {
          tmpX[[cntr]] <- sortXinds[i:min(n,i+2*Nn-1),]
          cntr <- cntr+1
          if (cntr>length(x0s)) { break }
        }
      }
    } else {
      if (min(x0s)<min(X)) { stop("Error: min(x0s)<min(X)") }
      if (max(x0s)>max(X)) { stop("Error: max(x0s)>max(X)") }
    }
    Y.hats <- Fp.hats <- Fpp.hats <- rep.int(NA,length(x0s)) #f.Y.hats <- 
    for (i in 1:length(x0s)) {
      if (kNN) {
        #first, pinpoint which indices (from original X,Y) are nearest Nn to each x0
        absdevs.raw <- abs(x0s[i]-tmpX[[i]][,1])
        absdevs <- ifelse(absdevs.raw>min(abs(x0s[i]-minX),abs(x0s[i]-maxX)),NA,absdevs.raw)
        if (sum(!is.na(absdevs))<10) {
          keepnum <- 10;  absdevs <- absdevs.raw  #ok if asymmetric, but only take 10 obs to minimize bias
        } else {keepnum <- Nn}
        absinds <- sort.list(absdevs,na.last=NA,method='quick')
        keepers <- absinds[1:min(length(absinds),keepnum)]
        h <- absdevs[!is.na(absdevs)][absinds][length(keepers)] #equivalent bandwidth, for scaling
        orig.keep <- tmpX[[i]][,2][!is.na(absdevs)][keepers]
      } else {
        h <- kde.h / 1.5
        while (h<max(X)-min(X)) {
          h <- h * 1.5
          orig.keep <- (X<=x0s[i]+h) & (X>=x0s[i]-h)
          if (length(unique(X[orig.keep]))>3) {
            break
          }
        }
      }
      #Set up local sample
      noninf <- (abs(X)!=Inf) & (abs(Y)!=Inf)
      Yloc <- Y[orig.keep & noninf]
      Xloc <- X[orig.keep & noninf] - x0s[i]  #center X s.t. x0=0, as in Chaudhuri (1991)
      if (length(unique(Xloc))<=3) {
        Fp.hats[i] <- Fpp.hats[i] <- 0
        next
      }
      #Yloc <- Yloc[order(Xloc)];  Xloc <- Xloc[order(Xloc)] #FOR MANUAL INSPECTION ONLY
      Xdes2 <- Xloc/h
      Xdes3 <- Xloc^2/h^2
      Xdes4 <- Xloc^3/h^3
      #Y.hats
      done <- FALSE
      if (!done && length(Yloc)>=20) { #only use cubic if "enough" data (threshold is somewhat arbitrary)
        done <- tryCatch({
          rq.ret <- rq(Yloc~Xdes2+Xdes3+Xdes4,tau=p) #constant incl. automatically
          TRUE}, error=function(ERR) {FALSE})
      }
      if (!done && length(Yloc)>=10) {
        done <- tryCatch({
          rq.ret <- rq(Yloc~Xdes2+Xdes3,tau=p)
          TRUE}, error=function(ERR) {FALSE})
      }
      if (!done) {
        useNA <- tryCatch({
          rq.ret <- rq(Yloc~Xdes2,tau=p)
          FALSE}, error=function(ERR) {TRUE})
        if (useNA) {
          Y.hats[i] <- Fp.hats[i] <- Fpp.hats[i] <- NA
          next
        }
    }
      Y.hats[i] <- predict(rq.ret,newdata=data.frame(Xdes2=0,Xdes3=0,Xdes4=0))
#       #f.Y.hats; density() is faster and seems fine
# #       f.Y.hats[i] <- kde(Yloc,eval.points=Y.hats[i],h=hns(Yloc,deriv.order=0))$estimate
#       f.Y.hats[i] <- density(Yloc,from=Y.hats[i],to=Y.hats[i],n=1)$y
      #Fp.hats, Fpp.hats
      Ybelowloc <- as.integer(Yloc<=Y.hats[i])
      if (FALSE) { #raw polynomial
        #in general, to est. order v derivative w/ degree p polynomial, want p-v=odd, not even
#         lmquad <- lm(Ybelowloc~Xdes2+Xdes3)$coefficients
        lmquad <- lm(I(Ybelowloc-p)~Xdes2+Xdes3+0)$coefficients #constrained to fit (0,0.25)
        if (length(Ybelowloc>=20)) {
#           lmcub <- lm(Ybelowloc~Xdes2+Xdes3+Xdes4)$coefficients
          lmcub <- lm(I(Ybelowloc-p)~Xdes2+Xdes3+Xdes4+0)$coefficients #constrained to fit (0,0.25)
        } else {
          lmcub <- lmquad
        }
        Fp.hats[i] <- lmquad[1]/h #1 if constrained, 2 if unconstrained
        Fpp.hats[i] <- 2*lmcub[2]/h^2 #2 if constrained, 3 if unconstrained
      } else { #use probit link
        #suppress <<glm.fit: fitted probabilities numerically 0 or 1 occurred>>
       suppressWarnings({
        if (length(Ybelowloc)<40) { #arbitrary threshold--but don't want cubic if small sample
#           tmp <- glm(Ybelowloc~Xdes2+Xdes3, family=binomial(link="probit"))$coefficients
          tmp <- glm(Ybelowloc~Xdes2+Xdes3+offset(qnorm(rep.int(p,length(Xdes2))))+0, family=binomial(link="probit"))$coefficients #constrained to fit (0,0.25)
        } else {
#           tmp <- glm(Ybelowloc~Xdes2+Xdes3+Xdes4, family=binomial(link="probit"))$coefficients
          tmp <- glm(Ybelowloc~Xdes2+Xdes3+Xdes4+offset(qnorm(rep.int(p,length(Xdes2))))+0, family=binomial(link="probit"))$coefficients #constrained
        }
       })
#         Fpp.hats[i] <- dnorm(tmp[1])*(2*tmp[3]-tmp[1]*tmp[2]^2)/h^2 #unconstrained
#         Fp.hats[i] <- dnorm(tmp[1])*tmp[2]/h #unconstrained
        Fpp.hats[i] <- dnorm(qnorm(p))*(2*tmp[2]-qnorm(p)*tmp[1]^2)/h^2 #constrained
        Fp.hats[i] <- dnorm(qnorm(p))*tmp[1]/h #constrained
      }
    }
  } else { #LOCAL==FALSE; use splines instead
    tmp.res <- qsreg(X,Y,alpha=p) #or try rqss
    #
    df.qsreg<- max(3,min(round(tmp.res$trace[1]),floor(length(Y)/3),length(unique(X))-6))
    tmp <- tryCatch(
      {desmat0 <- bs(X,df=df.qsreg,degree=3);tmp.fit.rq <- rq(Y~desmat0, tau=p);NULL}, 
  error=function(ME) {df.qsreg<-max(3,floor(df.qsreg/2)); desmat0 <- bs(X,df=df.qsreg,degree=3); tmp.fit.rq <- rq(Y~desmat0, tau=p); return(list(a=df.qsreg,b=tmp.fit.rq,c=desmat0))}
    )
    if (!is.null(tmp)) { df.qsreg<-tmp$a;tmp.fit.rq<-tmp$b;desmat0<-tmp$c } #tmp is not NULL if error code above ran
    desmat1 <- mybs(X,df=df.qsreg,degree=3,derivs=1)
    desmat2 <- mybs(X,df=df.qsreg,degree=3,derivs=2)
    tmp.app.fn0 <- approxfun(x=X,y=cbind(rep.int(1,length(X)),desmat0) %*% tmp.fit.rq$coef)
    Y.hats   <- tmp.app.fn0(x0s)
    #
    Fp.hats <- rep.int(0,length(Y.hats)) #estimated below
    Fpp.hats <- rep.int(0,length(Y.hats)) #estimated below
#     f.Y.hats <- fitted(npcdens(bws=npcdensbw(formula=Y~X,bwmethod="normal-reference",ckertype=KERNEL.TYPE,ckerorder=2),tau=p,cxkertype=KERNEL.TYPE,cxkorder=2,cykertype=KERNEL.TYPE,cykorder=2,exdat=x0s,eydat=Y.hats))
    for (i in 1:length(x0s)) {
      if (!is.na(Y.hats[i])) {
        Ybelow <- as.integer(Y<=Y.hats[i])
        tmp.fit.lm <- lm(Ybelow~desmat0) #unconstrained
        tmp.app.fn1 <- approxfun(x=X,y=cbind(rep.int(0,length(X)),desmat1) %*% tmp.fit.lm$coef)
        tmp.app.fn2 <- approxfun(x=X,y=cbind(rep.int(0,length(X)),desmat2) %*% tmp.fit.lm$coef)
#         tmp.fit.lm <- lm(I(Ybelow-p)~desmat0+0) #constrained
#         tmp.app.fn1 <- approxfun(x=X,y=desmat1 %*% tmp.fit.lm$coef)
#         tmp.app.fn2 <- approxfun(x=X,y=desmat2 %*% tmp.fit.lm$coef)
        Fp.hats[i] <- tmp.app.fn1(x0s[i])
        Fpp.hats[i] <- tmp.app.fn2(x0s[i])
      } else {
        Fp.hats[i] <- Fpp.hats[i] <- NA
      }
    }
  }
  #
  return(list(Y.hats=Y.hats,f.X.hats=f.X.hats,fp.X.hats=fp.X.hats,Fp.hats=Fp.hats,Fpp.hats=Fpp.hats))
}


########################
#                      #
# Bandwidths           #
#                      #
########################
quantile.inf.np.bandwidths <- function(p,n,X,x0s,f.X.hats,fp.X.hats,Fp.hats,Fpp.hats,ONESIDED,zz1,DENOM.MIN=0) {
    NUM.X0 <- length(x0s)
    hs <- matrix(0,NUM.X0,1);  h.mins <- matrix(0,NUM.X0,1)
    for (i in 1:NUM.X0) {
#       absdevs <- abs(X-x0s[i]) #transform X to abs distances from x0
#       absdevs.sort <- sort(absdevs) 
      h.max.bd <- min(max(X)-x0s[i],x0s[i]-min(X)) #if h were bigger, would hit boundary, inducing bias
      #Objects already estimated
        f.X.hat <- f.X.hats[i];  fp.X.hat <- fp.X.hats[i]
        Fp.hat  <-  Fp.hats[i];  Fpp.hat  <- Fpp.hats[i]
      if (ONESIDED==0) { #i.e. two-sided inference
        #CPE-optimal bandwidth
        hs[i] <- hhat.2s.fn(p=p,n=n,z=zz1,Fpp=Fpp.hat,Fp=Fp.hat,f.X=f.X.hat,fp.X=fp.X.hat,DENOM.MIN=DENOM.MIN)
        #As n grows, move toward MSE-optimal bandwidth, or rather n^(-1/20) undersmoothing thereof. Assumes d=1.
          # kap <- (2/((2+1)*(4+1)))*pmax(0,1-4*n^(-1/5))*(0.95)
          hs[i] <- hs[i] * pmax(1,n/1000)^(5/60) #(1/3)*(2+n^kap)
        #if running against boundary, use h.max.bd
        hs[i] <- min(hs[i],h.max.bd)
      } else { #ONESIDED!=0, i.e. one-sided CI
        stop('Not yet implemented: one-sided inference')
        hs[i] <- hhat.1s.fn(p,n,zz1,Fpp.hat,Fp.hat,f.X.hat,fp.X.hat,DENOM.MIN=DENOM.MIN)
        # COMPARED TO BOUNDARY........
      } #end of if/else for ONESIDED
    } #end of for-loop over x0s[i]
  return(hs)
}


########################
#                      #
# Bandwidth correction #
#                      #
########################
#Bandwidth correction: don't let x0s[2]'s window extend left of x0s[1]'s, nor right of x0s[3]'s
quantile.inf.np.bandwidth.correct <- function(hs,x0s,HCOR) {
  if (length(hs)==1) return(hs)
  noNAind <- !is.na(hs)
  hs.old <- rep.int(0,length(hs))
  while (any(hs[noNAind] != hs.old[noNAind])) {
    hs.old <- hs
    for (i in 2:length(hs)) {
      if (noNAind[i] && noNAind[i-1]) {
        hs[i] <- min(hs[i], hs[i-1] + HCOR*(x0s[i]-x0s[i-1]))
      }
      iend <- length(hs)-i+1
      if (noNAind[iend] && noNAind[iend+1]) {
        hs[iend] <- min(hs[iend], hs[iend+1] + HCOR*(x0s[iend+1]-x0s[iend]))
      }
    }
  }
  return(hs)
}


########################
#                      #
# Compute CIs          #
#                      #
########################
quantile.inf.np.compute.CIs <- function(Y,X,x0s,p,alpha,ONESIDED,METHOD.TYPE,PSI,hs,BETA.BLOCK.SIZE=10^4,BETA.BLOCKS=5, NORM.APPROX.CAL=FALSE, NORM.APPROX.Q=FALSE, NUM.INT=TRUE) {
  	NUM.X0 <- length(x0s)
    CIls <- CIus <- Nns <- rep(NA,NUM.X0)
    dim(CIls) <- dim(CIus) <- c(NUM.X0,1)
    if (METHOD.TYPE=='joint' && length(p)>1) {
      CIls <- CIus <- matrix(NA,nrow=NUM.X0,ncol=length(p))
    }
    if (METHOD.TYPE=='qte') {
      Nns <- matrix(NA,nrow=NUM.X0,ncol=2)
      # SHOULD DO ANYWAY IN QUANTILE.INF(): if (length(p)!=1) NORM.APPROX <- FALSE
    }
    ok <- rep.int(NA,NUM.X0)
    for (i in 1:NUM.X0) {
      x0 <- x0s[i]
      if (METHOD.TYPE=='qte') {
        Y.Ch <- list(t=Y$t[(x0-hs[i,1]<=X$t) & (X$t<=x0+hs[i,1])], c=Y$c[(x0-hs[i,2]<=X$c) & (X$c<=x0+hs[i,2])])
        Nns[i,] <- c(length(Y.Ch$t),length(Y.Ch$c))
      } else {
        Y.Ch <- Y[(x0-hs[i]<=X) & (X<=x0+hs[i])] #rem: d=1
        Nns[i] <- length(Y.Ch)
      }
      useNA <- tryCatch( {
        tmp <- quantile.inf(X=Y.Ch,p=p,ONESIDED=ONESIDED,ALPHA=alpha,PSI=PSI,METHOD.TYPE=METHOD.TYPE,NORM.APPROX.CAL=NORM.APPROX.CAL,NORM.APPROX.Q=NORM.APPROX.Q,SPACING.FLAG=TRUE, BETA.BLOCK.SIZE=BETA.BLOCK.SIZE, BETA.BLOCKS=BETA.BLOCKS, NUM.INT=NUM.INT)
        CIls[i,] <- tmp$CI$lo;  CIus[i,] <- tmp$CI$hi 
        ok[i] <- (tmp$methname %in% c('Hutson','GoldmanKaplan'))
        FALSE } ,
        error=function(ME) {
          if (grepl("have p as close to zero",ME$message,fixed=TRUE) || grepl("1/n<=p<1-1/n",ME$message,fixed=TRUE) || grepl("Order statistic does not exist",ME$message,fixed=TRUE)) {
            return(TRUE)
          } else {
            if (METHOD.TYPE=='qte') {
              stop(sprintf("Error running quantile.inf() with Nn.pt.t=%d,Nn.pt.c=%d,p=%5g,hs.pt.t(i)=%8g,hs.pt.c(i)=%8g,i=%d,alpha=%g,ONESIDED=%d,Y.Ch$t=%s,Y.Ch$c=%s: %s", Nns[i,1],Nns[i,2],p,hs[i,1],hs[i,2],i,alpha,ONESIDED,paste0(sprintf("%4.2f",Y.Ch$t),collapse="|"),paste0(sprintf("%4.2f",Y.Ch$c),collapse="|"),ME$message))
            } else {
              stop(sprintf("Error running quantile.inf() with Nn.pt=%d,p=%5g,hs.pt(i)=%8g,i=%d: %s", Nns[i],p,hs[i],i,ME$message))
            }
          }
        } )
        if (useNA) { CIls[i,] <- CIus[i,] <- ok[i] <- NA }
    }
    return(list(CIls=CIls,CIus=CIus,Nns=Nns,ok=ok))
}


###################
#                 #
# Main function   #
#                 #
###################
# Y is a vector of observations of a continuous outcome variable
# X is a vector of observations of a continuous conditioning variable, or a matrix where the first column is the continuous conditioning variable and the remaining columns are observations of discrete conditioning variables.  If METHOD.TYPE is 'qte' then second column of X should be the binary treatment indicator.
# x0s is a vector of points of interest for the continuous conditioning variable.  All possible combinations of the discrete conditioning variables in X (if present) will also be used; if only a subset are desired, then such filtering should be done before passing X to function quantile.inf.np
# p is the quantile(s), a scalar between 0 and 1 (e.g., 0.5 for the median), or a vector in the case of inference on multiple quantiles jointly or on linear combinations of quantiles
# ALPHA is between 0 and 1 such that desired coverage probability is (1-ALPHA); e.g., ALPHA=0.05 for 95% confidence intervals
# JOINT.FLAG=TRUE to calculate both pointwise and joint intervals (if length(x0s)>1); FALSE for only pointwise
# UNIFORM.FLAG=TRUE to compute uniform confidence band, if ONESIDED=0 and METHOD.TYPE='single'; requires "quantreg" package.
# ONESIDED=0 for two-sided confidence intervals, 1 for upper one-sided and -1 for lower one-sided
# KERNEL.TYPE is for nonparametric conditional density estimation (for plug-in bandwidth); can be gaussian (default), epanechnikov, or uniform
# METHOD.TYPE is one of four options: 'single', 'qte', 'joint', or 'lincom', which correspond respectively to conditional inference on a single quantile, conditional quantile treatment effect inference (including treatment effects on linear combinations of quantiles), conditional inference on multiple quantiles jointly, and conditional inference on a linear combination of quantiles (such as the interquartile range)
# PSI is a vector of scalar weights for METHOD.TYPE 'lincom' or 'qte', with length(PSI)=length(p) and normalized to PSI[1]=1
# BETA.BLOCK.SIZE and BETA.BLOCKS are passed directly to quantile.inf(), determining how many beta distribution draws are used in calibration
# NORM.APPROX.CAL=TRUE uses the normal approximation for calibrating alpha-tilde for METHOD.TYPE 'qte' and length(p)=1, and has no effect otherwise
# NORM.APPROX.Q=TRUE uses the normal approximation to find the fractional order statistic index
# NUM.INT=TRUE uses numerical integration rather than simulation (FALSE) for QTE with length(p)==1
# hs.pt, hs.joint: to supply your own bandwidths (not generally recommended).  Note for QTE, hs.pt (and/or hs.joint) should be a list with length(hs.pt$t)==length(hs.pt$c)==length(x0s), otherwise just vectors with length(hs.pt)==length(x0s).  Note: currently restricted to be the same regardless of the values of the discrete elements of X.
# LOCAL: for plug-in bandwidth estimations, use local polynomial approach if LOCAL==TRUE, otherwise B-spline approach if FALSE; if NA, then use local for larger n and spline for smaller n.
quantile.inf.np <- 
  function(Y,X,x0s,p=NULL,ALPHA=0.05,JOINT.FLAG=TRUE,UNIFORM.FLAG=FALSE,ONESIDED=0,KERNEL.TYPE='gaussian',METHOD.TYPE='single',PSI=NULL,BETA.BLOCK.SIZE=10^4,BETA.BLOCKS=5, hs.pt=NULL, hs.joint=NULL, NORM.APPROX.CAL=NA, NORM.APPROX.Q=FALSE, NUM.INT=TRUE, LOCAL=TRUE) {
  
  #Input validity checking
  if (is.na(NORM.APPROX.CAL)) NORM.APPROX.CAL <- (0.2<min(p)) && (max(p)<0.8)
  if (!is.logical(LOCAL)) { stop("LOCAL must be TRUE, FALSE, or NA.") }
  validated <- quantile.inf.np.validate.inputs(Y,X,x0s,p,ALPHA,JOINT.FLAG,ONESIDED,KERNEL.TYPE,METHOD.TYPE,PSI,hs.pt,hs.joint)
  if (!is.logical(UNIFORM.FLAG) || is.na(UNIFORM.FLAG)) stop("UNIFORM.FLAG must be either TRUE or FALSE.")
  if (UNIFORM.FLAG) {
    if (ONESIDED==0 && METHOD.TYPE=='single') {
      tryCatch(require("quantreg"),
               warning=function(WWW) {warning("Cannot load package quantreg; will not compute uniform band."); UNIFORM.FLAG <- FALSE})
    } else {
      warning("Ignoring UNIFORM.FLAG=TRUE: only implemented for ONESIDED=0 and METHOD.TYPE='single'")
      UNIFORM.FLAG <- FALSE
    }
  }
  
  #Initialize variables
  hs.pt.raw <- hs.pt;  hs.joint.raw <- hs.joint
  Y.raw <- Y
  n <- length(Y)
  X.raw.raw <- X.raw <- X
  if (METHOD.TYPE=='qte') {
    tmp <- dim(X.raw)
    if (tmp[1]!=n) { X.raw<-t(X.raw); X.raw.raw<-t(X.raw.raw) }
    X.treatment.raw <- X[,2] #treatment indicator
    X.raw <- X[,-2] #w/o treatment
    Y <- list(t=Y.raw[!!X.treatment.raw],c=Y.raw[!X.treatment.raw])
    if (dim(X.raw.raw)[2]==2) { X <- list(t=X.raw[!!X.treatment.raw], c=X.raw[!X.treatment.raw]) }
  }
  X.discrete <- NULL;  X.discrete.n <- 1
  if (is.matrix(X.raw) && (dim(X.raw)[1]>1 || dim(X.raw)[2]>1)) {
    tmp <- dim(X.raw)
    if (tmp[1]!=n) { X.raw<-t(X.raw); X.raw.raw<-t(X.raw.raw) }
    X.cts.raw <- X.raw[,1]
    X <- X.cts.raw
    X.discrete.raw <- X.raw[,2:dim(X.raw)[2]]
    tmp.ndreg <- dim(X.raw)[2] - 1 #number of discrete regressors
    if (tmp.ndreg>0) {
      if (tmp.ndreg==1) {dim(X.discrete.raw)<-c(n,1)}
      X.discrete.unique <- unique(X.discrete.raw)
      if (tmp.ndreg==1) {dim(X.discrete.unique) <- c(length(X.discrete.unique),1)} #account for single discrete X
      if (length(X.discrete.unique)==tmp.ndreg) {dim(X.discrete.unique)<-c(1,length(X.discrete.unique))} #multiple discrete vars, but same value for all obs
      X.discrete.n <- dim(X.discrete.unique)[1]
    }
  }
  NUM.X0 <- length(x0s)

  zz1 <- qnorm(1-ALPHA,0,1)
  zz2 <- qnorm(1-ALPHA/2,0,1)
  p.hat <- p[which.min(p*(1-p))]
  if (METHOD.TYPE=='qte') { tmp <- list(t=sort(Y$t),c=sort(Y$c)) } else { tmp <- sort(Y) }
  alpha.hat <- switch(METHOD.TYPE,'single'=ALPHA,'joint'=ALPHA,'lincom'=,'qte'=suppressWarnings(quantile.inf.calibrate(X=tmp, p=p, ALPHA=ALPHA, ONESIDED=ONESIDED, METHOD.TYPE=METHOD.TYPE, PSI=PSI, NORM.APPROX=NORM.APPROX.CAL)))
  zz1.hat <- qnorm(1-alpha.hat,0,1)
  zz2.hat <- qnorm(1-alpha.hat/2,0,1)
  
#   pdf(NULL) #don't draw plot() until after calculations

  # Loop over different cells of discrete (non-treatment) conditioning variables
  x0s.raw <- x0s;  NUM.X0.RAW <- NUM.X0
  rets <- vector("list",X.discrete.n)
  for (i.discX in 1:X.discrete.n)  {
    # Get subset of Y and X if discrete conditioners
    if (X.discrete.n>1) {
      X.disc.i <- X.discrete.unique[i.discX,]
      tmp.ind <- apply(X.discrete.raw, 1, FUN=function(x) identical(x,X.disc.i))
      if (METHOD.TYPE=='qte') {
        Treatment <- X.treatment.raw[tmp.ind]
        Y <- list(t=Y.raw[tmp.ind & Treatment], c=Y.raw[tmp.ind & !Treatment])
        X <- list(t=X.cts.raw[tmp.ind & Treatment], c=X.cts.raw[tmp.ind & !Treatment])
      } else {
        Y <- Y.raw[tmp.ind]
        X <- X.cts.raw[tmp.ind]
      }
    }

    # Get sample size, check if large enough
    if (METHOD.TYPE=='qte') { 
      n <- list(t=length(Y$t),c=length(Y$c)) 
      if (min(p)<1/min(n$t,n$c) || max(p)>=1-1/min(n$t,n$c)) { #note: 1/0=Inf
        warning("Not enough observations for inference when discrete X vector has value = ",sprintf("%g,",X.disc.i)," so inference based on extreme value theory is recommended instead.")
        next
      }
    } else {
      n <- length(Y) 
      if (min(p)<1/n || max(p)>=1-1/n) {
        warning("Not enough observations for inference when discrete X vector has value = ",sprintf("%g,",X.disc.i)," so inference based on extreme value theory is recommended instead.")
        next
      }
    }

    #remove some x0s if necessary
    if (METHOD.TYPE=='qte') {
      x0s.raw.inds <- (x0s.raw>=max(min(X$t),min(X$c))) & (x0s.raw<=min(max(X$t),max(X$c)))
      x0s <- x0s.raw[x0s.raw.inds]
    } else {
      x0s.raw.inds <- (x0s.raw>=min(X)) & (x0s.raw<=max(X))
      x0s <- x0s.raw[x0s.raw.inds]
    }
    NUM.X0 <- length(x0s);  NUM.X0.raw <- length(x0s.raw)
    
    # Set up Bonferroni
    if (JOINT.FLAG) {
      ALPHA.BONF <- ALPHA/NUM.X0
      zz1.bonf <- qnorm(1-ALPHA.BONF,0,1)
      zz2.bonf <- qnorm(1-ALPHA.BONF/2,0,1)
      if (METHOD.TYPE=='qte') { tmp <- list(t=sort(Y$t),c=sort(Y$c)) } else { tmp <- sort(Y) }
      alpha.bonf.hat <- switch(METHOD.TYPE,'single'=ALPHA.BONF,'joint'=ALPHA.BONF,'lincom'=,'qte'=suppressWarnings(quantile.inf.calibrate(X=tmp, p=p, ALPHA=ALPHA.BONF, ONESIDED=ONESIDED, METHOD.TYPE=METHOD.TYPE, PSI=PSI, NORM.APPROX=NORM.APPROX.CAL)))
      zz1.bonf.hat <- qnorm(1-alpha.bonf.hat,0,1)
      zz2.bonf.hat <- qnorm(1-alpha.bonf.hat/2,0,1)
    }
    
    # Set up uniform band
    if (UNIFORM.FLAG) {
      ALPHA.UNIF <- NA
      if (n>4e3) {
        warning("Skipping uniform bands when subsample size n>4000")
      } else {
        if (n>1e3) {
          warning("Uniform band computation may take up to a few minutes for subsamples of size 1000<n<=4000")
        }
        g <- function(lam,y,x) AIC(rqss(y ~ qss(x, lambda = lam)),k = -1) #k=-1 does BIC (look @ code for AIC.rqss)
        lamstar <- optimize(g, interval = c(0.001, 0.5), x = X, y = Y) 
        rqss.fit <- rqss(Y~qss(X, lambda=lamstar$min), tau=p)
        cvu <- rqss.unif.cv(rqss.obj=rqss.fit, coverage=1-ALPHA, newdata=x0s)
        ALPHA.UNIF <- 2*pnorm(-cvu)
        zz2.unif <- qnorm(1-ALPHA.UNIF/2,0,1)
      }
    } 
    
    # Initialize res.gk to hold results
    NAvec <- rep.int(NA,NUM.X0.raw)
    NAmat <- rep.int(NA,length(p)*NUM.X0.raw)
    dim(NAmat) <- c(NUM.X0.raw,length(p))
    if (METHOD.TYPE=='qte') {
      res.gk <- list(CI.pt.lo=NAvec,CI.pt.hi=NAvec,h.pt.t=NAvec,h.pt.c=NAvec,N.pt.t=NAvec,N.pt.c=NAvec,ok.pt=NAvec)
    } else if (METHOD.TYPE=='joint') {
      res.gk <- list(CI.pt.lo=NAmat,CI.pt.hi=NAmat,h.pt=NAvec,N.pt=NAvec,ok.pt=NAvec)
    } else {
      res.gk <- list(CI.pt.lo=NAvec,CI.pt.hi=NAvec,h.pt=NAvec,N.pt=NAvec,ok.pt=NAvec)
    }
    if (JOINT.FLAG) {
      if (METHOD.TYPE=='qte') {
        res.gk <- c(res.gk,list(CI.bonf.lo=NAvec,CI.bonf.hi=NAvec,h.bonf.t=NAvec,h.bonf.c=NAvec,N.bonf.t=NAvec,N.bonf.c=NAvec,ok.bonf=NAvec))
      } else if (METHOD.TYPE=='joint') {
        res.gk <- c(res.gk,list(CI.bonf.lo=NAmat,CI.bonf.hi=NAmat,h.bonf=NAvec,N.bonf=NAvec,ok.bonf=NAvec))
      } else {
        res.gk <- c(res.gk,list(CI.bonf.lo=NAvec,CI.bonf.hi=NAvec,h.bonf=NAvec,N.bonf=NAvec,ok.bonf=NAvec))
      }
    }
    if (UNIFORM.FLAG) {
      res.gk <- c(res.gk,list(CI.unif.lo=NAvec,CI.unif.hi=NAvec,h.unif=NAvec,N.unif=NAvec,ok.unif=NAvec))
    }

    # Estimate objects for plug-in bandwidth
    if (UNIFORM.FLAG || is.null(hs.pt.raw) || (is.null(hs.joint.raw) && JOINT.FLAG)) {
      if (METHOD.TYPE=='qte') {
        list[Y.hats.t,f.X.hats.t,fp.X.hats.t,Fp.hats.t,Fpp.hats.t] <- quantile.inf.np.est.objs(Y=Y$t,X=X$t,x0s=x0s,p=p.hat,KERNEL.TYPE=KERNEL.TYPE,LOCAL=LOCAL)
        list[Y.hats.c,f.X.hats.c,fp.X.hats.c,Fp.hats.c,Fpp.hats.c] <- quantile.inf.np.est.objs(Y=Y$c,X=X$c,x0s=x0s,p=p.hat,KERNEL.TYPE=KERNEL.TYPE,LOCAL=LOCAL)
      } else {
        list[Y.hats,f.X.hats,fp.X.hats,Fp.hats,Fpp.hats] <- quantile.inf.np.est.objs(Y=Y,X=X,x0s=x0s,p=p.hat,KERNEL.TYPE=KERNEL.TYPE,LOCAL=LOCAL)
      }
    }

    # Calculate, correct, and store bandwidths
    if (METHOD.TYPE=='qte') {
      hs.pt.t <- hs.pt.c <- rep.int(NA,NUM.X0) #*not* .raw
      if (is.null(hs.pt.raw)) {
        if (ONESIDED==0) { zz <- zz2.hat } else { zz <- zz1.hat }
        hs.pt.t <- quantile.inf.np.bandwidths(p=p.hat,n=n$t,X=X$t,x0s=x0s,f.X.hats=f.X.hats.t,fp.X.hats=fp.X.hats.t,Fp.hats=Fp.hats.t,Fpp.hats=Fpp.hats.t,ONESIDED=ONESIDED,zz1=zz)
        hs.pt.c <- quantile.inf.np.bandwidths(p=p.hat,n=n$c,X=X$c,x0s=x0s,f.X.hats=f.X.hats.c,fp.X.hats=fp.X.hats.c,Fp.hats=Fp.hats.c,Fpp.hats=Fpp.hats.c,ONESIDED=ONESIDED,zz1=zz)
        # Adjustment to CPE-optimal rate: see Cor. 4 in Kaplan (2012)
        if (ONESIDED==0) {
          hs.pt.t <- n$t^(-2/((6*2+7*1)*(2+1)))*hs.pt.t
        } else {
          hs.pt.c <- n$c^(4*2/((6*2+7*1)*(2*2+3*1)))*hs.pt.c
        }
      } else {
        hs.pt.t <- hs.pt.raw$t[x0s.raw.inds];  hs.pt.c <- hs.pt.raw$c[x0s.raw.inds]
      }
      res.gk$h.pt.t[x0s.raw.inds] <- quantile.inf.np.bandwidth.correct(hs=hs.pt.t,x0s=x0s,HCOR=HCOR)
      res.gk$h.pt.c[x0s.raw.inds] <- quantile.inf.np.bandwidth.correct(hs=hs.pt.c,x0s=x0s,HCOR=HCOR)
      if (JOINT.FLAG) {
        if (is.null(hs.joint.raw)) {
          if (ONESIDED==0) { zz <- zz2.bonf.hat } else { zz <- zz1.bonf.hat }
          hs.bonf.t <- quantile.inf.np.bandwidths(p=p.hat,n=n$t,X=X$t,x0s=x0s,f.X.hats=f.X.hats.t,fp.X.hats=fp.X.hats.t,Fp.hats=Fp.hats.t,Fpp.hats=Fpp.hats.t,ONESIDED=ONESIDED,zz1=zz)
          hs.bonf.c <- quantile.inf.np.bandwidths(p=p.hat,n=n$c,X=X$c,x0s=x0s,f.X.hats=f.X.hats.c,fp.X.hats=fp.X.hats.c,Fp.hats=Fp.hats.c,Fpp.hats=Fpp.hats.c,ONESIDED=ONESIDED,zz1=zz)
          # Adjustment to CPE-optimal rate: see Cor. 4 in Kaplan (2012)
          if (ONESIDED==0) {
            hs.bonf.t <- n$t^(-2/((6*2+7*1)*(2+1)))*hs.bonf.t
          } else {
            hs.bonf.c <- n$c^(4*2/((6*2+7*1)*(2*2+3*1)))*hs.bonf.c
          }
        } else {
          hs.bonf.t <- hs.joint.raw$t[x0s.raw.inds];  hs.bonf.c <- hs.joint.raw$c[x0s.raw.inds]
        }
        res.gk$h.bonf.t[x0s.raw.inds] <- quantile.inf.np.bandwidth.correct(hs=hs.bonf.t,x0s=x0s,HCOR=HCOR)
        res.gk$h.bonf.c[x0s.raw.inds] <- quantile.inf.np.bandwidth.correct(hs=hs.bonf.c,x0s=x0s,HCOR=HCOR)
      }
    } else {
      if (is.null(hs.pt.raw)) {
        if (ONESIDED==0) { zz <- zz2.hat } else { zz <- zz1.hat }
        hs.pt <- quantile.inf.np.bandwidths(p=p.hat,n=n,X=X,x0s=x0s,f.X.hats=f.X.hats,fp.X.hats=fp.X.hats,Fp.hats=Fp.hats,Fpp.hats=Fpp.hats,ONESIDED=ONESIDED,zz1=zz)
        # Adjustment to CPE-optimal rate: see Cor. 4 in Kaplan (2012).  No adj needed for 'joint' since same optimal rate as for 'single'
        if (METHOD.TYPE=='lincom') {
          if (ONESIDED==0) {
            hs.pt <- n^(-2/((6*2+7*1)*(2+1)))*hs.pt
          } else {
            hs.pt <- n^(4*2/((6*2+7*1)*(2*2+3*1)))*hs.pt
          }
        }
      } else {
        hs.pt <- hs.pt.raw[x0s.raw.inds]
      }
      res.gk$h.pt[x0s.raw.inds] <- quantile.inf.np.bandwidth.correct(hs=hs.pt,x0s=x0s,HCOR=HCOR)
      if (JOINT.FLAG) {
        if (is.null(hs.joint.raw)) {
          if (ONESIDED==0) { zz <- zz2.bonf.hat } else { zz <- zz1.bonf.hat }
          hs.bonf <- quantile.inf.np.bandwidths(p=p.hat,n=n,X=X,x0s=x0s,f.X.hats=f.X.hats,fp.X.hats=fp.X.hats,Fp.hats=Fp.hats,Fpp.hats=Fpp.hats,ONESIDED=ONESIDED,zz1=zz)
          if (METHOD.TYPE=='lincom') {
            if (ONESIDED==0) {
              hs.bonf <- n^(-2/((6*2+7*1)*(2+1)))*hs.bonf
            } else {
              hs.bonf <- n^(4*2/((6*2+7*1)*(2*2+3*1)))*hs.bonf
            }
          }
        } else {
          hs.bonf <- hs.joint.raw[x0s.raw.inds]
        }
        res.gk$h.bonf[x0s.raw.inds] <- quantile.inf.np.bandwidth.correct(hs=hs.bonf,x0s=x0s,HCOR=HCOR)
      }
      if (UNIFORM.FLAG) {
        hs.unif <- quantile.inf.np.bandwidths(p=p.hat,n=n,X=X,x0s=x0s,f.X.hats=f.X.hats,fp.X.hats=fp.X.hats,Fp.hats=Fp.hats,Fpp.hats=Fpp.hats,ONESIDED=ONESIDED,zz1=zz2.unif)
        res.gk$h.unif[x0s.raw.inds] <- quantile.inf.np.bandwidth.correct(hs=hs.unif,x0s=x0s,HCOR=HCOR)
      }
    }

    # Compute GoldmanKaplan CIs
    tmph <- res.gk$h.pt[x0s.raw.inds]
    if (METHOD.TYPE=='qte') { tmph <- data.frame(t=res.gk$h.pt.t[x0s.raw.inds], c=res.gk$h.pt.c[x0s.raw.inds]) }
    list[res.gk$CI.pt.lo[x0s.raw.inds],res.gk$CI.pt.hi[x0s.raw.inds],tmpN,res.gk$ok.pt[x0s.raw.inds]] <- tryCatch(quantile.inf.np.compute.CIs(Y=Y,X=X,x0s=x0s,p=p,alpha=ALPHA,ONESIDED=ONESIDED,METHOD.TYPE=METHOD.TYPE,PSI=PSI,hs=tmph,NORM.APPROX.CAL=NORM.APPROX.CAL,NORM.APPROX.Q=NORM.APPROX.Q,NUM.INT=NUM.INT),error=function(QE) {warning("Uncaught error while computing CIs: ", QE$message, call. = FALSE); NA } ) #set all to NULL if error, output message as warning & continue on
    if (any(is.na(tmpN))) next
    if (METHOD.TYPE=='qte') { res.gk$N.pt.t[x0s.raw.inds] <- tmpN[,1]; res.gk$N.pt.c[x0s.raw.inds] <- tmpN[,2] } else { res.gk$N.pt[x0s.raw.inds] <- tmpN }
    if (JOINT.FLAG) {
      tmph <- res.gk$h.bonf[x0s.raw.inds]
      if (METHOD.TYPE=='qte') { tmph <- data.frame(t=res.gk$h.bonf.t[x0s.raw.inds], c=res.gk$h.bonf.c[x0s.raw.inds]) }
      list[res.gk$CI.bonf.lo[x0s.raw.inds],res.gk$CI.bonf.hi[x0s.raw.inds],tmpN,res.gk$ok.bonf[x0s.raw.inds]] <- tryCatch(quantile.inf.np.compute.CIs(Y=Y,X=X,x0s=x0s,p=p,alpha=ALPHA.BONF,ONESIDED=ONESIDED,METHOD.TYPE=METHOD.TYPE,PSI=PSI,hs=tmph,NORM.APPROX.CAL=NORM.APPROX.CAL,NORM.APPROX.Q=NORM.APPROX.Q,NUM.INT=NUM.INT),error=function(QE) { cat(sprintf("min(p)=%g,max(p)=%g,ALPHA.BONF=%g,ONESIDED=%g,METHOD.TYPE=%s,NORM.APPROX.CAL=%s,NORM.APPROX.Q=%s,NUM.INT=%s",min(p),max(p),ALPHA.BONF,ONESIDED,METHOD.TYPE,NORM.APPROX.CAL,NORM.APPROX.Q,NUM.INT),file="",append=T,sep='\n');flush.console(); warning("Uncaught error while computing joint CIs: ",QE$message,call.=F); NA } )
      if (any(is.na(tmpN))) next
      if (METHOD.TYPE=='qte') { res.gk$N.bonf.t[x0s.raw.inds] <- tmpN[,1]; res.gk$N.bonf.c[x0s.raw.inds] <- tmpN[,2] } else { res.gk$N.bonf[x0s.raw.inds] <- tmpN }
    }
    if (UNIFORM.FLAG) {
      tmph <- res.gk$h.unif[x0s.raw.inds]
      list[res.gk$CI.unif.lo[x0s.raw.inds],res.gk$CI.unif.hi[x0s.raw.inds],tmpN,res.gk$ok.unif[x0s.raw.inds]] <- tryCatch(quantile.inf.np.compute.CIs(Y=Y,X=X,x0s=x0s,p=p,alpha=ALPHA.UNIF,ONESIDED=ONESIDED,METHOD.TYPE=METHOD.TYPE,PSI=PSI,hs=tmph,NORM.APPROX.CAL=NORM.APPROX.CAL,NORM.APPROX.Q=NORM.APPROX.Q,NUM.INT=NUM.INT),error=function(QE) { warning("Uncaught error while computing uniform band: ",QE$message,call.=F); NA } )
      if (any(is.na(tmpN))) next
      res.gk$N.unif[x0s.raw.inds] <- tmpN
    }
    
    # Store results for return
    tmp <- list(Y.cts=Y,X.cts=X,x0s=x0s.raw,p=p,ALPHA=ALPHA, CI.pointwise.lows=res.gk$CI.pt.lo, CI.pointwise.highs=res.gk$CI.pt.hi, pointwise.method.Lstat=res.gk$ok.pt)
    if (METHOD.TYPE=='qte') {
      tmp <- c(tmp,list(bandwidth.pointwise.treat=res.gk$h.pt.t, bandwidth.pointwise.control=res.gk$h.pt.c, effective.N.pointwise.treat=res.gk$N.pt.t, effective.N.pointwise.control=res.gk$N.pt.c))
    } else {
      tmp <- c(tmp,list(bandwidth.pointwise=res.gk$h.pt, effective.N.pointwise=res.gk$N.pt))
    }
    if (JOINT.FLAG) {
      tmp.bonf <- list(CI.joint.lows=res.gk$CI.bonf.lo, CI.joint.highs=res.gk$CI.bonf.hi, joint.method.Lstat=res.gk$ok.bonf)
      if (METHOD.TYPE=='qte') {
        tmp.bonf <- c(tmp.bonf, list(bandwidth.joint.treat=res.gk$h.bonf.t, bandwidth.joint.control=res.gk$h.bonf.c, effective.N.joint.treat=res.gk$N.bonf.t, effective.N.joint.control=res.gk$N.bonf.c))
      } else {
        tmp.bonf <- c(tmp.bonf, list(bandwidth.joint=res.gk$h.bonf, effective.N.joint=res.gk$N.bonf))
      }
      tmp <- c(tmp,tmp.bonf)
    }
    if (UNIFORM.FLAG) {
      tmp.unif <- list(CI.uniform.lows=res.gk$CI.unif.lo, CI.uniform.highs=res.gk$CI.unif.hi, uniform.method.Lstat=res.gk$ok.unif)
      tmp.unif <- c(tmp.unif, list(bandwidth.uniform=res.gk$h.unif, effective.N.uniform=res.gk$N.unif))
      tmp <- c(tmp,tmp.unif)
    }
    if (X.discrete.n>1) { tmp$X.discrete <- X.disc.i }
    rets[[i.discX]] <- tmp
  } #end of loop over X.discrete.n

  for (i in 1:X.discrete.n) {
  	if (is.null(rets[[i]])) { next }
    if (any(is.na(c(rets[[i]]$CI.pointwise.lows,rets[[i]]$CI.pointwise.highs))) || (JOINT.FLAG && any(is.na(c(rets[[i]]$CI.pointwise.lows,rets[[i]]$CI.pointwise.highs))))) { warning("Effective sample size too small at some x0: NA values stored.  For those x0, inference based on extreme value theory may be necessary to avoid large bias from larger bandwidth (effective sample size).",call.=F); break }
  }

#   dev.off() #ok to plot normally again--turn off pdf(file=NULL)
  return(rets)
} #end of main function quantile.inf.np()


###################
#                 #
# Plot function 1 #
#                 #
###################
# NOTE: you may need to open device (X11(), quartz(), etc.) *before* calling this function.
plot.quantile.inf.np <- function(rets,plot.data=TRUE,CI.interpolate=FALSE,CI.type='pointwise',x.axis.title="X",y.axis.title="Y",title=NULL, ylim=NULL, xlim=NULL, legpos='bottomright', datacex=0.4, CIlwd=2) {
  # Validate argument values
  if (is.null(CI.type) || !(CI.type %in% c('pointwise','joint','both'))) stop(sprintf("Argument CI.type must be 'pointwise', 'joint', or 'both'; current value is %s",CI.type))
  if (!is.logical(plot.data) || is.na(plot.data)) stop("Argument plot.data must be either TRUE or FALSE")
  if (!is.logical(CI.interpolate) || is.na(CI.interpolate)) stop("Argument CI.interpolate must be either TRUE or FALSE")

  Y.cts.all <- unlist(lapply(rets,function(l)l$Y.cts))
  X.cts.all <- unlist(lapply(rets,function(l)l$X.cts))
  for (i.ret in 1:length(rets)) {
    #Load values from rets
    reti <- rets[[i.ret]]
    pt.lo <- reti$CI.pointwise.lows
    pt.hi <- reti$CI.pointwise.highs
    jt.lo <- reti$CI.joint.lows
    jt.hi <- reti$CI.joint.highs
    if ((is.null(jt.lo) || is.null(jt.hi)) && CI.type %in% c('joint','both')) stop(sprintf("CI.type is either 'joint' or 'both', but rets[[%d]]$CI.joint.lows (or highs) is NULL",i.ret))
    X.cts <- reti$X.cts;  Y.cts <- reti$Y.cts
    x0s <- reti$x0s;  ALPHA <- reti$ALPHA;  p <- reti$p

    # Draw graph: data+Lstat
    subtitle=NULL
    if (length(rets)>1) {
      subtitle <- paste0("Discrete X = ",paste0(sprintf("%g",reti$X.discrete),collapse=","))
    }
#     quartz()
    par(family="serif")
    #plot (not drawn) first w/ Y.cts.all and X.cts.all to keep same scaling for all graphs over different discrete combos
    if (is.null(title)) title <- sprintf("Lstat %g%% CI for %g-quantile",(1-ALPHA)*100,p)
    tmp <- rep.int(NA,4*length(x0s))
    if (CI.type=='pointwise' || CI.type=='both') {
      tmp[1:length(x0s)] <- pt.lo;  tmp[(length(x0s)+1):(2*length(x0s))] <- pt.hi
    }
    if (CI.type=='joint' || CI.type=='both') {
      tmp[(2*length(x0s)+1):(3*length(x0s))] <- jt.lo
      tmp[(3*length(x0s)+1):(4*length(x0s))] <- jt.hi
    }
    plot(x=c(X.cts.all,x0s,x0s,x0s,x0s),y=c(Y.cts.all,tmp),
         type="n",pch=1,col=1, 
         main=title, sub=subtitle, xlab=x.axis.title, ylab=y.axis.title, 
         cex.main=2.0, mgp=c(2.1,0.5,0), cex.axis=1.5, cex.lab=2, 
         xlim=xlim, ylim=ylim) 
    #
    if (plot.data) {
      if (is.list(X.cts)) {
        points(X.cts$t,Y.cts$t,type='p',pch=84,col=1,cex=datacex)
        points(X.cts$c,Y.cts$c,type='p',pch=67,col=1,cex=datacex)
      } else {
        points(X.cts,Y.cts,type='p',pch=16,col=1,cex=datacex)
      }
    }
    tmp <- ifelse(CI.interpolate,'o','p')
    if (CI.type=='pointwise' || CI.type=='both') {
      points(x0s,pt.lo,type=tmp,col=3,pch=3,lwd=CIlwd,lty=3)
      points(x0s,pt.hi,type=tmp,col=3,pch=3,lwd=CIlwd,lty=3)
    }
    if (CI.type=='joint' || CI.type=='both') {
      points(x0s,jt.lo,type=tmp,col=4,pch=4,lwd=CIlwd,lty=3)
      points(x0s,jt.hi,type=tmp,col=4,pch=4,lwd=CIlwd,lty=3)
    }
    #
    CIlty <- 0;  if (CI.interpolate) {CIlty <- 3} 
    legtxt <- NULL;  legcol <- NULL;  legpch <- NULL;  leglty <- NULL
    if (plot.data) {
      if (is.list(X.cts)) {
        legtxt<-c(legtxt,'Data (treated)','Data (control)');legcol<-c(legcol,1,1);legpch<-c(legpch,84,67);leglty<-c(leglty,NA,NA)
      } else {
        legtxt<-c(legtxt,'Data');legcol<-c(legcol,1);legpch<-c(legpch,16);leglty<-c(leglty,0)
      }
    }
    if (CI.type=='pointwise' || CI.type=='both') {legtxt<-c(legtxt,'Pointwise');legcol<-c(legcol,3);legpch<-c(legpch,3);leglty<-c(leglty,CIlty)}
    if (CI.type=='joint' || CI.type=='both') {legtxt<-c(legtxt,'Joint');legcol<-c(legcol,4);legpch<-c(legpch,4);leglty<-c(leglty,CIlty)}
    legend(legpos, legend=legtxt, inset=0.01, col=legcol, pch=legpch, lty=leglty, lwd=CIlwd, cex=1.5, y.intersp=0.9, bg='white')
  } #end for loop over length(rets), i.e. number of discrete combos
}



###################
#                 #
# Plot function 2 #
#                 #
###################
#Single graph with multiple CIs (per discrete X value)
# NOTE: you may need to open device (X11, quartz, etc.) *before* calling this function.
plot2.quantile.inf.np <- function(rets,plot.data=FALSE,CI.interpolate=TRUE,CI.type='pointwise',x.axis.title="X",y.axis.title="Y",title=NULL, ylim=NULL, xlim=NULL, legpos='bottomright', datacex=0.4, CIlwd=2) {
  # Validate arguments
  if (is.null(CI.type) || !(CI.type %in% c('pointwise','joint'))) stop("CI.type may only be 'pointwise' or 'joint' for plot2.quantile.inf.np.  (For plot.quantile.inf.np, 'both' is also accepted.)")
  if (!is.logical(plot.data) || is.na(plot.data)) stop("Argument plot.data must be either TRUE or FALSE")
  if (!is.logical(CI.interpolate) || is.na(CI.interpolate)) stop("Argument CI.interpolate must be either TRUE or FALSE")

  Y.cts.all <- unlist(lapply(rets,function(l)l$Y.cts)) #ok w/ NULL values, and includes $t and $c
  X.cts.all <- unlist(lapply(rets,function(l)l$X.cts))
  tryCatch({
    Y.cts.all.t <- unlist(lapply(rets,function(l)l$Y.cts$t))
    X.cts.all.t <- unlist(lapply(rets,function(l)l$X.cts$t))
    Y.cts.all.c <- unlist(lapply(rets,function(l)l$Y.cts$c))
    X.cts.all.c <- unlist(lapply(rets,function(l)l$X.cts$c))
  }, error=function(UE)NULL)
  #
  pt.lo <- vector("list",length(rets))
  pt.hi <- vector("list",length(rets))
  jt.lo <- vector("list",length(rets))
  jt.hi <- vector("list",length(rets))
  #
  p <- NULL;  ALPHA <- NULL
  x0s <- vector("list",length(rets))
  for (tmpi in 1:length(rets)) {
    pt.lo[[tmpi]] <- rets[[tmpi]]$CI.pointwise.lows
    pt.hi[[tmpi]] <- rets[[tmpi]]$CI.pointwise.highs
    jt.lo[[tmpi]] <- rets[[tmpi]]$CI.joint.lows
    jt.hi[[tmpi]] <- rets[[tmpi]]$CI.joint.highs
    x0s[[tmpi]]<-rets[[tmpi]]$x0s
    if (!is.null(rets[[tmpi]]$p) && !is.na(rets[[tmpi]]$p)) { p <- rets[[tmpi]]$p }
    if (!is.null(rets[[tmpi]]$ALPHA) && !is.na(rets[[tmpi]]$ALPHA)) { ALPHA <- rets[[tmpi]]$ALPHA }
  }
  if (is.null(p)) stop("No value of p found in rets other than NULL or NA")
  if (is.null(ALPHA)) stop("No value of p found in rets other than NULL or NA")
  if (length(p)>1) stop("Not yet implemented: plotting when length(p)>1")

  #Draw graph
  noninf <- (abs(Y.cts.all)!=Inf & !is.na(Y.cts.all))
  X.cts.all <- X.cts.all[noninf]
  Y.cts.all <- Y.cts.all[noninf]
  if (CI.type=='pointwise') { plo<-pt.lo; phi<-pt.hi } else { plo<-jt.lo; phi<-jt.hi }
  if (!plot.data) { #used only for scaling plot in initial plot() call
    tmp <- c(unlist(plo),unlist(phi))
    Y.cts.all <- range(tmp[abs(tmp)!=Inf & !is.na(tmp)])
    X.cts.all <- range(x0s)
  }
#   quartz()
  par(family="serif")
  #plot (not drawn) first w/ Y.cts.all and X.cts.all to keep same scaling for all graphs over different discrete combos
  if (is.null(title)) title<-sprintf("Lstat %g%% CI for %g-quantile",(1-ALPHA)*100,p)
  plot(x=X.cts.all,y=Y.cts.all,type="n",pch=1,col=1, main=title, sub=NULL, xlab=x.axis.title, ylab=y.axis.title, cex.main=2.0, mgp=c(2.1,0.5,0), cex.axis=1.5, cex.lab=2, xlim=xlim,ylim=ylim) 
  #
  legtxt<-NULL; legcol<-NULL; legpch<-NULL; leglty<-NULL
  if (plot.data) {
    if (is.list(X.cts)) {
      points(X.cts$t,Y.cts$t,type='p',pch=84,col=1,cex=datacex)
      points(X.cts$c,Y.cts$c,type='p',pch=67,col=1,cex=datacex)
      legtxt<-c(legtxt,'Data (treated)','Data (control)'); legcol<-c(legcol,1,1); legpch<-c(legpch,84,67); leglty<-c(leglty,NA,NA)
    } else {
      points(X.cts.all,Y.cts.all,type='p',pch=16,col=1,cex=datacex)
      legtxt<-c(legtxt,'Data'); legcol<-c(legcol,1); legpch<-c(legpch,16); leglty<-c(leglty,NA)
    }
  }
  for (tmpi in 1:length(rets)) {
    tmp <- ifelse(CI.interpolate,'o','p')
    if (CI.interpolate) leglty<-c(leglty,3) else leglty<-c(leglty,NA)
    points(x0s[[tmpi]],plo[[tmpi]],type=tmp,col=tmpi+1,pch=tmpi+1,lwd=CIlwd,lty=3)
    points(x0s[[tmpi]],phi[[tmpi]],type=tmp,col=tmpi+1,pch=tmpi+1,lwd=CIlwd,lty=3)
    legtxt<-c(legtxt,paste0(CI.type," ",as.character(rets[[tmpi]]$X.discrete))); legcol<-c(legcol,tmpi+1); legpch<-c(legpch,tmpi+1) 
  }
  #
  legend(legpos, legend=legtxt, inset=0.01, col=legcol, pch=legpch, lty=leglty, lwd=CIlwd, cex=1.5, y.intersp=0.9, bg='white')
}

#
# Hotelling (1939) tube uniform cv: from code for rqss.plot
#
# rqss.obj: return object from rqss()
# coverage: e.g., coverage=0.95 for 95% uniform confidence band
# newdata: vector of x0 values
rqss.unif.cv <- function(rqss.obj,coverage=NULL,newdata=NULL) {
  rdf <- rqss.obj$n - rqss.obj$edf
  V <- summary(rqss.obj, cov = TRUE)$Vqss[[1]]
  xdat <- rqss.obj$qss[[1]]$xyz[, 1]
  if (is.null(newdata)) {
    eps <- 0.01
    newdata <- seq(min(xdat) + eps, max(xdat) - eps, length=400)
  }
  G <- predict.qss1(rqss.obj$qss[[1]], data.frame(x=newdata))
  ones <- as.matrix.csr(matrix(1, nrow(G$D), 1))
  D <- cbind(ones, G$D)
#   S <- as.matrix(D %*% V %*% t(D))
#   se <- sqrt(diag(S))
  E <- eigen(as.matrix(V))
  B <- E$vectors %*% diag(sqrt(pmax(0, E$values))) %*% t(E$vectors)
  D <- as.matrix(D)
  BX1 <- B %*% t(D[-1, ])
  BX1 <- BX1/sqrt(apply(BX1^2, 2, sum))
  BX0 <- B %*% t(D[-nrow(D), ])
  BX0 <- BX0/sqrt(apply(BX0^2, 2, sum))
  kappa <- sum(sqrt(apply((BX1 - BX0)^2, 2, sum)))
  cvu <- critval(kappa, alpha = 1 - coverage, rdf = rdf)
  return(cvu)
#   #   B <- summary(rqss.obj$qss[[i]], V[[i]], ...)
#   yhats <- G$y + rqss.obj$coef["(Intercept)"]
#   return(list(cvu=cvu, x=G$x, lo=yhats-se*cvu, hi=yhats+se*cvu))
}

# From S. Weisberg, http://users.stat.umn.edu/~sandy/courses/8053/Data/mybs.R
mybs <- function(x, df, knots, degree = 3, intercept = FALSE,
    Boundary.knots = range(x), derivs) {
# This is a slight modification of the 'bs' function.  If derivs=0, the
# default, it returns the b-spline basis as usual.  If degree is a positive
# interger greater than one it returns the derivatives of the b-spline basis
# with respect to x.
   if(missing(derivs))
   { require(splines)
     mf <- match.call()
     mf[[1]] <- as.name("bs")
     return(eval(mf))
   }
   else
   {
    if(!(as.integer(derivs) %in% c(1L, 2L))) stop("derivs must be 1 or 2")
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax))
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x >
            Boundary.knots[2L])
    }
    else outside <- FALSE
    if(outside == TRUE) stop(
 "Derivatives not implemented for boundry.knots that do not cover the data")
    ord <- 1 + (degree <- as.integer(degree))
    if (ord <= 1)
        stop("'degree' must be integer >= 1")
    if (!missing(df) && missing(knots)) {
        nIknots <- df - ord + (1 - intercept)
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used  ", ord -
                (1 - intercept))
        }
        knots <- if (nIknots > 0) {
            knots <- seq.int(from = 0, to = 1, length.out = nIknots +
                2)[-c(1, nIknots + 2)]
            stats::quantile(x[!outside], knots)
        }
    }
    Aknots <- sort(c(rep(Boundary.knots, ord), knots))
    if (any(outside)) {
        warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
        derivs <- 0:degree
        scalef <- gamma(1L:ord)
        basis <- array(0, c(length(x), length(Aknots) - degree -
            1L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree,
                "^"))
            tt <- spline.des(Aknots, rep(k.pivot, ord), ord,
                derivs)$design
            basis[ol, ] <- xl %*% (tt/scalef)
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree,
                "^"))
            tt <- spline.des(Aknots, rep(k.pivot, ord), ord,
                derivs)$design
            basis[or, ] <- xr %*% (tt/scalef)
        }
        if (any(inside <- !outside))
            basis[inside, ] <- spline.des(Aknots, x[inside],
                ord)$design
    }
    else basis <- spline.des(Aknots, x, ord,
              derivs=rep(as.integer(derivs), length(x)))$design
    if (!intercept)
        basis <- basis[, -1L, drop = FALSE]
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots,
        Boundary.knots = Boundary.knots, intercept = intercept,
        derivs=as.integer(derivs))
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("mybs", "bs", "basis", "matrix")
    return(basis)
}}

predict.mybs <- function (object, newx, ...)
{
    if (missing(newx))
        return(object)
    a <- c(list(x = newx), attributes(object)[c("degree", "knots",
        "Boundary.knots", "intercept", "derivs")])
    do.call("mybs", a)
}

# Cool list[a,b,c] <- multivalued.fn(...) functionality
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
   args <- as.list(match.call())
   args <- args[-c(1:2,length(args))]
   length(value) <- length(args)
   for(i in seq(along=args)) {
     a <- args[[i]]
     if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
   }
   x
}


# EOF
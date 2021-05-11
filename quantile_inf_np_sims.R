# David M. Kaplan, KaplanDM@Missouri.edu
# First: 8 Aug 2012; updated 31 Jan 2015, Aug 2015 (QY/FL)
# Simulations for CIs for (locally smoothed) nonparametric conditional quantiles
#pt=pointwise, bonf=Bonferroni, unif=use Hotelling tube

rm(list=ls()) #clear workspace
SAVE.FLAG <- TRUE ############

####################
#                  #
# Setup/init       #
#                  #
####################
# Load packages
library(quantreg) #needed if use rqss (Koenker method)
library(splines)  #rq/bs--obviously rqss is better, but just to check...
source("quantile_inf_np.R") #for GK method

#
# rqss uniform bands: from code for rqss.plot
#
rqss.unif.band <- function(rqss.obj,coverage=NULL,newdata=NULL) {
  rdf <- rqss.obj$n - rqss.obj$edf
  V <- summary(rqss.obj, cov = TRUE)$Vqss[[1]]
#   B <- summary(rqss.obj$qss[[i]], V[[i]], ...)
  xdat <- rqss.obj$qss[[1]]$xyz[, 1]
  if (is.null(newdata)) {
    eps <- 0.01
    newdata <- seq(min(xdat) + eps, max(xdat) - eps, length=400)
  }
  G <- predict.qss1(rqss.obj$qss[[1]], data.frame(x=newdata))
  ones <- as.matrix.csr(matrix(1, nrow(G$D), 1))
  D <- cbind(ones, G$D)
  S <- as.matrix(D %*% V %*% t(D))
  se <- sqrt(diag(S))
#   cv <- qt(1 - (1 - coverage)/2, rdf)
  E <- eigen(as.matrix(V))
  B <- E$vectors %*% diag(sqrt(pmax(0, E$values))) %*% t(E$vectors)
  D <- as.matrix(D)
  BX1 <- B %*% t(D[-1, ])
  BX1 <- BX1/sqrt(apply(BX1^2, 2, sum))
  BX0 <- B %*% t(D[-nrow(D), ])
  BX0 <- BX0/sqrt(apply(BX0^2, 2, sum))
  kappa <- sum(sqrt(apply((BX1 - BX0)^2, 2, sum)))
  cvu <- critval(kappa, alpha = 1 - coverage, rdf = rdf)
#   B <- summary(rqss.obj$qss[[i]], V[[i]], ...)
  yhats <- G$y + rqss.obj$coef["(Intercept)"]
  return(list(cvu=cvu, x=G$x, lo=yhats-se*cvu, hi=yhats+se*cvu))
}

#
# Simulation and DGP parameters, flags
#
VERBOSE.FLAG <- FALSE
QY.BC <- -1 #0: no bias correction. 1: usual. -1: conservative, like QY15 paper.
if (SAVE.FLAG) {
  OUTFILE <- paste0("qinfnp_",format(Sys.time(),"%Y_%m_%d"),".txt")
  NREPLIC <- 1000 
  NUM.X0 <- 47 #47, or 231 for evaluating uniform bands
  Findices <- c(010,011,012,013,014,015,016,017,001) #for pointwise/joint; note: need 001 last b/c changes NUM.X0
  # Findices <- c(010,012,014,016,001) #for uniform bands
  FL.ORIG.FLAG <- TRUE #Use implementation from their paper's sims, or following their paper?
  QY.MEAN.FLAG <- TRUE #whether to replace estimate of 2nd derivative of conditional median w/ that of 2nd derivative of conditional mean (like Yu and Jones, 1998, p. 231)
  RQSS.FLAG <- QY.FLAG <- TRUE
  BREPLIC <- 299 #for rq
  LOCAL <- TRUE
  ALPHA <- 0.05
  n <- 400 #rqss vignette: 400
  ps <- c(0.5,0.25) #0.5 or 0.25 in paper
  # ps <- 0.5 #for uniform bands
  CUBEADJ <- n^(1/12) #to get CPE-optimal rate derived from C91
  # Others:
  #Findices <- c(010,012,014,016);  p <- 0.25
  #Findices <- c(010,012,014,016);  NUM.X0 <- 231 #0.004-spacing, to test Uniform bands
} else { #for preliminary runs
  OUTFILE <- ""
  NREPLIC <- 50
  NUM.X0 <- 47
  Findices <- c(001)
  FL.ORIG.FLAG <- TRUE #Use "original" code from their simulations, vs. what their paper describes
  QY.MEAN.FLAG <- TRUE #whether to replace estimate of 2nd derivative of conditional median w/ that of 2nd derivative of conditional mean (like Yu and Jones, 1998, p. 231)
  RQSS.FLAG <- FALSE;  QY.FLAG <- FALSE
  BREPLIC <- 0 #for rq
  LOCAL <- TRUE
  ALPHA <- 0.05
  n <- 400
  ps <- c(0.25,0.5,0.75)
  ps <- 0.5
  CUBEADJ <- n^(1/12) 
}
#
TRIM.X0 <- 0.04 #trim off top and bottom TRIM portion of X0 points
#
ONESIDED <- 0 #0:two-sided; 1/-1:upper/lower one-sided
for (Findex in Findices) { #010:017=rqss vignette variants
  for (p in ps) {

#
# Start time output
#
STARTTIME<-Sys.time()
tmp<-sprintf("Start time is %s",format(STARTTIME, "%X, %A, %d %b %Y"))
if (OUTFILE=="") {
  cat(tmp)
} else {
  cat(tmp,file=OUTFILE,sep="\n",append=TRUE)
}

###################
#                 #
# DGPs            #
#                 #
###################
set.seed(112358) #for replicability
DESC<-sprintf("Findex=%03d",Findex)
x0.ps <- TRIM.X0 + (1-2*TRIM.X0)*c(0:(NUM.X0-1))/(NUM.X0-1) #X points of interest, as quantiles of X dist'n
if (Findex==001) {
  # Fan & Liu simulation Model #1 (Table 1 in Feb. 2015 version)
  n <- 500
  Xs <- matrix(rnorm(n*NREPLIC,0,1),nrow=NREPLIC,ncol=n)
  if (NUM.X0>100) x0s <- qnorm(x0.ps) else x0s <- (0:2)*3/4 #seq(from=0.1,to=0.9,length.out=NUM.X0)
  NUM.X0 <- length(x0s)
  Y.fn <- function(X) 2.5 + sin(2*X) + 2*exp(-16*X^2)
  Us <- 0.5 * matrix(rnorm(n*NREPLIC,0,1),nrow=NREPLIC,ncol=n)
  Ys <- Y.fn(Xs) + Us
  y0s <- Y.fn(x0s) + 0.5*qnorm(p)
  #plot.quantile.inf.np(quantile.inf.np(Ys[1,],Xs[1,],x0s,p=0.5)); lines(-300:300/100,Y.fn(-300:300/100))
} else if (Findex==002) {
  #Mine (MATLAB)
  X.RAD <- 20
  f.X.fn <- function(Xs) rep.int(1/(2*X.RAD),length(Xs)) #dunif(Xs,-X.RAD,X.RAD)
  fp.X.fn <- function(Xs) rep.int(0,length(Xs))
  Xs <- replicate(n,runif(NREPLIC,-X.RAD,X.RAD)) #NREPLIC rows, n cols
  x0s <- qunif(x0.ps,-X.RAD,X.RAD)
  # x0s.unif <- qunif(c(1:n)/(n+1),-X.RAD,X.RAD)
  Y.fn <- function(X) X^2
  Y.fn.expr.txt <- expression(x^2)
  Qp.fn <- function(X) 2*X; Qpp.fn <- function(X) rep.int(2,length(X))
  Us <- replicate(n,rnorm(NREPLIC,0,1)) -qnorm(p,0,1) #p-quantile of U is zero
  Ys <- Y.fn(Xs) + Us
  y0s <- Y.fn(x0s) #since p-quantile of U is zero; was xi_p in MATLAB
  # y0s.unif <- Y.fn(x0s.unif)
  f.Y.fn <- function(p,Xs) rep.int(dnorm(qnorm(p)),length(Xs))
} else if (Findex==003) {
    # mine
    Xs <- matrix(runif(n*NREPLIC),nrow=NREPLIC,ncol=n)
    x0s <- 1:3/4;  NUM.X0 <- length(x0s)
    tmp <- 20
    Y.fn <- function(X) ifelse(abs(X-1/2)<(1/tmp), 1-tmp*abs(X-1/2), 0)
    Us <- matrix(rnorm(n*NREPLIC),nrow=NREPLIC,ncol=n)
    Ys <- Y.fn(Xs) + 0.2*Us
    y0s <- Y.fn(x0s) + 0.2*qnorm(p)
    #plot.quantile.inf.np(quantile.inf.np(Ys[1,],Xs[1,],x0s,p=0.5)); lines(-300:300/100,Y.fn(-300:300/100))
#################################################
} else if (floor(Findex/10)==01) {
  #From rqss vignette (from Ruppert, Wand, and Carroll (2003, \S17.5.1))
  Xs <- replicate(n,runif(NREPLIC,0,1)) #NREPLIC rows, n cols; X~U(0,1)
  f.X.fn <- function(Xs) rep.int(1,length(Xs))
  fp.X.fn <- function(Xs) rep.int(0,length(Xs))
  x0s <- qunif(x0.ps,0,1)
  # x0s.unif <- qunif(c(1:n)/(n+1),0,1)
  Y.fn <- function(X) (X*(1-X))^(1/2)*sin(2*pi*(1+2^(-7/5))/(X+2^(-7/5)))
  Y.fn.expr.txt <- expression(sqrt(x(1-x))*sin*bgroup("(",2*pi*(1+2^{-7/5}) / (x+2^{-7/5}),")"))
  Y.fn.expr <- quote((X*(1-X))^(1/2)*sin(2*pi*(1+2^(-7/5))/(X+2^(-7/5))))
  Qp.fn <- function(X) eval(D(Y.fn.expr,'X'))
  Qpp.fn <- function(X) eval(D(D(Y.fn.expr,'X'),'X'))
  sig0 <- 0.2
  #Homoskedastic if Findex is even, heteroskedastic if odd
  if (((Findex-10*floor(Findex/10)) %% 2) == 0) {
    sig.fn <- function(X) rep.int(sig0,length(X)) #homoskedastic
    sigp.fn <-function(X) rep.int(0,length(X)) 
    sig.fn.expr <- bquote(.(sig0)*1)
  } else {
    sig.fn <- function(X) sig0*(1+X) #heteroskedastic
    sigp.fn <-function(X) rep.int(sig0,length(X)) 
    sig.fn.expr <- bquote(.(sig0)*(1+X))
  }
  #Different distribution shapes for 010/011 vs. 012/013 vs. ...
  if (floor((Findex-10*floor(Findex/10))/2)==0) {
    Us <- replicate(n,rnorm(NREPLIC,0,1)) - qnorm(p,0,1) #Gaussian (010/011)
    f.Y.fn <- function(p,Xs) dnorm(qnorm(p,0,sig.fn(Xs)),0,sig.fn(Xs))
  } else if (floor((Findex-10*floor(Findex/10))/2)==1) {
    Us <- replicate(n,rt(NREPLIC,3)) - qt(p,3) #t3 (012/013)
    f.Y.fn <- function(p,Xs) dt(qt(p,3),3) / sig.fn(Xs)
  } else if (floor((Findex-10*floor(Findex/10))/2)==2) {
    Us <- replicate(n,rt(NREPLIC,1)) - qt(p,1) #t1 (014/015)
    f.Y.fn <- function(p,Xs) dt(qt(p,1),1) / sig.fn(Xs)
  } else if (floor((Findex-10*floor(Findex/10))/2)==3) {
    Us <- replicate(n,rchisq(NREPLIC,3)) - qchisq(p,3) #chi2(3) (016/017)
    f.Y.fn <- function(p,Xs) dchisq(qchisq(p,3),3) / sig.fn(Xs)
  } else if (floor((Findex-10*floor(Findex/10))/2)==4) {
  	U.rate <- 1/2
    Us <- replicate(n,rexp(NREPLIC,U.rate)) - qexp(p,U.rate) #exp(U.rate) (018/019)
    f.Y.fn <- function(p,Xs) dexp(qexp(p,U.rate),U.rate) / sig.fn(Xs)
  }
  Ys <- Y.fn(Xs)+sig.fn(Xs)*Us
  y0s <- Y.fn(x0s) #no add'l effect of heterosk. since p-quantile of U is zero
} else if (floor(Findex/10)==02) {
  #Like 01 but shifted right, so squiggly part not also at boundary
  TMP.SHIFT<-TRIM.X0-0.0025 #if exact, then x0s[1] has Qp=-Inf; 1/401=.0025
  Xs <- replicate(n,runif(NREPLIC,0,1)) #NREPLIC rows, n cols; X~U(0,1)
  f.X.fn <- function(Xs) rep.int(1,length(Xs))
  fp.X.fn <- function(Xs) rep.int(0,length(Xs))
  x0s <- qunif(x0.ps,0,1)
  Y.fn.orig <- function(X) (X*(1-X))^(1/2)*sin(2*pi*(1+2^(-7/5))/(X+2^(-7/5)))
  Qp.fn.orig <- function(X) eval(D(expression((X*(1-X))^(1/2)*sin(2*pi*(1+2^(-7/5))/(X+2^(-7/5)))),'X'))
  Qpp.fn.orig <- function(X) eval(D(D(expression((X*(1-X))^(1/2)*sin(2*pi*(1+2^(-7/5))/(X+2^(-7/5)))),'X'),'X'))
  Y.fn <- function(X) ifelse(X>=TMP.SHIFT, Y.fn.orig(X-TMP.SHIFT), -Y.fn.orig(TMP.SHIFT-X))
  Y.fn.expr.txt <- expression(sqrt(x(1-x))*sin*bgroup("(",2*pi*(1+2^{-7/5}) / (x+2^{-7/5}),")"))
  Qp.fn <- function(X) ifelse(X>=TMP.SHIFT, Qp.fn.orig(X-TMP.SHIFT), Qp.fn.orig(TMP.SHIFT-X))
  Qpp.fn <- function(X) ifelse(X>=TMP.SHIFT, Qpp.fn.orig(X-TMP.SHIFT), Qpp.fn.orig(TMP.SHIFT-X))
  sig0 <- 0.2
  #Homoskedastic if Findex is even, heteroskedastic if odd
  if (((Findex-10*floor(Findex/10)) %% 2) == 0) {
    sig.fn <- function(X) rep.int(sig0,length(X)) #homoskedastic
  } else {
    sig.fn <- function(X) sig0*(1+abs(X-TMP.SHIFT)) #heteroskedastic
  }
  #Different distribution shapes for 010/011 vs. 012/013 vs. ...
  if (floor((Findex-10*floor(Findex/10))/2)==0) {
    Us <- replicate(n,rnorm(NREPLIC,0,1)) - qnorm(p,0,1) #Gaussian (010/011)
    f.Y.fn <- function(p,Xs) dnorm(qnorm(p,0,sig.fn(Xs)), 0, sig.fn(Xs))
  } else if (floor((Findex-10*floor(Findex/10))/2)==1) {
    Us <- replicate(n,rt(NREPLIC,3)) - qt(p,3) #t3 (012/013)
    f.Y.fn <- function(p,Xs) dt(qt(p,3),3) / sig.fn(Xs)
  } else if (floor((Findex-10*floor(Findex/10))/2)==2) {
    Us <- replicate(n,rt(NREPLIC,1)) - qt(p,1) #t1 (014/015)
    f.Y.fn <- function(p,Xs) dt(qt(p,1),1) / sig.fn(Xs)
  } else if (floor((Findex-10*floor(Findex/10))/2)==3) {
    Us <- replicate(n,rchisq(NREPLIC,3)) - qchisq(p,3) #chi2(3) (016/017)
    f.Y.fn <- function(p,Xs) dchisq(qchisq(p,3),3) / sig.fn(Xs)
  }
  Ys <- Y.fn(Xs)+sig.fn(Xs)*Us
  y0s <- Y.fn(x0s) #no add'l effect of heterosk. since p-quantile of U is zero
} else if (floor(Findex/100)==1) {  #flat, X~N(0,1)
  #X dist'n
  f.X.fn <- function(Xs) dnorm(Xs,0,1)
  fp.X.fn <- function(Xs) -Xs*dnorm(Xs,0,1)
  Xs <- replicate(n,rnorm(NREPLIC,0,1)) #NREPLIC rows, n cols
  x0s <- qnorm(x0.ps,0,1)
  #Y.fn
  if (floor((Findex-100*floor(Findex/100))/10)==3) { #flat Y.fn
    Y.fn <- function(X) rep.int(0,length(X))
    Y.fn.expr.txt <- expression(0)
    Qp.fn <- function(X) rep.int(0,length(X))
    Qpp.fn <- function(X) rep.int(0,length(X))
  } else if (floor((Findex-100*floor(Findex/100))/10)==04) { #\Phi(x)
    Y.fn <- function(X) pnorm(X,0,1)
    Y.fn.expr.txt <- expression(Phi(x))
    Qp.fn <- function(X) dnorm(X,0,1)
    Qpp.fn <- function(X) -X*dnorm(X,0,1)
  }
  #Homoskedastic if Findex is even, heteroskedastic if odd
  sig0 <- 0.2
  if (((Findex-10*floor(Findex/10)) %% 2) == 0) {
    sig.fn <- function(X) rep.int(sig0,length(X)) #homoskedastic
  } else {
    sig.fn <- function(X) sig0*(1+abs(X)) #heteroskedastic
  }
  #Different error distribution shapes
  if (any(Findex-10*floor(Findex/10)==c(0,1))) { #Cauchy
      Us <- replicate(n,rt(NREPLIC,1)) -qt(p,1) #p-quantile of U is zero
      f.Y.fn <- function(p,Xs) rep.int(dt(qt(p,1),1),length(Xs)) / sig.fn(Xs)
  } else if (any(Findex-10*floor(Findex/10)==c(2,3))) { #Exp(1)
      Us <- replicate(n,rexp(NREPLIC,1)) -qexp(p,1)
      f.Y.fn <- function(p,Xs) rep.int(dexp(qexp(p,1),1),length(Xs)) / sig.fn(Xs)
  } else if (any(Findex-10*floor(Findex/10)==c(4,5))) { #Logn(0,1)
      Us <- replicate(n,rlnorm(NREPLIC,0,1)) -qlnorm(p,0,1)
      f.Y.fn <- function(p,Xs) rep.int(dlnorm(qlnorm(p,0,1),0,1),length(Xs)) / sig.fn(Xs)
  } else if (any(Findex-10*floor(Findex/10)==c(6,7))) { #N(0,1)
      Us <- replicate(n,rnorm(NREPLIC,0,1)) -qnorm(p,0,1)
      f.Y.fn <- function(p,Xs) rep.int(dnorm(qnorm(p,0,1),0,1),length(Xs)) / sig.fn(Xs)
  } else stop(sprintf("Not implemented: Findex=%g",Findex))
  #
  Ys <- Y.fn(Xs) + sig.fn(Xs)*Us
  y0s <- Y.fn(x0s) #since p-quantile of U is zero; was xi_p in MATLAB
}
if (NREPLIC==1) {
  Xs <- matrix(Xs,1);  Ys <- matrix(Ys,1);  Us <- matrix(Us,1)
}

ALPHA.BONF <- ALPHA/length(x0s)

#
# Filename suffix
#
FILENAME.SUFFIX <- sprintf("_F%03d_a%02g_n%d_p%02d_numx%d_trimx%02d_nrep%d",Findex,100*ALPHA,n,as.integer(round(100*p)),NUM.X0,as.integer(100*TRIM.X0),NREPLIC)



###################
#                 #
# Precomputation  #
#                 #
###################
sqrtn <- sqrt(n)
zz1 <- qnorm(1-ALPHA,0,1)
zz2 <- qnorm(1-ALPHA/2,0,1)
phizz2 <- dnorm(zz2,0,1)




###################
#                 #
# Simulation loop #
#                 #
###################
#
#Variables to hold results for rqss and npqreg CIs
#
res.rqss.ci.pt <- res.rqss.ci.bonf <- res.rqss.ci.unif <- 
  res.rq.ci.pt <- res.rq.ci.bonf <- res.rq.ci.unif <- 
  res.QYu.ci.pt <- res.QYu.ci.bonf <- res.QYu.ci.unif <- 
  res.QYg.ci.pt <- res.QYg.ci.bonf <- res.QYg.ci.unif <- 
  res.FLSu.ci.pt <- res.FLSu.ci.bonf <- res.FLSu.ci.unif <- 
  res.FLSb.ci.pt <- res.FLSb.ci.bonf <- res.FLSb.ci.unif <- 
  res.gk.ci.pt <- res.gk.ci.bonf <- res.gk.ci.unif <- matrix(NA,NREPLIC,2*NUM.X0)
res.gk.hs.pt <- res.gk.hs.bonf <- res.gk.hs.unif <- 
  res.gk.Ns.pt <- res.gk.Ns.bonf <- res.gk.Ns.unif <- matrix(NA,NREPLIC,NUM.X0) #mat.or.vec(NREPLIC,NUM.X0)
res.gk.ok.pt   <- matrix(as.logical(NA),NREPLIC,NUM.X0)
res.gk.ok.bonf <- res.gk.ok.unif <- matrix(as.logical(NA),NREPLIC,1)
Hotelling.alphas <- matrix(NA,NREPLIC,1)
#
LOOP.STARTTIME <- Sys.time()
set.seed(112358) #for replicability (bootstrap)
for (irep in 1:NREPLIC) {
  cat(sprintf("rep %d of %d (%s since loop start at %s)\r",irep,NREPLIC,capture.output(print(Sys.time()-LOOP.STARTTIME)),LOOP.STARTTIME))
  flush.console() #otherwise, sometimes won't print until loop is done...
  #Set up data for this replication
  X <- Xs[irep,]
  Y <- Ys[irep,]
  Y <- Y[order(X)];  X <- sort(X) #required later

  #
  # Koenker, rqss, following vignette (p. 10)
  #
  if (RQSS.FLAG) {
    g <- function(lam,y,x) AIC(rqss(y ~ qss(x, lambda = lam), tau=p),k = -1) #k=-1 does BIC (look @ code for AIC.rqss)
    suppressWarnings({
      lamstar <- optimize(g, interval = c(0.001, 0.5), x = X, y = Y) 
      rqss.fit <- rqss(Y~qss(X, lambda=lamstar$min), tau=p)
    })
    #
    tmp <- predict.rqss(rqss.fit,newdata=data.frame(X=x0s),interval="confidence",level=1-ALPHA)
    tmp <- t(tmp[,2:3]) #rem: read down columns, then across (need right order)
    res.rqss.ci.pt[irep,] <- tmp
    #
    tmp <- predict.rqss(rqss.fit,newdata=data.frame(X=x0s),interval="confidence",level=1-ALPHA.BONF)
    tmp <- t(tmp[,2:3]) #rem: read down columns, then across (need right order)
    res.rqss.ci.bonf[irep,] <- tmp
    #
    #compare to tmp<-plot(rqss.fit,bands="both"), tmp[[1]]$blo[100,2] and tmp[[1]]$bhi[100,2], etc.
    tmp <- rqss.unif.band(rqss.obj=rqss.fit,coverage=1-ALPHA,newdata=x0s)
    ti <- 2*(1:NUM.X0)
    res.rqss.ci.unif[irep,] <- c(tmp$lo,tmp$hi)[order(c(ti-1,ti))]
    Hotelling.alphas[irep] <- 2*(1-pnorm(tmp$cvu))
  } else {
    Hotelling.alphas[irep] <- NA
  }
  
  #
  # Goldman-Kaplan
  #
  useNA <- tryCatch({
    ret <- quantile.inf.np(Y=Y,X=X,x0s=x0s,p=p,ALPHA=ALPHA,JOINT.FLAG=TRUE,ONESIDED=ONESIDED,LOCAL=LOCAL)
    FALSE}, error=function(ERR) {
      warning(ERR$message)
      next
  })
  res.gk.ok.pt[irep,] <- ret[[1]]$pointwise.method.Lstat #If Nn too small to use Hutson, discard later
  res.gk.ci.pt[irep,] <- c(ret[[1]]$CI.pointwise.lows,ret[[1]]$CI.pointwise.highs)[order(c(2*(1:NUM.X0)-1,2*(1:NUM.X0)))]
  #Store bandwidths used
  res.gk.hs.pt[irep,] <- ret[[1]]$bandwidth.pointwise #tmpdf$pt
  res.gk.Ns.pt[irep,] <- ret[[1]]$effective.N.pointwise
  res.gk.ci.bonf[irep,] <- c(ret[[1]]$CI.joint.lows,ret[[1]]$CI.joint.highs)[order(c(2*(1:NUM.X0)-1,2*(1:NUM.X0)))]
  res.gk.hs.bonf[irep,] <- ret[[1]]$bandwidth.joint #tmpdf$bonf #hs.bonf
  res.gk.Ns.bonf[irep,] <- ret[[1]]$effective.N.joint
  res.gk.ok.bonf[irep] <- prod(ret[[1]]$joint.method.Lstat)
  if (RQSS.FLAG) { #p=1/2 => h invariant to ALPHA
    rethot <- quantile.inf.np(Y=Y,X=X,x0s=x0s,p=p,ALPHA=Hotelling.alphas[irep],JOINT.FLAG=FALSE,ONESIDED=ONESIDED,hs.pt=res.gk.hs.pt[irep,], LOCAL=LOCAL)
    res.gk.ok.unif[irep] <- prod(rethot[[1]]$pointwise.method.Lstat)
    res.gk.hs.unif[irep,] <- rethot[[1]]$bandwidth.pointwise #tmpdf$unif #hs.unif
    res.gk.Ns.unif[irep,] <- rethot[[1]]$effective.N.pointwise
    res.gk.ci.unif[irep,] <- c(rethot[[1]]$CI.pointwise.lows,rethot[[1]]$CI.pointwise.highs)[order(c(2*(1:NUM.X0)-1,2*(1:NUM.X0)))]
  }
  
  #
  # local poly rq, Chaudhuri (1991)
  #
  CVpt <- qnorm(1-ALPHA/2)
  CVbonf <- qnorm(1-ALPHA.BONF/2)
  if (RQSS.FLAG) CVunif <- qnorm(1-Hotelling.alphas[irep]/2)
  for (ix in 1:NUM.X0) {
    h <- CUBEADJ*res.gk.hs.pt[irep,ix]
    locind <- (x0s[ix]-h <= X) & (X <= x0s[ix]+h)
    Xloc <- X[locind]-x0s[ix];  Yloc <- Y[locind]
    if (length(Yloc)<4) {
      res.rq.ci.pt[irep,c(2*ix-1,2*ix)] <- res.rq.ci.bonf[irep,c(2*ix-1,2*ix)] <- res.rq.ci.unif[irep,c(2*ix-1,2*ix)] <- NA
      next
    }
    Xdes <- cbind(1,Xloc/h,Xloc^2/h^2,Xloc^3/h^3)
    useNA <- tryCatch({
      # suppressWarnings(
      rq.ret <- rq(Yloc~Xdes[,2]+Xdes[,3]+Xdes[,4],tau=p)
        # ) #never here: <<In rq.fit.sfn(x, y, tau = tau, rhs = rhs, control = control,  ... : \\ tiny diagonals replaced with Inf when calling blkfct>>
      deg <- 3
      FALSE},
      error=function(rqERR) { return(TRUE) }
    )
    if (useNA) {
      res.rq.ci.pt[irep,c(2*ix-1,2*ix)] <- res.rq.ci.bonf[irep,c(2*ix-1,2*ix)] <- res.rq.ci.unif[irep,c(2*ix-1,2*ix)] <- NA
      next
    }
    est <- rq.ret$coefficients[1]
    tryCatch({
      # suppressWarnings(
        SE <- summary(rq.ret,se='boot',R=BREPLIC)$coefficients[1,2]
        # ) #never here...: <<In rq.fit.sfn(x, y, tau = tau, rhs = rhs, control = control,  ... : \\ tiny diagonals replaced with Inf when calling blkfct>>
      res.rq.ci.pt[irep,c(2*ix-1,2*ix)] <- est + c(-1,1) * CVpt * SE
      res.rq.ci.bonf[irep,c(2*ix-1,2*ix)] <- est + c(-1,1) * CVbonf * SE
      if (RQSS.FLAG) {
        res.rq.ci.unif[irep,c(2*ix-1,2*ix)] <- est + c(-1,1) * CVunif * SE
      }
      }, error=function(ERR){
        warning(ERR$message)
        res.rq.ci.pt[irep,c(2*ix-1,2*ix)] <- res.rq.ci.bonf[irep,c(2*ix-1,2*ix)] <- res.rq.ci.unif[irep,c(2*ix-1,2*ix)] <- NA
    })
  }
  
  # Qu and Yoon (2015)
  if (QY.FLAG) {
    if (!require("np")) stop("Need package np to run Qu and Yoon (2015).")
    if (!require("fields")) stop("Need package fields to run Qu and Yoon (2015).")
    if (!require("quantreg")) stop("Need package quantreg to run Qu and Yoon (2015).")
    # estimate objects for bandwidth: fX=c(fX.c,fX.t), fYX=c(fYX.c,fYX.t), Qpp=c(Qpp.c,Qpp.t)
    # "50" means for conditional median, for use in Cor 1 optimal bandwidth
    # Otherwise, for use in estimating Thm 2 scale factor
    g50 <- function(lam,y,x) AIC(rqss(y~qss(x, lambda=lam),tau=1/2),k = -1) #k=-1 does BIC (look @ code for AIC.rqss)
    suppressWarnings({
      Q50 <- tryCatch({
        tmplamstar <- optimize(g50, interval = c(0.001, 0.5), x=X, y=Y) 
        tmp <- rqss(Y~qss(X, lambda=tmplamstar$min), tau=1/2)
        predict.rqss(tmp,newdata=data.frame(X=x0s))
      }, error=function(err) rep.int(NA,NUM.X0)) #warning=function(err) rep.int(NA,NUM.X0)
    }) #sometimes <<In rq.fit.sfn(x, y, tau = tau, rhs = rhs, control = control,  ... : \\ tiny diagonals replaced with Inf when calling blkfct>>
    #
    g <- function(lam,y,x) AIC(rqss(y~qss(x, lambda=lam),tau=p),k = -1) #k=-1 does BIC (look @ code for AIC.rqss)
    suppressWarnings({
      Q <- tryCatch({
        tmplamstar <- optimize(g, interval = c(0.001, 0.5), x=X, y=Y) 
        tmp <- rqss(Y~qss(X, lambda=tmplamstar$min), tau=p)
        predict.rqss(tmp,newdata=data.frame(X=x0s))
      }, error=function(err) rep.int(NA,NUM.X0)) #warning=function(err) rep.int(NA,NUM.X0)
    }) #sometimes <<In rq.fit.sfn(x, y, tau = tau, rhs = rhs, control = control,  ... : \\ tiny diagonals replaced with Inf when calling blkfct>>
    #
    fYX50 <- tryCatch(npcdens(npcdensbw(formula=Y~X,bwmethod='normal-reference'),txdat=X,tydat=Y,exdat=x0s,eydat=Q50)$condens,
                      warning=function(w)rep.int(NA,NUM.X0),
                      error=function(err)rep.int(NA,NUM.X0))
    fYX <- tryCatch(npcdens(npcdensbw(formula=Y~X,bwmethod='normal-reference'),txdat=X,tydat=Y,exdat=x0s,eydat=Q)$condens,
                    warning=function(w)rep.int(NA,NUM.X0),
                    error=function(err)rep.int(NA,NUM.X0))
    #
    kde.h <- hns(X,deriv.order=0)
    fX <- kde(X,eval.points=x0s,h=kde.h)$estimate
    #
    # suppressWarnings({
    tmp <- qsreg(X,Y,alpha=1/2);  Qpp50 <- predict(tmp,derivative=2,x=x0s)
    tmp <- qsreg(X,Y,alpha= p );  Qpp <- predict(tmp,derivative=2,x=x0s)
    # }) #not here either... <<In rq.fit.sfn(x, y, tau = tau, rhs = rhs, control = control,  ... : \\ tiny diagonals replaced with Inf when calling blkfct>>
    # ALTERNATIVE (QY): use Yu and Jones (1998) based on h_{mean}
    if (QY.MEAN.FLAG) { #like Yu and Jones's (1998, p. 231) argument that Qpp50 is roughly Mpp
      tmp <- smooth.spline(x=X,y=Y);  Mpp <- predict(tmp,x=x0s,deriv=2)$y
    }
    #
    # Fan & Liu w/ Yang-Stute symmetric k-NN estimator
    require(KernSmooth) #their code uses KernSmooth::dpill to get bandwidth
    FnX <- 1:n/n #since already sorted at top of NREPLIC loop
    FnX0s <- rep.int(NA,length(x0s))
    for (tmp in 1:length(x0s)) FnX0s[tmp] <- sum(X<x0s[tmp])/n
    #
    tmp <- qsreg(FnX,Y,alpha= p );  FnQpp <- predict(tmp,derivative=2,x=FnX0s)
    # kde.Fnh <- hns(FnX,deriv.order=0)
    FnfX <- rep.int(1,length(FnX0s)) #kde(FnX,eval.points=FnX0s,h=kde.Fnh)$estimate
    Fng <- function(lam,y,x) AIC(rqss(y~qss(x, lambda=lam),tau=p),k = -1) #k=-1 does BIC (look @ code for AIC.rqss)
    suppressWarnings({
      FnQ <- tryCatch({
        tmplamstar <- optimize(Fng, interval = c(0.001, 0.5), x=FnX, y=Y) 
        tmp <- rqss(Y~qss(FnX, lambda=tmplamstar$min), tau=p)
        predict.rqss(tmp,newdata=data.frame(FnX=FnX0s))
      }, error=function(err) rep.int(NA,NUM.X0)) #warning=function(err) rep.int(NA,NUM.X0)
    }) #sometimes <<In rq.fit.sfn(x, y, tau = tau, rhs = rhs, control = control,  ... : \\ tiny diagonals replaced with Inf when calling blkfct>>
    FnfYX <- tryCatch(npcdens(npcdensbw(formula=Y~FnX,bwmethod='normal-reference'),
                              txdat=FnX,tydat=Y,exdat=FnX0s,eydat=FnQ)$condens,
                      warning=function(w)rep.int(NA,NUM.X0),
                      error=function(err)rep.int(NA,NUM.X0))
    FLdelta <- 1/5+1/20
    for (ktype in 1:3) { #1=uniform, 2=Gaussian, 3=bisquare (as in FL)
      if (ktype==1) { #uniform
        Kfn <- function(u)ifelse((-1/2<=u)&(u<=1/2),1,0);  intKsq <- 1;  mu2K <- 1/12
      } else if (ktype==2) { #Gaussian
        Kfn <- dnorm;  intKsq <- 1/(2*sqrt(pi));  mu2K <- 1
      } else if (ktype==3) { #bisquare
        Kfn <- function(u) ifelse((-1<=u)&(u<=1),(15/16)*(1-u^2)^2,0) 
        intKsq <- 5/7;  mu2K <- 1/7;  intKpsq <- 15/7
      } else stop("Only ktype 1,2,3 currently supported.")
      #
      if (ktype==3) {
        UNIFcv <- (log(2)-log(abs(log(1-ALPHA))))/sqrt(2*FLdelta*log(n)) + 
          sqrt(2*FLdelta*log(n)) + log(intKpsq/(4*pi*intKsq)) / sqrt(2*FLdelta*log(n)) #FL eqn (20)
      } else FLuh.unif <- FLul.unif <- NA
      for (ix in 1:NUM.X0) {
        # bias from Theorem 1: dtau=c(dtau.c,dtau.t)
        dtau <- (1/2)*Qpp[ix]*mu2K
        #  bandwidth from (16) for median and (17) for other quantiles
        if (QY.MEAN.FLAG) {
          h50 <- n^(-1/5) * ((0.5*(1-0.5)*1*intKsq)/(fX[ix]*fYX50[ix]^2*Mpp[ix]^2*mu2K^2))^(1/5)
        } else {
          h50 <- n^(-1/5) * ((0.5*(1-0.5)*1*intKsq)/(fX[ix]*fYX50[ix]^2*Qpp50[ix]^2*mu2K^2))^(1/5)
        }
        hp <- h50 * (2*p*(1-p)/(pi*dnorm(qnorm(p))^2))^(1/5) # (17), like Yu and Jones (1998, p. 231)
        x0 <- x0s[ix]
        tmp <- tryCatch({
          kwgts <- Kfn((X-x0)/hp)
          # use bandwidths to get local linear QR estimators w/ rq()
          # suppressWarnings({
          tmp <- rq(Y~X,weights=kwgts,tau=p)
          ahat <- predict(tmp,newdata=data.frame(X=x0))
          # }) #not here either... <<In rq.fit.sfn(x, y, tau = tau, rhs = rhs, control = control,  ... : \\ tiny diagonals replaced with Inf when calling blkfct>>
          # use formulas derived from Gaussian process limits and conservative bias correction
          sig <- 1/(sqrt(n*hp*fX[ix])*fYX[ix])
          V <- p*(1-p)*intKsq*sig^2
          # CI: est +/- bias +/- z*stdev
          BC <- c(0,0)
          B <- dtau * hp^2
          if (QY.BC==1) BC <- c(B,B) else if (QY.BC==-1) BC <- c(max(0,B),min(0,B))
          c(ahat - BC[1] - qnorm(1-ALPHA/2)*sqrt(V),
            ahat - BC[2] - qnorm(ALPHA/2)*sqrt(V))
        }, error=function(err)c(NA,NA))
        CI.lo <- ifelse(is.numeric(tmp[1]),tmp[1],NA)
        CI.hi <- ifelse(is.numeric(tmp[2]),tmp[2],NA)
        if (ktype==1) {
          res.QYu.ci.pt[irep,2*ix-1] <- CI.lo
          res.QYu.ci.pt[irep,2*ix] <- CI.hi
          if (!is.na(CI.lo)) {
            res.QYu.ci.bonf[irep,2*ix-1] <- ahat - BC[1] - qnorm(1-ALPHA.BONF/2)*sqrt(V)
            res.QYu.ci.bonf[irep,2*ix] <- ahat - BC[2] - qnorm(ALPHA.BONF/2)*sqrt(V)
            if (RQSS.FLAG) {
              res.QYu.ci.unif[irep,2*ix-1] <- ahat - BC[1] - qnorm(1-Hotelling.alphas[irep]/2)*sqrt(V)
              res.QYu.ci.unif[irep,2*ix] <- ahat - BC[2] - qnorm(Hotelling.alphas[irep]/2)*sqrt(V)
            }
          }
        } else if (ktype==2) { 
          res.QYg.ci.pt[irep,2*ix-1] <- CI.lo
          res.QYg.ci.pt[irep,2*ix] <- CI.hi
          if (!is.na(CI.lo)) {
            res.QYg.ci.bonf[irep,2*ix-1] <- ahat - BC[1] - qnorm(1-ALPHA.BONF/2)*sqrt(V)
            res.QYg.ci.bonf[irep,2*ix] <- ahat - BC[2] - qnorm(ALPHA.BONF/2)*sqrt(V)
            if (RQSS.FLAG) {
              res.QYg.ci.unif[irep,2*ix-1] <- ahat - BC[1] - qnorm(1-Hotelling.alphas[irep]/2)*sqrt(V)
              res.QYg.ci.unif[irep,2*ix] <- ahat - BC[2] - qnorm(Hotelling.alphas[irep]/2)*sqrt(V)
            }
          }
        }
        
        # Fan & Liu w/ Yang-Stute symmetric k-NN estimator
        #       FnX <- 1:n/n #since already sorted at top of NREPLIC loop
        #       FnX0s <- rep.int(NA,length(x0s))
        #       for (tmp in 1:length(x0s)) FnX0s[tmp] <- sum(X<x0s[tmp])/n #can't use ix--already in ix loop
        #       # same Q and fYX as before, but not Qpp
        #       tmp <- qsreg(FnX,Y,alpha= p );  FnQpp <- predict(tmp,derivative=2,x=FnX0s)
        # Bandwidth (per FL sims) is n^(-1/20) times Yu and Jones (1998) eqn (6), which is same as Qu and Yoon (2015) eqn (16); but need to transform X by EDF first.
        YJh <- n^(-1/5) * ((p*(1-p)*intKsq)/(FnfX[ix]*FnfYX[ix]^2*FnQpp[ix]^2*mu2K^2))^(1/5)
        FLdelta <- 1/5+1/20
        FLh <- n^(1/5-FLdelta) * YJh # undersmoothing by n^(-1/20) like in Fan & Liu simulations
        FLh <- dpill(FnX,Y)*(pi/2)^(1/5)*n^(-1/20) # what they *actually* do (R code from Ruixuan Liu)
        FLsig <- sqrt(p*(1-p)*intKsq/(n*FLh)) #Fan & Liu eqn (18)
        FLul.pt <- p - qnorm(1-ALPHA/2)*FLsig
        FLuh.pt <- p + qnorm(1-ALPHA/2)*FLsig
        FLul.jt <- p - qnorm(1-ALPHA.BONF/2)*FLsig
        FLuh.jt <- p + qnorm(1-ALPHA.BONF/2)*FLsig
        if (ktype==3) { #FL eqn (21)
          FLul.unif <- p - UNIFcv*FLsig
          FLuh.unif <- p + UNIFcv*FLsig
        } else FLuh.unif <- FLul.unif <- NA
        kwgts <- Kfn((FnX-FnX0s[ix])/FLh)
        if (ktype==1) { #uniform kernel
          res.FLSu.ci.pt[irep,c(2*ix-1,2*ix)] <- tryCatch({
            if (FLul.pt<=0 || FLuh.pt>=1) stop("no")
            Yloc <- sort(Y[as.logical(kwgts)]) #for unif kernel, kwgts=1 or 0
            c(quantile.inf.interp(Yloc,FLul.pt),
              quantile.inf.interp(Yloc,FLuh.pt))
          }, error=function(err) c(NA,NA))
          res.FLSu.ci.bonf[irep,c(2*ix-1,2*ix)] <- tryCatch({
            if (FLul.jt<=0 || FLuh.jt>=1) stop("no")
            Yloc <- sort(Y[as.logical(kwgts)]) #for unif kernel, kwgts=1 or 0
            c(quantile.inf.interp(Yloc,FLul.jt),
              quantile.inf.interp(Yloc,FLuh.jt))
          }, error=function(err) c(NA,NA))
        } else if (ktype==3) {
          # done below...
        } #end of ktype==3
      }
    }
  }
  
  # Using Fan & Liu code
  source("Fan_Liu_2015_CI_code.R")
  for (ix in 1:length(x0s)) {
    # Pointwise...
    res.FLSb.ci.pt[irep,2*ix-1:0] <- 
      FL.CI.fn(X=X,Y=Y,x0=x0s[ix],tau=p,ALPHA=ALPHA,UNIF=FALSE,ORIGINAL=FL.ORIG.FLAG)
    if (any(is.na(res.FLSb.ci.pt[irep,c(2*ix-1,2*ix)]))) res.FLSb.ci.pt[irep,2*ix-1:0] <- NA
    # Joint...
    res.FLSb.ci.bonf[irep,2*ix-1:0] <- 
      FL.CI.fn(X=X,Y=Y,x0=x0s[ix],tau=p,ALPHA=ALPHA.BONF,UNIF=FALSE,ORIGINAL=FL.ORIG.FLAG)
    if (any(is.na(res.FLSb.ci.bonf[irep,2*ix-1:0]))) res.FLSb.ci.bonf[irep,2*ix-1:0] <- NA
    # Uniform...
    res.FLSb.ci.unif[irep,2*ix-1:0] <- 
      FL.CI.fn(X=X,Y=Y,x0=x0s[ix],tau=p,ALPHA=ALPHA,UNIF=TRUE,ORIGINAL=FL.ORIG.FLAG)
    if (any(is.na(res.FLSb.ci.unif[irep,2*ix-1:0]))) res.FLSb.ci.unif[irep,2*ix-1:0] <- NA
  }

} #end NREPLIC loop

# tryCatch({for (i in 1:100) {dev.off()}},error=function(DOE) {})
# dev.off() #ok to plot normally again


###################
#                 #
# Debugging       #
#                 #
###################
# No bugs detected


###################
#                 #
# Results         #
#                 #
###################
pcIncrs <- (4*0.5/sqrt(n))*c(-40:40)/10
pcAlts <- matrix(1,length(pcIncrs),1)%x%t(y0s) + matrix(1,1,length(y0s))%x%pcIncrs
#
# Compute coverage probabilities (under true H0)
#
#rqss
if (RQSS.FLAG) {
  res.rqss.pt.cover <- (res.rqss.ci.pt[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.rqss.ci.pt[,seq(2,2*length(y0s),by=2)])
  res.rqss.cp.pt <- colMeans(res.rqss.pt.cover,na.rm=TRUE)
  res.rqss.bonf.cover <- (res.rqss.ci.bonf[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.rqss.ci.bonf[,seq(2,2*length(y0s),by=2)])
  res.rqss.cp.bonf.pt <- colMeans(res.rqss.bonf.cover,na.rm=TRUE)
  res.rqss.cp.bonf <- mean(apply(res.rqss.bonf.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
  res.rqss.unif.cover <- (res.rqss.ci.unif[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.rqss.ci.unif[,seq(2,2*length(y0s),by=2)])
  res.rqss.cp.unif.pt <- colMeans(res.rqss.unif.cover,na.rm=TRUE)
  res.rqss.cp.unif <-mean(apply(res.rqss.unif.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
} else {
  res.rqss.cp.pt <- rep.int(NA,NUM.X0)
  res.rqss.cp.bonf <- res.rqss.cp.unif <- NA
}
#
#rq (local poly)
res.rq.pt.cover <- (res.rq.ci.pt[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.rq.ci.pt[,seq(2,2*length(y0s),by=2)])
res.rq.cp.pt <- colMeans(res.rq.pt.cover,na.rm=TRUE)
#
res.rq.bonf.cover <- (res.rq.ci.bonf[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.rq.ci.bonf[,seq(2,2*length(y0s),by=2)])
res.rq.cp.bonf.pt <- colMeans(res.rq.bonf.cover,na.rm=TRUE)
res.rq.cp.bonf <-mean(apply(res.rq.bonf.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
if (RQSS.FLAG) {
  res.rq.unif.cover <- (res.rq.ci.unif[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.rq.ci.unif[,seq(2,2*length(y0s),by=2)])
  res.rq.cp.unif.pt <- colMeans(res.rq.unif.cover,na.rm=TRUE)
  res.rq.cp.unif <-mean(apply(res.rq.unif.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
} else res.rq.cp.unif <- NA
#
#Qu and Yoon (2015)
res.QYu.pt.cover <- (res.QYu.ci.pt[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.QYu.ci.pt[,seq(2,2*length(y0s),by=2)])
res.QYu.cp.pt <- colMeans(res.QYu.pt.cover,na.rm=TRUE)
#
res.QYu.bonf.cover <- (res.QYu.ci.bonf[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.QYu.ci.bonf[,seq(2,2*length(y0s),by=2)])
res.QYu.cp.bonf.pt <- colMeans(res.QYu.bonf.cover,na.rm=TRUE)
res.QYu.cp.bonf <-mean(apply(res.QYu.bonf.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
if (RQSS.FLAG) {
  res.QYu.unif.cover <- (res.QYu.ci.unif[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.QYu.ci.unif[,seq(2,2*length(y0s),by=2)])
  res.QYu.cp.unif.pt <- colMeans(res.QYu.unif.cover,na.rm=TRUE)
  res.QYu.cp.unif <-mean(apply(res.QYu.unif.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
} else res.QYu.cp.unif <- NA
#
res.QYg.pt.cover <- (res.QYg.ci.pt[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.QYg.ci.pt[,seq(2,2*length(y0s),by=2)])
res.QYg.cp.pt <- colMeans(res.QYg.pt.cover,na.rm=TRUE)
#
res.QYg.bonf.cover <- (res.QYg.ci.bonf[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.QYg.ci.bonf[,seq(2,2*length(y0s),by=2)])
res.QYg.cp.bonf.pt <- colMeans(res.QYg.bonf.cover,na.rm=TRUE)
res.QYg.cp.bonf <-mean(apply(res.QYg.bonf.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
if (RQSS.FLAG) {
  res.QYg.unif.cover <- (res.QYg.ci.unif[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.QYg.ci.unif[,seq(2,2*length(y0s),by=2)])
  res.QYg.cp.unif.pt <- colMeans(res.QYg.unif.cover,na.rm=TRUE)
  res.QYg.cp.unif <-mean(apply(res.QYg.unif.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
} else res.QYg.cp.unif <- NA
#
# Fan and Liu (2013)
res.FLSu.pt.cover <- (res.FLSu.ci.pt[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.FLSu.ci.pt[,seq(2,2*length(y0s),by=2)])
res.FLSu.cp.pt <- colMeans(res.FLSu.pt.cover,na.rm=TRUE)
#
res.FLSu.bonf.cover <- (res.FLSu.ci.bonf[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.FLSu.ci.bonf[,seq(2,2*length(y0s),by=2)])
res.FLSu.cp.bonf.pt <- colMeans(res.FLSu.bonf.cover,na.rm=TRUE)
res.FLSu.cp.bonf <-mean(apply(res.FLSu.bonf.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
res.FLSu.cp.unif <- NA
#
res.FLSb.pt.cover <- (res.FLSb.ci.pt[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.FLSb.ci.pt[,seq(2,2*length(y0s),by=2)])
res.FLSb.cp.pt <- colMeans(res.FLSb.pt.cover,na.rm=TRUE)
#
res.FLSb.bonf.cover <- (res.FLSb.ci.bonf[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.FLSb.ci.bonf[,seq(2,2*length(y0s),by=2)])
res.FLSb.cp.bonf.pt <- colMeans(res.FLSb.bonf.cover,na.rm=TRUE)
res.FLSb.cp.bonf <-mean(apply(res.FLSb.bonf.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
res.FLSb.unif.cover <- (res.FLSb.ci.unif[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.FLSb.ci.unif[,seq(2,2*length(y0s),by=2)])
res.FLSb.cp.unif.pt <- colMeans(res.FLSb.unif.cover,na.rm=TRUE)
res.FLSb.cp.unif <-mean(apply(res.FLSb.unif.cover,1,prod,na.rm=FALSE),na.rm=TRUE)
#
#GK
res.gk.pt.cover <- (res.gk.ci.pt[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.gk.ci.pt[,seq(2,2*length(y0s),by=2)])
res.gk.cp.pt <- colMeans(res.gk.pt.cover, na.rm=TRUE)
res.gk.bonf.cover <- (res.gk.ci.bonf[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.gk.ci.bonf[,seq(2,2*length(y0s),by=2)])
res.gk.cp.bonf.pt <- colMeans(res.gk.bonf.cover, na.rm=TRUE)
res.gk.cp.bonf <- mean(apply(res.gk.bonf.cover,1,prod,na.rm=FALSE), na.rm=TRUE)
if (RQSS.FLAG) {
  res.gk.unif.cover <- (res.gk.ci.unif[,seq(1,2*length(y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(y0s))) & ((matrix(1,NREPLIC,1)%x%t(y0s))<res.gk.ci.unif[,seq(2,2*length(y0s),by=2)])
  res.gk.cp.unif.pt <- colMeans(res.gk.unif.cover, na.rm=TRUE)
  res.gk.cp.unif <- mean(apply(res.gk.unif.cover,1,prod,na.rm=FALSE), na.rm=TRUE)
} else res.gk.cp.unif <- NA

#
# Compute exclusion probabilities (under H1s)
#
#rqss
res.rqss.ep.alts <- matrix(NA,length(pcIncrs),2+length(y0s))
if (RQSS.FLAG) {
  for (i in 1:length(pcIncrs)) {
    tmp.y0s <- pcAlts[i,]
    tmp.pt.cover <- (res.rqss.ci.pt[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.rqss.ci.pt[,seq(2,2*length(tmp.y0s),by=2)])
    res.rqss.ep.alts[i,-(1:2)] <- 1-colMeans(tmp.pt.cover, na.rm=TRUE)
    tmp.bonf.cover <- (res.rqss.ci.bonf[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.rqss.ci.bonf[,seq(2,2*length(tmp.y0s),by=2)])
    res.rqss.ep.alts[i, 2] <- 1-mean(apply(tmp.bonf.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
    tmp.unif.cover <- (res.rqss.ci.unif[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.rqss.ci.unif[,seq(2,2*length(tmp.y0s),by=2)])
    res.rqss.ep.alts[i, 1] <- 1-mean(apply(tmp.unif.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
  }
}
#
#rq (local cubic)
res.rq.ep.alts <- matrix(NA,length(pcIncrs),2+length(y0s))
for (i in 1:length(pcIncrs)) {
  tmp.y0s <- pcAlts[i,]
  tmp.pt.cover <- (res.rq.ci.pt[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.rq.ci.pt[,seq(2,2*length(tmp.y0s),by=2)])
  res.rq.ep.alts[i,-(1:2)] <- 1-colMeans(tmp.pt.cover, na.rm=TRUE)
  tmp.bonf.cover <- (res.rq.ci.bonf[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.rq.ci.bonf[,seq(2,2*length(tmp.y0s),by=2)])
  res.rq.ep.alts[i, 2] <- 1-mean(apply(tmp.bonf.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
  if (RQSS.FLAG) {
    tmp.unif.cover <- (res.rq.ci.unif[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.rq.ci.unif[,seq(2,2*length(tmp.y0s),by=2)])
    res.rq.ep.alts[i, 1] <- 1-mean(apply(tmp.unif.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
  }
}
#
# Qu and Yoon (2015)
res.QYu.ep.alts <- matrix(NA,length(pcIncrs),2+length(y0s))
for (i in 1:length(pcIncrs)) {
  tmp.y0s <- pcAlts[i,]
  tmp.pt.cover <- (res.QYu.ci.pt[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.QYu.ci.pt[,seq(2,2*length(tmp.y0s),by=2)])
  res.QYu.ep.alts[i,-(1:2)] <- 1-colMeans(tmp.pt.cover, na.rm=TRUE)
  tmp.bonf.cover <- (res.QYu.ci.bonf[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.QYu.ci.bonf[,seq(2,2*length(tmp.y0s),by=2)])
  res.QYu.ep.alts[i, 2] <- 1-mean(apply(tmp.bonf.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
  if (RQSS.FLAG) {
    tmp.unif.cover <- (res.QYu.ci.unif[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.QYu.ci.unif[,seq(2,2*length(tmp.y0s),by=2)])
    res.QYu.ep.alts[i, 1] <- 1-mean(apply(tmp.unif.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
  }
}
#
res.QYg.ep.alts <- matrix(NA,length(pcIncrs),2+length(y0s))
for (i in 1:length(pcIncrs)) {
  tmp.y0s <- pcAlts[i,]
  tmp.pt.cover <- (res.QYg.ci.pt[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.QYg.ci.pt[,seq(2,2*length(tmp.y0s),by=2)])
  res.QYg.ep.alts[i,-(1:2)] <- 1-colMeans(tmp.pt.cover, na.rm=TRUE)
  tmp.bonf.cover <- (res.QYg.ci.bonf[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.QYg.ci.bonf[,seq(2,2*length(tmp.y0s),by=2)])
  res.QYg.ep.alts[i, 2] <- 1-mean(apply(tmp.bonf.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
  if (RQSS.FLAG) {
    tmp.unif.cover <- (res.QYg.ci.unif[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.QYg.ci.unif[,seq(2,2*length(tmp.y0s),by=2)])
    res.QYg.ep.alts[i, 1] <- 1-mean(apply(tmp.unif.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
  }
}
#
# Fan and Liu (2013)
res.FLSu.ep.alts <- matrix(NA,length(pcIncrs),2+length(y0s))
for (i in 1:length(pcIncrs)) {
  tmp.y0s <- pcAlts[i,]
  tmp.pt.cover <- (res.FLSu.ci.pt[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.FLSu.ci.pt[,seq(2,2*length(tmp.y0s),by=2)])
  res.FLSu.ep.alts[i,-(1:2)] <- 1-colMeans(tmp.pt.cover, na.rm=TRUE)
  tmp.bonf.cover <- (res.FLSu.ci.bonf[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.FLSu.ci.bonf[,seq(2,2*length(tmp.y0s),by=2)])
  res.FLSu.ep.alts[i, 2] <- 1-mean(apply(tmp.bonf.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
}
#
res.FLSb.ep.alts <- matrix(NA,length(pcIncrs),2+length(y0s))
for (i in 1:length(pcIncrs)) {
  tmp.y0s <- pcAlts[i,]
  tmp.pt.cover <- (res.FLSb.ci.pt[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.FLSb.ci.pt[,seq(2,2*length(tmp.y0s),by=2)])
  res.FLSb.ep.alts[i,-(1:2)] <- 1-colMeans(tmp.pt.cover, na.rm=TRUE)
  tmp.bonf.cover <- (res.FLSb.ci.bonf[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.FLSb.ci.bonf[,seq(2,2*length(tmp.y0s),by=2)])
  res.FLSb.ep.alts[i, 2] <- 1-mean(apply(tmp.bonf.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
  tmp.unif.cover <- (res.FLSb.ci.unif[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.FLSb.ci.unif[,seq(2,2*length(tmp.y0s),by=2)])
  res.FLSb.ep.alts[i, 1] <- 1-mean(apply(tmp.unif.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
}
#
#gk
res.gk.ep.alts<-matrix(NA,length(pcIncrs),2+length(y0s))
for (i in 1:length(pcIncrs)) {
    tmp.y0s <- pcAlts[i,]
    tmp.pt.cover <- (res.gk.ci.pt[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.gk.ci.pt[,seq(2,2*length(tmp.y0s),by=2)])
    res.gk.ep.alts[i,-(1:2)] <- 1-colMeans(tmp.pt.cover, na.rm=TRUE)
    tmp.bonf.cover <- (res.gk.ci.bonf[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.gk.ci.bonf[,seq(2,2*length(tmp.y0s),by=2)])
    res.gk.ep.alts[i, 2] <- 1-mean(apply(tmp.bonf.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
    if (RQSS.FLAG) {
      tmp.unif.cover <- (res.gk.ci.unif[,seq(1,2*length(tmp.y0s),by=2)]<(matrix(1,NREPLIC,1)%x%t(tmp.y0s))) & ((matrix(1,NREPLIC,1)%x%t(tmp.y0s))<res.gk.ci.unif[,seq(2,2*length(tmp.y0s),by=2)])
      res.gk.ep.alts[i, 1] <- 1-mean(apply(tmp.unif.cover,1,prod, na.rm=FALSE),na.rm=TRUE)
    }
}

#
# Format coverage probs to file
#
tmp0 <- ">>> Coverage Probabilities <<<"
#
tmp1 <- paste0(sprintf("%-11s & Unif  & Joint & %s","MethodName",paste0(sprintf('&x0=%6.3f',x0s),collapse=""))," \\\\")
#
tmp2a <- sprintf("%-11s & %5.3f & %5.3f &", "rqss", res.rqss.cp.unif, res.rqss.cp.bonf)
tmp2b <- paste0(sprintf(" & %7.3f",res.rqss.cp.pt),collapse="")
tmp2 <- paste0(tmp2a,tmp2b,"  \\\\")
#
tmp3a <- sprintf("%-11s & %5.3f & %5.3f &", "L-stat", res.gk.cp.unif, res.gk.cp.bonf)
tmp3b <- paste0(sprintf(" & %7.3f",res.gk.cp.pt),collapse="")
tmp3 <- paste0(tmp3a,tmp3b,"  \\\\")
#
tmp4a <- sprintf("%-11s & %5.3f & %5.3f &", "rq", res.rq.cp.unif, res.rq.cp.bonf)
tmp4b <- paste0(sprintf(" & %7.3f",res.rq.cp.pt),collapse="")
tmp4 <- paste0(tmp4a,tmp4b,"  \\\\")
#
tmp5a <- sprintf("%-11s & %5.3f & %5.3f &", "QYu", res.QYu.cp.unif, res.QYu.cp.bonf)
tmp5b <- paste0(sprintf(" & %7.3f",res.QYu.cp.pt),collapse="")
tmp5  <- paste0(tmp5a,tmp5b,"  \\\\")
#
tmp6a <- sprintf("%-11s & %5.3f & %5.3f &", "QYg", res.QYg.cp.unif, res.QYg.cp.bonf)
tmp6b <- paste0(sprintf(" & %7.3f",res.QYg.cp.pt),collapse="")
tmp6  <- paste0(tmp6a,tmp6b,"  \\\\")
#
tmp7a <- sprintf("%-11s & %5.3f & %5.3f &", "FLSu", res.FLSu.cp.unif, res.FLSu.cp.bonf)
tmp7b <- paste0(sprintf(" & %7.3f",res.FLSu.cp.pt),collapse="")
tmp7  <- paste0(tmp7a,tmp7b,"  \\\\")
#
tmp8a <- sprintf("%-11s & %5.3f & %5.3f &", "FLSb", res.FLSb.cp.unif, res.FLSb.cp.bonf)
tmp8b <- paste0(sprintf(" & %7.3f",res.FLSb.cp.pt),collapse="")
tmp8  <- paste0(tmp8a,tmp8b,"  \\\\")
#
tmp9 <- ""
tmp.all <- sprintf("%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n",tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9)
cp.table <- tmp.all

#
# Format joint exclusion probabilities (under H1s) to file
#
tmp0 <- ">>> Joint Power <<<\n"
tmp1 <- "&\\multicolumn{9}{c}{Deviation from null at all points}\\\\\n"
pcSeq <- seq(from=21,by=5,length.out=9) #max(which(pcIncrs<=-0.2))
tmp2 <- paste0("Method Name",paste0(sprintf(' & %6.3f',pcIncrs[pcSeq]),collapse="")," \\\\\n")
#
tmp3 <- paste0(" GK (Bonf) ",paste0(sprintf(' & %6.3f',res.gk.ep.alts[pcSeq,2]),collapse="")," \\\\\n")
tmp4 <- paste0("GK (Hotel.)",paste0(sprintf(' & %6.3f',res.gk.ep.alts[pcSeq,1]),collapse="")," \\\\\n")
tmp5 <- paste0("rqss (Bonf)",paste0(sprintf(' & %6.3f',res.rqss.ep.alts[pcSeq,2]),collapse="")," \\\\\n")
tmp6 <- paste0("rqss (unif)",paste0(sprintf(' & %6.3f',res.rqss.ep.alts[pcSeq,1]),collapse="")," \\\\\n")
tmp7 <- paste0(" rq (Bonf) ",paste0(sprintf(' & %6.3f',res.rq.ep.alts[pcSeq,2]),collapse="")," \\\\\n")
tmp8 <- paste0(" rq (Hot.) ",paste0(sprintf(' & %6.3f',res.rq.ep.alts[pcSeq,1]),collapse="")," \\\\\n")
tmp9 <- paste0(" QYu (Bonf) ",paste0(sprintf(' & %6.3f',res.QYu.ep.alts[pcSeq,2]),collapse="")," \\\\\n")
tmp10 <- paste0(" QYu (Hot.) ",paste0(sprintf(' & %6.3f',res.QYu.ep.alts[pcSeq,1]),collapse="")," \\\\\n")
tmp11 <- paste0(" QYg (Bonf) ",paste0(sprintf(' & %6.3f',res.QYg.ep.alts[pcSeq,2]),collapse="")," \\\\\n")
tmp12 <- paste0(" QYg (Hot.) ",paste0(sprintf(' & %6.3f',res.QYg.ep.alts[pcSeq,1]),collapse="")," \\\\\n")
tmp13 <- paste0(" FLSu (Bonf) ",paste0(sprintf(' & %6.3f',res.FLSu.ep.alts[pcSeq,2]),collapse="")," \\\\\n")
tmp14 <- paste0(" FLSb (Bonf) ",paste0(sprintf(' & %6.3f',res.FLSb.ep.alts[pcSeq,2]),collapse="")," \\\\\n")
tmp15 <- paste0(" FLSb (Band) ",paste0(sprintf(' & %6.3f',res.FLSb.ep.alts[pcSeq,1]),collapse="")," \\\\\n")
pwr.table <- c(tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15)

#
# Output to console/file
#
desc.str <- sprintf("%s,n=%d,p=%5.3f,ALPHA=%5.3f,NREPLIC=%d,BREPLIC=%d\nx0.ps=[%s]",DESC,n,p,ALPHA,NREPLIC,BREPLIC,paste0(sprintf("%5.3f",x0.ps),collapse=","))
if (OUTFILE=="") {
  cat(desc.str,sep='\n')
  tmp <- quantile(Hotelling.alphas,c(0.25,0.50,0.75),na.rm=TRUE)
  cat(sprintf("Hotelling.alphas quartiles: %g, %g, %g",tmp[1],tmp[2],tmp[3]),sep='\n')
  if (TRUE) {
    cat(cp.table)
    cat(pwr.table)
  }
} else {
  cat(desc.str,file=OUTFILE,sep="\n",append=TRUE)
  tmp <- quantile(Hotelling.alphas,c(0.25,0.50,0.75))
  cat(sprintf("Hotelling.alphas quartiles: %g, %g, %g",tmp[1],tmp[2],tmp[3]),
      sep='\n',file=OUTFILE,append=TRUE)
  if (TRUE) {
    cat(cp.table,file=OUTFILE,sep="",append=TRUE)  #already have \n
    cat(pwr.table,file=OUTFILE,sep="",append=TRUE) #already have \n
  }
}
# 

#
# h,N
#
med.h.pt <- apply(res.gk.hs.pt,2,median)
tmp.out <- paste0("gk med. hs: ",paste0(sprintf(" & %5.3f",med.h.pt),collapse=""))
if (OUTFILE=="") {
  cat(tmp.out,sep="\n")
} else {
  cat(tmp.out,file=OUTFILE,sep="\n",append=TRUE)
}
#
iqr.h.pt <- apply(res.gk.hs.pt,2,IQR)
tmp.out <- paste0("gk IQR(hs): ",paste0(sprintf(" & %5.3f",iqr.h.pt),collapse=""))
if (OUTFILE=="") {
  cat(tmp.out,sep="\n")
} else {
  cat(tmp.out,file=OUTFILE,sep="\n",append=TRUE)
}
#
med.N.pt <- apply(res.gk.Ns.pt,2,median)
tmp.out <- paste0("gk med. Ns: ",paste0(sprintf(" & %g",med.N.pt),collapse=""))
if (OUTFILE=="") {
  cat(tmp.out,sep="\n")
} else {
  cat(tmp.out,file=OUTFILE,sep="\n",append=TRUE)
}


###################
#                 #
# Graphs          #
#                 #
###################

#X prep
tmp.xr.full <- range(Xs)
tmp.xx.full <- seq(from=tmp.xr.full[1],to=tmp.xr.full[2],by=(tmp.xr.full[2]-tmp.xr.full[1])/n)
#
tmp.xr <- range(x0s)
tmp.xx <- seq(from=tmp.xr[1],to=tmp.xr[2],by=(tmp.xr[2]-tmp.xr[1])/n)



############# GRAPHS FOR PAPER ###################

LWD <- 4

COLS <- c("#A6BDDB","#67A9CF","#3690C0","#02818A","#016450") #from RColorBrewer::brewer.pal(7,"PuBuGn")[-(1:2)]
COLS <- c("#A6BDDB","#3690C0","#02818A","#016C59","#014636") #from RColorBrewer::brewer.pal(7,"PuBuGn")[c(4,6:9)]
COLS <- rev(COLS)

#
# True function
#
if (SAVE.FLAG && FALSE) {
  pdf(file=sprintf("F%03d_true_fn.pdf",Findex),pointsize=12, width=10, height=7)
  par(family="serif",mar=c(5.0,6.0,6.0,2.1))
  plot(tmp.xx.full,Y.fn(tmp.xx.full),type="l",lwd=4,xlab="X",ylab=expression(Q[Y*symbol("|")*X](p*symbol(";")*x)),xlim=c(0,1), mgp=c(3,1,0),col=1,lty=1,main=Y.fn.expr.txt, cex.main=2,cex.lab=2,cex.axis=2)
  points(x0s,Y.fn(x0s),pch=1,cex=3)
  dev.off()
}

#
# Graph: Pointwise CP (only)
#
if (SAVE.FLAG) pdf(file=sub(".txt",paste0("_ptCP",FILENAME.SUFFIX,".pdf"),OUTFILE)) #else draw to current device
par(family="serif",mar=c(5.0,6.0,6.0,2.1))
plot(x=x0s, y=rep.int(100*(1-ALPHA),NUM.X0), col=1,lwd=2,lty=1,type="l", ylim=c(50,100), cex.main=2.0, main="Pointwise Coverage Probability", xlab="X", ylab="CP (%)", mgp=c(3.0,1.0,0), cex.axis=1.5, cex.lab=2.0)
  lines(x0s,100*res.FLSb.cp.pt,lwd=LWD,lty=6,col=COLS[5],pch=NA)
  lines(x0s, 100*res.QYg.cp.pt,lwd=LWD,lty=5,col=COLS[4],pch=NA)
  lines(x0s,  100*res.rq.cp.pt,lwd=LWD,lty=4,col=COLS[3],pch=NA) #was:col=4
  lines(x0s,100*res.rqss.cp.pt,lwd=LWD,lty=3,col=COLS[2],pch=NA) #was:col=3
  lines(x0s,100*res.gk.cp.pt,lwd=LWD,lty=1,col=COLS[1],pch=NA) #was:col=2
  # lines(x0s,100*res.gkbc.cp.pt,lwd=LWD,lty=3,col=2,pch=NA)
  #
  tmp<-c(expression(1-alpha),"L-stat","rqss","boot","QYg","FLb")
  if (Findex==001) {
    legend("bottomright", tmp, inset=c(0.01,0.00), bty='n',
           col=c("#000000",COLS[1:5]), pch=NA, lty=c(1,1,3:6), lwd=c(2,LWD,rep(LWD,4)), 
           cex=1.8, x.intersp=0.5, y.intersp=1.0, ncol=2, text.width=rep(c(0.2,0.15),each=3))
  } else {
      legend("bottomright", tmp, inset=c(0.01,0.00), bty='n',
             col=c("#000000",COLS[1:5]), pch=NA, lty=c(1,1,3:6), lwd=c(2,LWD,rep(LWD,4)), 
             cex=1.8, x.intersp=0.5, y.intersp=1.0, ncol=2, text.width=rep(c(0.13,0.1),each=3))
  }
if (SAVE.FLAG) dev.off()

# #
# # Censored version of ptwise pwr: need at least 1-2*ALPHA CP
# #
# if (FALSE) {
# if (SAVE.FLAG) pdf(file=sub(".txt",paste0("_ptPWRcens",FILENAME.SUFFIX,".pdf"),OUTFILE))
# ind1 <- which(pcIncrs==-0.1)
# ind2 <- which(pcIncrs==0.1)
# par(family="serif",mar=c(5.0,6.0,6.0,2.1))
# THRESH<-1-2*ALPHA
# plot(x=x0s, y=rep.int(1,NUM.X0), col=1,lwd=4,lty=1,type="n", ylim=c(0,100), cex.main=2.0, main=sprintf("Pointwise Power (Omitting CP<%g%%)",100*THRESH),sub=sprintf("Deviation: %g or %g",pcIncrs[ind1],pcIncrs[ind2]), xlab="X", ylab="Rejection Probability (%)", mgp=c(3.0,1.0,0), cex.axis=1.5, cex.lab=2.0, cex.sub=1.0)
# #
# tmp1 <- ifelse(res.gk.cp.pt<THRESH,NA,100*res.gk.ep.alts[ind1,3:dim(res.gk.ep.alts)[2]])
# tmp2 <- ifelse(res.gk.cp.pt<THRESH,NA,100*res.gk.ep.alts[ind2,3:dim(res.gk.ep.alts)[2]])
# tmp <- (tmp1+tmp2)/2 #NA if either is NA
# lines(type="l", x0s,tmp, lwd=3,lty=1,col=2,pch=1)
# lines(type="p", x0s,tmp, lwd=1,lty=1,col=2,pch=1,cex=2)
# #
# tmp1 <- ifelse(res.rqss.cp.pt<THRESH,NA,100*res.rqss.ep.alts[ind1,3:dim(res.rqss.ep.alts)[2]])
# tmp2 <- ifelse(res.rqss.cp.pt<THRESH,NA,100*res.rqss.ep.alts[ind2,3:dim(res.rqss.ep.alts)[2]])
# tmp <- (tmp1+tmp2)/2 #NA if either is NA
# lines(type="l", x0s,tmp, lwd=3,lty=2,col=3,pch=2)
# lines(type="p", x0s,tmp, lwd=1,lty=2,col=3,pch=2,cex=2)
# #
# tmp1 <- ifelse(res.rq.cp.pt<THRESH,NA,100*res.rq.ep.alts[ind1,3:dim(res.rq.ep.alts)[2]])
# tmp2 <- ifelse(res.rq.cp.pt<THRESH,NA,100*res.rq.ep.alts[ind2,3:dim(res.rq.ep.alts)[2]])
# tmp <- (tmp1+tmp2)/2 #NA if either is NA
# lines(type="l", x0s,tmp, lwd=3,lty=4,col=4,pch=6)
# lines(type="p", x0s,tmp, lwd=1,lty=4,col=4,pch=6,cex=2)
# #
# legend("top", c("L-stat","rqss","local cubic"), inset=0.0, bty='n', col=c(2,3,4), pch=c(1,2,6), lty=c(1,2,4), lwd=2, x.intersp=0.5, y.intersp=0.8, cex=1.8, ncol=3,
#        text.width=c(0.18,0.2,0.1))
# if (SAVE.FLAG) dev.off()
# }

#
# SECOND version of ptwise pwr: uncensored
#
if (SAVE.FLAG) pdf(file=sub(".txt",paste0("_ptPWR2",FILENAME.SUFFIX,".pdf"),OUTFILE))
ind1 <- which.min(abs(pcIncrs - (-0.1)))
ind2 <- which.min(abs(pcIncrs-0.1))
par(family="serif",mar=c(5.0,6.0,6.0,2.1))
plot(x=x0s, y=rep.int(1,NUM.X0), col=1,lwd=4,lty=1,type="n", ylim=c(0,100), cex.main=2.0, main="Pointwise Power",sub=sprintf("Deviation: %g or %g",pcIncrs[ind1],pcIncrs[ind2]), xlab="X", ylab="Rejection Probability (%)", mgp=c(3.0,1.0,0), cex.axis=1.5, cex.lab=2.0, cex.sub=1.0)
#
tmp <- (100*res.FLSb.ep.alts[ind1,3:dim(res.FLSb.ep.alts)[2]]+100*res.FLSb.ep.alts[ind2,3:dim(res.FLSb.ep.alts)[2]])/2 #avg. of -0.1 and +0.1 pwr
lines(x0s,tmp, type="l", lwd=LWD,lty=6,col=COLS[5],pch=6)
#
tmp <- (100*res.QYg.ep.alts[ind1,3:dim(res.QYg.ep.alts)[2]]+100*res.QYg.ep.alts[ind2,3:dim(res.QYg.ep.alts)[2]])/2 #avg. of -0.1 and +0.1 pwr
lines(x0s,tmp, type="l", lwd=LWD,lty=5,col=COLS[4],pch=2)
#
tmp <- (100*res.rq.ep.alts[ind1,3:dim(res.rq.ep.alts)[2]]+100*res.rq.ep.alts[ind2,3:dim(res.rq.ep.alts)[2]])/2 #avg. of -0.1 and +0.1 pwr
lines(x0s,tmp, type="l", lwd=LWD,lty=4,col=COLS[3],pch=6)
#
tmp <- (100*res.rqss.ep.alts[ind1,3:dim(res.rqss.ep.alts)[2]]+100*res.rqss.ep.alts[ind2,3:dim(res.rqss.ep.alts)[2]])/2 #avg. of -0.1 and +0.1 pwr
lines(x0s,tmp, type="l", lwd=LWD,lty=3,col=COLS[2],pch=2)
#
tmp <- (100*res.gk.ep.alts[ind1,3:dim(res.gk.ep.alts)[2]]+100*res.gk.ep.alts[ind2,3:dim(res.gk.ep.alts)[2]])/2 #avg. of +0.1 and -0.1 pwr
lines(x0s,tmp, type="l", lwd=LWD+1,lty=1,col=COLS[1],pch=1)
#
if (Findex==001) {
  legend("top", c("L-stat","rqss","boot","QYg","FLb",""), inset=0.01, bty='n',
         col=COLS[1:5], pch=NA, lty=c(1,3:6,NA), lwd=LWD,
         x.intersp=0.5, y.intersp=1.0, cex=1.8, ncol=3,
         text.width=rep(c(0.225,0.15,0.375),each=2))
} else {
  legend("top", c("L-stat","rqss","boot","QYg","FLb",""), inset=0.01, bty='n',
         col=COLS[1:5], pch=NA, lty=c(1,3:6,NA), lwd=LWD,
         x.intersp=0.5, y.intersp=1.0, cex=1.8, ncol=3,
         text.width=rep(c(0.15,0.1,0.25),each=2))
}
if (SAVE.FLAG) dev.off()

#
# Graph of joint power curves
#
if (SAVE.FLAG) pdf(file=sub(".txt",paste0("_jtPWR",FILENAME.SUFFIX,".pdf"),OUTFILE))
pcSeq <- 11:71 #pretty much all 100% outside this for n=400
par(family="serif",mar=c(5.0,6.0,6.0,2.1))
plot(x=pcIncrs[pcSeq], y=rep.int(100*ALPHA,length(pcIncrs[pcSeq])), col=1,lwd=2,lty=1,type="l", ylim=c(0,100), cex.main=2.0, main="Joint Power Curves", xlab="Deviation of Null from Truth", ylab="Rejection Probability (%)", mgp=c(3.0,1.0,0), cex.lab=2, cex.axis=1.5) #, sub="Entire Curve Shifted Vertically By Deviation"
#
lines(pcIncrs[pcSeq],100*res.FLSb.ep.alts[pcSeq,2],lwd=LWD,lty=6,col=COLS[5],pch=NA)
lines(pcIncrs[pcSeq], 100*res.QYg.ep.alts[pcSeq,2],lwd=LWD,lty=5,col=COLS[4],pch=NA)
lines(pcIncrs[pcSeq],  100*res.rq.ep.alts[pcSeq,2],lwd=LWD,lty=4,col=COLS[3],pch=NA)
lines(pcIncrs[pcSeq],100*res.rqss.ep.alts[pcSeq,2],lwd=LWD,lty=3,col=COLS[2],pch=NA)
lines(pcIncrs[pcSeq],100*res.gk.ep.alts[pcSeq,2],lwd=LWD,lty=1,col=COLS[1],pch=NA) #2=Bonferroni; 1=unif
#
tmp <- c(expression(alpha),"L-stat","rqss","boot","QYg","FLb")
# legend("bottomright", tmp, inset=c(0.00,0.06), bty='n',
#        col=c("#000000",COLS[1:5]), pch=NA, lty=c(1,1,3:6), lwd=c(2,LWD,rep(LWD,4)), cex=1.8, x.intersp=0.5, y.intersp=0.8) #box.lwd=0, box.col='white',bg='white',
if (Findex==001) {
  legend("top", tmp, inset=0.01, bty='n',
         col=c("#000000",COLS[1:5]), pch=NA, lty=c(1,1,3:6), lwd=c(1,rep(LWD,5)),
         x.intersp=0.5, y.intersp=1.0, cex=1.8, ncol=3,
         text.width=rep(c(0.1,0.06,0.15),each=2))
} else {
  legend("bottomleft", tmp[1:3], inset=c(0.02,0.06), bty='n',
         col=c("#000000",COLS[1:2]), pch=NA, lty=c(1,1,3), lwd=c(2,LWD,LWD), cex=1.8, x.intersp=0.5, y.intersp=1.0) #box.lwd=0, box.col='white',bg='white',
  legend("bottomright", tmp[-(1:3)], inset=c(0.02,0.06), bty='n',
         col=c(COLS[3:5]), pch=NA, lty=c(4:6), lwd=c(rep(LWD,3)), cex=1.8, x.intersp=0.5, y.intersp=1.0) #box.lwd=0, box.col='white',bg='white',
}
if (SAVE.FLAG) dev.off()

#
# Graph of uniform power curves
#
if (RQSS.FLAG) {
  if (SAVE.FLAG) pdf(file=sub(".txt",paste0("_unifPWR",FILENAME.SUFFIX,".pdf"),OUTFILE))
  pcSeq <- 11:71 #pretty much all 100% outside this for n=400
  par(family="serif",mar=c(5.0,6.0,6.0,2.1))
  plot(x=pcIncrs[pcSeq], y=rep.int(100*ALPHA,length(pcIncrs[pcSeq])), col=1,lwd=2,lty=1,type="l", ylim=c(0,100), cex.main=2.0, main="Uniform Power Curves", xlab="Deviation of Null from Truth", ylab="Rejection Probability (%)", mgp=c(3.0,1.0,0), cex.lab=2, cex.axis=1.5) #, sub="Entire Curve Shifted Vertically By Deviation"
  #
  lines(pcIncrs[pcSeq],100*res.FLSb.ep.alts[pcSeq,1],lwd=LWD,lty=6,col=COLS[5],pch=NA)
  lines(pcIncrs[pcSeq], 100*res.QYg.ep.alts[pcSeq,1],lwd=LWD,lty=5,col=COLS[4],pch=NA)
  lines(pcIncrs[pcSeq],  100*res.rq.ep.alts[pcSeq,1],lwd=LWD,lty=4,col=COLS[3],pch=NA)
  lines(pcIncrs[pcSeq],100*res.rqss.ep.alts[pcSeq,1],lwd=LWD,lty=3,col=COLS[2],pch=NA)
  lines(pcIncrs[pcSeq],100*res.gk.ep.alts[pcSeq,1],lwd=LWD,lty=1,col=COLS[1],pch=NA) #2=Bonferroni; 1=unif
  #
  tmp <- c(expression(alpha),"L-stat","rqss","boot","QYg","FLb") 
#   legend("bottomright", tmp, inset=c(0.00,0.06), bty='n',
#          col=c("#000000",COLS[1:5]), pch=NA, lty=c(1,1,3:6), lwd=c(2,LWD,rep(LWD,4)), cex=1.8, x.intersp=0.5, y.intersp=0.8) #box.lwd=0, box.col='white',bg='white',
  if (Findex==001) {
    legend("top", tmp, inset=0.01, bty='n',
           col=c("#000000",COLS[1:5]), pch=NA, lty=c(1,1,3:6), lwd=c(1,rep(LWD,5)),
           x.intersp=0.5, y.intersp=1.0, cex=1.8, ncol=3,
           text.width=rep(c(0.1,0.06,0.15),each=2))
  } else {
    legend("bottomleft", tmp[1:3], inset=c(0.02,0.06), bty='n',
           col=c("#000000",COLS[1:2]), pch=NA, lty=c(1,1,3), lwd=c(2,LWD,LWD), cex=1.8, x.intersp=0.5, y.intersp=1.0) #box.lwd=0, box.col='white',bg='white',
    legend("bottomright", tmp[-(1:3)], inset=c(0.02,0.06), bty='n',
           col=c(COLS[3:5]), pch=NA, lty=c(4:6), lwd=c(rep(LWD,3)), cex=1.8, x.intersp=0.5, y.intersp=1.0) #box.lwd=0, box.col='white',bg='white',
  }
  if (SAVE.FLAG) dev.off()
}




###################
#                 #
# Cleanup         #
#                 #
###################
tmpt <- as.numeric(Sys.time()-STARTTIME,units="secs")
tmps <- sprintf("Total time elapsed=%d seconds (i.e., %dmin; i.e., %4.2fhrs)\n\n\n",as.integer(tmpt),as.integer(tmpt/60),tmpt/3600)
cat(tmps)
cat(tmps,file=OUTFILE,sep="",append=TRUE)

} # LOOP OVER ps
} # LOOP OVER Findices

###################
#                 #
# Timing          #
#                 #
###################
# Note: for n=6e3, rqss runs for ~10min and then gets error "cannot allocate vector of size 9.9 Gb" in diag(X %*% summary(object, cov = TRUE)$V %*% t(X))
comp.time.fn <- function(reps=1,RQSS=TRUE) {
  require("quantreg");  source("quantile_inf_np.R")
  p <- 0.5;  BREPLIC <- 299;  NREPLIC <- reps;  ALPHA <- 0.05
  CV <- qnorm(1-ALPHA/2)
  cat(sprintf("NREPLIC=%d,BREPLIC=%d,p=%g",NREPLIC,BREPLIC,p),sep='\n')
  rxfn <- function(n) runif(n)
  ryfn <- function(n) runif(n)
  ns <- c(4e2,1e3,4e3,1e4,1e5,1e6)
  time.rqss.base <- rep.int(0,length(ns))
  rqss.fits <- vector("list",length(ns)*NREPLIC)
  dim(rqss.fits) <- c(NREPLIC,length(ns))
  NUM.X0s <- 10^(1:3)
  for (inum in 1:length(NUM.X0s)) {
    NUM.X0 <- NUM.X0s[inum]
    cat(sprintf("NUM.X0=%d",NUM.X0),sep='\n')
    time.gk <- time.bs <- rep.int(0,length(ns))
    time.rqss.pred <- rep.int(0,length(ns))
    for (nind in 1:length(ns)) {
      n <- ns[nind]
      CUBEADJ <- n^(1/12)
      set.seed(112358) #doesn't really matter, but good habit...
      for (irep in 1:NREPLIC) {
        X <- rxfn(n);  Y <- ryfn(n)
        x0s <- quantile(x=X,probs=seq(from=0.25,to=0.75,length.out=NUM.X0))
        # RQSS
        if (RQSS && n<=4e3) {
          if (inum==1) {
            tmpbase <- system.time({
              g <- function(lam,y,x) AIC(rqss(y ~ qss(x, lambda=lam), tau=p),k = -1)
              suppressWarnings({
                lamstar <- optimize(g, interval = c(0.001, 0.5), x = X, y = Y) 
                rqss.fits[[irep,nind]] <- rqss(Y~qss(X, lambda=lamstar$min), tau=p)
              })
            })[1]
            time.rqss.base[nind] <- time.rqss.base[nind] + tmpbase
          }
          tmppred <- system.time({
            suppressWarnings({
              tmp <- predict.rqss(rqss.fits[[irep,nind]],newdata=data.frame(X=x0s),interval="confidence",level=1-ALPHA)
            })
          })[1]
          time.rqss.pred[nind] <- time.rqss.pred[nind] + tmppred
        }
        
        # GK
        tmptot <- system.time({
          hs <- quantile.inf.np(Y=Y,X=X,x0s=x0s,p=p,ALPHA=ALPHA,JOINT.FLAG=FALSE,ONESIDED=0,LOCAL=TRUE)[[1]]$bandwidth.pointwise
        })[1]
        tmpgkonly <- system.time({
          quantile.inf.np(Y=Y,X=X,x0s=x0s,p=p,ALPHA=ALPHA,JOINT.FLAG=FALSE,ONESIDED=0,hs.pt=hs)
        })[1]
        time.gk[nind] <- time.gk[nind] + tmptot
        
        # local cubic
        tmpbs <- system.time({
          cubehs <- apply(cbind(x0s-min(X),max(X)-x0s,CUBEADJ*hs),1,min)
          CI.lo <- CI.hi <- rep.int(NA,NUM.X0)
          for (ix in 1:NUM.X0) {
            h <- cubehs[ix]
            locind <- (x0s[ix]-h <= X) & (X <= x0s[ix]+h)
            Xloc <- X[locind]-x0s[ix];  Yloc <- Y[locind]
            if (length(Yloc)<4) {
              warning("length(Yloc)<4")
              CI.lo[ix] <- CI.hi[ix] <- NA
              next
            }
            Xdes <- cbind(1,Xloc/h,Xloc^2/h^2,Xloc^3/h^3)
            useNA <- tryCatch({
              rq.ret <- rq(Yloc~Xdes[,2]+Xdes[,3]+Xdes[,4],tau=p)
              FALSE},
              error=function(rqERR) { return(TRUE) }
            )
            if (useNA) {
              warning("estimation error in rq")
              CI.lo[ix] <- CI.hi[ix] <- NA
              next
            }
            est <- rq.ret$coefficients[1]
            tryCatch({
              SE <- summary(rq.ret,se='boot',R=BREPLIC)$coefficients[1,2]
              CI.lo[ix] <- est - CV * SE
              CI.hi[ix] <- est + CV * SE 
            }, error=function(ERR){
              warning(ERR$message)
              CI.lo[ix] <- CI.hi[ix] <- NA
            })
          } # ix in x0s loop
        })[1]
        time.bs[nind] <- time.bs[nind] + (tmptot-tmpgkonly) + tmpbs
      }
    }
    cat(paste0("$L$-stat              ",sprintf("& %4d ",NUM.X0),
               paste0(sprintf("& %8.2f ",time.gk/NREPLIC),collapse=''),' \\\\',collapse=''),sep='\n')
    cat(paste0("Local cubic bootstrap ",sprintf("& %4d ",NUM.X0),
               paste0(sprintf("& %8.2f ",time.bs/NREPLIC),collapse=''),' \\\\',collapse=''),sep='\n')
    cat(paste0("\texttt{rqss}         ",sprintf("& %4d ",NUM.X0),
               paste0(sprintf("& %8.2f ",(time.rqss.base+time.rqss.pred)/NREPLIC),collapse=''),' \\\\',collapse=''),sep='\n')
  }
  cat(sprintf("n=%d",ns),sep=' '); cat("\n")
  cat(sprintf("rxfn=%s, ryfn=%s",as.character(as.expression(body(rxfn))),as.character(as.expression(body(ryfn)))))
}

# Feedback: KaplanDM@Missouri.edu
# Code for Kaplan (2014) "Nonparametric inference on quantile marginal effects"
# qme.inf() to calculate confidence intervals for quantile marginal effects (QME)
# plot.qme.inf() to plot results
# Other functions herein are for internal use only
# Examples: qme_inf_examples.R
# Simulations: qme_inf_sims.R


###################
#                 #
# Load packages   #
#                 #
###################
# sink("NUL")
# Load code from https://github.com/kaplandm/R
source("quantile_inf.R")
source("quantile_inf_np.R")
# sink()

###################
#                 #
# Main function   #
#                 #
###################
# Y is a vector of observations of a continuous scalar outcome variable
# X is a vector of observations of a continuous scalar regressor of interest
# To condition on discrete control variables, simply run this method on the subsample of observations with the discrete value(s) of interest.
# Conditioning on continuous control variables is not yet implemented.
# x0s is a vector of values of X at which to build CIs for the quantile marginal effect.
# p is the quantile of interest, a scalar between 0 and 1 (e.g., 0.5 for the median).
# ALPHA is between 0 and 1 such that desired pointwise coverage probability is (1-ALPHA); e.g., ALPHA=0.05 for 95% confidence intervals.
# ONESIDED=0 for two-sided confidence intervals, 1 for upper one-sided and -1 for lower one-sided.
# KERNEL.TYPE is for nonparametric conditional density estimation (for plug-in bandwidth); can be gaussian (default), epanechnikov, or uniform.
# BETA.BLOCK.SIZE and BETA.BLOCKS are passed directly to quantile.inf(), determining how many beta distribution draws are used in calibration.
# NORM.APPROX=TRUE uses the normal approximation for calibrating alpha-tilde.
# hs gives the option to use your own pre-computed bandwidths, which should be 1/4 the total local window size (x0-2h to x0 for lower window and x0 to x0+2h for upper)
# HMULT is a constant by which to multiply the baseline plug-in bandwidth.  It is not used when hs is non-null.
# GAMMA.EST should be TRUE to estimate the nuisance parameter PDF ratios, or FALSE to use the approximate value of 1.
# LOCAL: TRUE to use local polynomial bandwidth plug-in object estimation (faster), FALSE for B-spline-based; NA means local for larger n, spline for smaller n.
# Return values: x0s (same as argument); hs (bandwidth for each x0); NL and NR (local sample sizes to left and right of each x0); CI.lo and CI.hi (QME confidence interval endpoints); methname (GoldmanKaplan, Kaplan, bootstrap-perc-t, or NA) for which method at each x0
qme.inf <- 
  function(Y,X,x0s,p=NULL,ALPHA=0.05,ONESIDED=0,KERNEL.TYPE='gaussian',BETA.BLOCK.SIZE=10^4,BETA.BLOCKS=5, NORM.APPROX=FALSE,hs=NULL,HMULT=1,GAMMA.EST=FALSE, QTE.GK.ONLY=FALSE,QTE.NO.BS=TRUE, QTE.NO.GK=FALSE, LOCAL=TRUE, DENOM.MIN=NULL) {
    # Argument validation and value initialization
    validated <- quantile.inf.np.validate.inputs(Y=Y,X=X,x0s=x0s,p=p,ALPHA=ALPHA,JOINT.FLAG=FALSE,ONESIDED=ONESIDED,KERNEL.TYPE=KERNEL.TYPE,METHOD.TYPE='single',PSI=NULL,hs.pt=hs,hs.joint=NULL)
    if (!is.null(hs) && length(hs)!=length(x0s)) { stop("Length of hs must be same as of x0s.") }
    if (is.null(HMULT)) { HMULT <- 1 }
    if (!is.logical(LOCAL)) { stop("LOCAL must be TRUE, FALSE, or NA.") }
    for (varname in c('QTE.GK.ONLY','QTE.NO.BS','QTE.NO.GK')) {
      if (!is.logical(get(varname)) || is.na(get(varname))) stop(sprintf("%s must be TRUE or FALSE.",varname))
    }
    GAMMA <- NULL
    if (!is.logical(GAMMA.EST) || is.na(GAMMA.EST)) {
      stop("GAMMA.EST must be TRUE or FALSE.")
    } else if (!GAMMA.EST) {
      GAMMA <- list(t=1,c=1)
    }
    n <- length(Y)
    if (is.null(DENOM.MIN)) { DENOM.MIN <- (1e6/n)^(1/3) }
    d <- 1 #currently no continuous control variables allowed
    zz1 <- qnorm(1-ALPHA,0,1)
    if (p<1/n || p>=1-1/n) {
      stop("Not enough observations for inference; inference based on extreme value theory is recommended instead.")
    }
    
    #Scale-invariance (for sake of qsreg)
    X.sd <- sd(X) * sqrt(12)
    X <- X/X.sd #need to account for this later
    x0s <- x0s/X.sd
    if (!is.null(hs)) { hs <- hs/X.sd }
    
    # Compute bandwidth from quantile_inf_np.R
    if (is.null(hs)) {
      list[Y.hats,f.X.hats,fp.X.hats,Fp.hats,Fpp.hats] <- 
        quantile.inf.np.est.objs(Y=Y,X=X,x0s=x0s,p=p,KERNEL.TYPE=KERNEL.TYPE,LOCAL=LOCAL)
      hs.pt <- quantile.inf.np.bandwidths(p=p,n=n,X=X,x0s=x0s,f.X.hats=f.X.hats,fp.X.hats=fp.X.hats,Fp.hats=Fp.hats,Fpp.hats=Fpp.hats,ONESIDED=0,zz1=zz1,DENOM.MIN=DENOM.MIN)
      hs.pt <- hs.pt * n^(4/((18+7*d)*(2+d))) #different optimal rate, per Kaplan (2014)
      #shrink away from CPE-optimal toward (almost) MSE-optimal for larger n.
        kap <- (24/((18+7*d)*(6+d)))*pmax(0,1-4*n^(-1/5))*(0.95)
        hs.pt <- hs.pt * (1/3)*(2+n^kap)
      hs.pt <- quantile.inf.np.bandwidth.correct(hs=hs.pt,x0s=x0s,HCOR=HCOR)
      hs.pt <- hs.pt / 2  #since together windows span x0-2h to x0+2h, not just +/-h
      hs.pt <- hs.pt * HMULT
      maxhs <- apply(cbind(x0s-min(X),max(X)-x0s),1,min) #assuming scalar regressor
      hs.pt <- apply(cbind(maxhs/2,hs.pt),1,min)
    } else {
      hs.pt <- hs
    }
    
    # Gather local samples from [x0-2h,x0) and [x0,x0+2h] and run unconditional QTE inference
    NL <- NR <- CI.lo <- CI.hi <- methname  <- array(NA,length(x0s))
    for (ix in 1:length(x0s)) {
      if (is.na(hs.pt[ix])) { CI.lo[ix] <- CI.hi[ix] <- NA; next }
      indL <- which(((x0s[ix]-2*hs.pt[ix])<=X) & (X<x0s[ix]))
      indR <- which((x0s[ix]<=X) & (X<=(x0s[ix]+2*hs.pt[ix])))
      YL <- Y[indL];  YR <- Y[indR]
      NL[ix] <- length(YL);  NR[ix] <- length(YR)
      if (NL[ix]<3 || NR[ix]<3) { 
        CI.lo[ix] <- CI.hi[ix] <- methname[ix] <- NA 
      } else {
        # Use Kaplan instead of GK?
        #The following is alpha-tilde with normal approx and gamma=1, to get rough idea:
        tmptheta <- (1+sqrt(NR[ix]) / sqrt(NL[ix]))/sqrt(1+NR[ix]/NL[ix])
        tmpa <- 2*pnorm(qnorm(ALPHA/2)/tmptheta)
        tmp <- TRUE;  t.fn <- function(u,n) {(u<=1/(n+1))} #2?
        edge <- any(t.fn(quantile.inf.CIul(p=p,n=NL[ix],a=tmpa/2,APPROX=tmp),NL[ix]),
                    t.fn(1-quantile.inf.CIuh(p=p,n=NL[ix],a=tmpa/2,APPROX=tmp),NL[ix]),
                    t.fn(quantile.inf.CIul(p=p,n=NR[ix],a=tmpa/2,APPROX=tmp),NR[ix]),
                    t.fn(1-quantile.inf.CIuh(p=p,n=NR[ix],a=tmpa/2,APPROX=tmp),NR[ix]))
        if (QTE.GK.ONLY && edge) {
          useNA <- TRUE
        } else {
          useNA <- tryCatch( {
            tmp <- quantile.inf(X=list(t=YR,c=YL),p=p,ONESIDED=ONESIDED,ALPHA=ALPHA,METHOD.TYPE='qte',NORM.APPROX=NORM.APPROX,SPACING.FLAG=TRUE, BETA.BLOCK.SIZE=BETA.BLOCK.SIZE, BETA.BLOCKS=BETA.BLOCKS,GAMMA=GAMMA,NUM.INT=TRUE, QTE.GK.ONLY=QTE.GK.ONLY, QTE.NO.BS=QTE.NO.BS, QTE.NO.GK=(QTE.NO.GK||edge))
            CI.lo[ix] <- tmp$CI$lo / (2*hs.pt[ix])
            CI.hi[ix] <- tmp$CI$hi / (2*hs.pt[ix])
            methname[ix] <- tmp$methname
            FALSE } ,
            error=function(ME) {
              if (grepl("have p as close to zero",ME$message,fixed=TRUE) ||
                    grepl("1/n<=p<1-1/n",ME$message,fixed=TRUE)) {
                return(TRUE)
              } else {
                stop(sprintf("Error running quantile.inf() with NL[ix]=%d,NR[ix]=%d,p=%5g,hs.pt[ix]=%8g,ix=%d: %s", NL[ix],NR[ix],p,hs.pt[ix],ix,ME$message))
              }
            }
          )
        }
        if (useNA) { CI.lo[ix] <- CI.hi[ix] <- methname[ix] <- NA }
      }
    }
    
    # At each x0, return: x0, bandwidth, local sample sizes (2/ea), CI endpoints
    return(list(ALPHA=ALPHA,p=p,x0s=x0s*X.sd,hs=hs.pt*X.sd,NL=NL,NR=NR,CI.lo=CI.lo/X.sd,CI.hi=CI.hi/X.sd,methname=methname))
}

###################
#                 #
# Plot function   #
#                 #
###################
# NOTE: you may need to open device (X11, quartz, etc.) *before* calling this function.
plot.qme.inf <- function(qme.inf.obj,title=NULL,subtitle=NULL,x.axis.title='X',y.axis.title='Y',xlim=NULL,ylim=NULL,CIlwd=3) {
  par(family="serif")
  if (is.null(title)) { title=sprintf('IDEAL %g%% CI for %4.2f-QME',100*(1-qme.inf.obj$ALPHA),qme.inf.obj$p) }
  plot(c(qme.inf.obj$x0s,qme.inf.obj$x0s),c(qme.inf.obj$CI.lo,qme.inf.obj$CI.hi),type='n',
       main=title, sub=subtitle, xlab=x.axis.title, ylab=y.axis.title, cex.main=2.0, mgp=c(2.1,0.5,0), cex.axis=1.5, cex.lab=2, xlim=xlim, ylim=ylim)
  points(qme.inf.obj$x0s,qme.inf.obj$CI.lo,type='o',col=1,pch=1,lwd=CIlwd,lty=2)
  points(qme.inf.obj$x0s,qme.inf.obj$CI.hi,type='o',col=1,pch=1,lwd=CIlwd,lty=2)
}

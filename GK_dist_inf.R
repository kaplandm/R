# Feedback: David M. Kaplan  kaplandm.github.io
#
# Better version: distcomp Stata command
# https://kaplandm.github.io/#distcomp
# Install within Stata with command
#   net from https://kaplandm.github.io/stata
# and follow instructions.
#
# Warning: for use with continuous distributions/data; may not work correctly with discrete X and/or Y.  (Stata distcomp is better.)
#
# Stata data: try package  readstata13  (now supports [at least] version 15...)
# Or use "distcomp" Stata command
#
# Main (external) functions below: 
#   GK.dist.inf : uniform confidence bands, hypothesis testing, p-values
#   GK.dist.inf.plot.1s  and  GK.dist.inf.plot.2s : plotting uniform bands
# Other functions are just helper functions, though may also be called externally if helpful.
#
# Author (original): David M. Kaplan
# First version: July 8, 2013; last update Oct. 2024 (improved GK.dist.inf.rej.2s.r when duplicate sample values)
# Implementation of Dirichlet-based distributional inference
#   1-sample and 2-sample
#   1-sided  and 2-sided
#   global/GOF and multiple testing
# References (see https://kaplandm.github.io and links therein)
# "Comparing distributions by multiple testing across quantiles or CDF values"
#   Goldman & Kaplan, 2018 Journal of Econometrics, https://doi.org/10.1016/j.jeconom.2018.04.003
# "distcomp: Comparing distributions"
#   Kaplan, 2019 Stata Journal (& Stata command), https://doi.org/10.1177/1536867x19893626
# "Fractional order statistic approximation for nonparametric conditional quantile inference"
#   Goldman & Kaplan, 2017 Journal of Econometrics, https://doi.org/10.1016/j.jeconom.2016.09.015
# "Calibration for simultaneity: (re)sampling methods for simultaneous inference with applications to function estimation and functional data"
#   Buja & Rolke, 2006 WP, http://stat.wharton.upenn.edu/~buja/PAPERS/paper-sim.pdf


prob.loaded <- exists("quantile.inf")
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
# Load 2-sample lookup tables (1-sample uses formula)
#
#lookup (2-sample): ALPHA, nx, ny, alpha.tilde
GK.dist.inf.2s.lookup <- NULL
success <-
  tryCatch({GK.dist.inf.2s.lookup <- read.csv("https://raw.githubusercontent.com/kaplandm/R/main/GK_dist_inf_2s_lookup.txt");TRUE},
           error=function(w) FALSE)
if (!success) {
  success <-
    tryCatch({GK.dist.inf.2s.lookup <- read.csv("GK_dist_inf_2s_lookup.txt");TRUE},
             error=function(w) FALSE)
  if (!success) stop("Failed to load GK_dist_inf_2s_lookup.txt from web or local file. It should be available for download at https://github.com/kaplandm/R")
}

SIGFIGa <- 8

#
# Main (external) function: delegate to *.1s (if Y=NULL) or *.2s (otherwise),
# and to p-value function if alpha=NULL and CI/test otherwise.
#
# X: vector with observations from first (and perhaps only) sample
# Y: *either* vector with observations from second sample (see GK.dist.inf.2s function below) *or* NULL for 1-sample inference (see GK.dist.inf.1s function below)
# alpha: *either* NULL to compute p-value *or* compute uniform confidence band with confidence level (1-alpha), e.g. alpha=0.05 for 95% band, and corresponding hypothesis test with significance level alpha (e.g. 5%).  Note that when alpha=NULL, only the scalar p-value is returned, whereas otherwise a large list is (see GK.dist.inf.1s and GK.dist.inf.2s for details on the returned object).  Note: p-values do not use pretest or stepdown power improvements, so they may actually be smaller than reported (i.e., the returned p-value may be conservative). 
# PLOT.FLAG: if alpha is not NULL, then plot the uniform confidence bands if PLOT.FLAG=TRUE and do not plot if FALSE.
# tt, ttY: keep as NULL for now.
# FRAC: keep as FRAC=1 for now.
# H0.pval: for 1-sample p-value, null hypothesis CDF evaluated at the points in X (e.g., H0.pval=pnorm(x) to test against standard normal distribution); see GK.dist.inf.1s.pval below.  (For 2-sample p-value, null hypothesis is equality.)
# VERBOSE, MIN.DRAWS: arguments pertaining to just-in-time simulation for p-value function *or* 1-sample testing with PRETEST.FLAG or STEPDOWN.FLAG.  Computation of p-value should take well under 30 seconds; if not, decrease MIN.DRAWS for exploratory analysis.  For final analysis, a value higher than the default is recommended for increased accuracy, such as 2e5.
# PARALLEL: for 2-sample p-values, can optionally set to an integer greater than 1 corresponding to the number of CPUs available (e.g., Intel quad core processor has 4 CPUs).
# ONESIDED: 0 for two-sided, 
#          -1 for lower one-sided (in CDF space, i.e., H0:F()>=F0() or H0:Fx()>=Fy()), 
#          +1 for upper one-sided (in CDF space, i.e., H0:F()<=F0() or H0:Fx()<=Fy())
#          Note: the -1/+1 is the opposite in terms of quantile functions 
#                since, e.g., F()<=F0() is equivalent to Q()>=Q0()
# PRETEST.FLAG: for 1-sided hypothesis testing only, whether to pre-test for quantiles at which the null hypothesis is not binding to increase power.  Note: requires just-in-time simulation, which may take a while to compute; using PRETEST.FLAG=FALSE for exploratory/initial analysis and PRETEST.FLAG=TRUE for the final run is recommended.
# STEPDOWN.FLAG: whether to follow a Holm-type stepdown procedure to improve power while still controlling (strong) FWER.  If there are no rejections, it doesn't matter whether STEPDOWN.FLAG is TRUE or FALSE; if there are, just-in-time simulation is required for "stepping down," which may take a while to compute.  Using FALSE for exploratory/initial analysis and TRUE for the final run is recommended.
# H0.invCDF.fn: for one-sample hypothesis testing with either PRETEST.FLAG or STEPDOWN.FLAG, the null hypothesis inverse CDF function, i.e. the quantile function (e.g., qnorm).  Note: if you only have a CDF function, and it is continuous (which the method requires anyway), then you may construct a quantile function using uniroot(), like (example for standard normal): qfn <- function(q) uniroot(function(x)q-pnorm(x),lower=-100,upper=100)$root  (it need not be vectorized)
# 
# Two-sample rejections may be calculated by passing the returned object to GK.dist.inf.rej.2s.r as in
# ret <- GK.dist.inf(...);  H0r.rejections <- GK.dist.inf.rej.2s.r(ret)
GK.dist.inf <- function(X=NULL,Y=NULL,alpha=0.05,PLOT.FLAG=FALSE,tt=NULL,ttY=NULL,FRAC=1,
                        H0.pval=NULL,VERBOSE=FALSE,MIN.DRAWS=1e4,PARALLEL=1,
                        ONESIDED=0,PRETEST.FLAG=FALSE,STEPDOWN.FLAG=FALSE,H0.invCDF.fn=NULL,QTEONLY.FLAG=FALSE) {
  if (is.null(X)) { stop("X must contain vector of observed data.") }
  if (!is.logical(PLOT.FLAG) || !is.logical(PRETEST.FLAG) || !is.logical(STEPDOWN.FLAG)) stop("Arguments PLOT.FLAG, PRETEST.FLAG, and STEPDOWN.FLAG must all be either TRUE, FALSE, or NA.")
  if (PRETEST.FLAG || STEPDOWN.FLAG) {
    tryCatch(source('https://raw.githubusercontent.com/kaplandm/R/main/quantile_inf.R'),
             error=function(err)tryCatch(source("quantile_inf.R"),
                                         error=function(er2)stop("To use PRETEST.FLAG=TRUE or STEPDOWN.FLAG=TRUE, must have quantile_inf.R in the working directory (getwd()); couldn't load from GitHub either, but should be available at available at https://github.com/kaplandm/R")))
    if (PRETEST.FLAG && ONESIDED==0) {
      warning("PRETEST.FLAG=TRUE has no effect when ONESIDED==0.")
      PRETEST.FLAG <- FALSE
    }
  }
  if (!is.null(alpha) && ONESIDED!=0) { #same for one-sample or two-sample
    if (!is.null(Y)) { # && alpha %in% c(0.01,0.05,0.1,0.2)/2) {
      # Don't adjust, so can use lookup table (much faster)
    } else {
      alpha <- (2*alpha-alpha^2)/2 #asymptotically exact, but conservative in finite samples; e.g., consider one-sided 99% uniform confidence band.
    }
  }
  if (is.null(Y)) { #1-sample
    if (is.null(alpha)) {
      if (is.null(H0.pval)) { stop("H0.pval must contain vector of null hypothesis CDF values (evaluated at the observed X values).") }
      if (length(X)!=length(H0.pval)) { stop("X and H0.pval must be the same length.") }
      return(GK.dist.inf.1s.pval(X=sort(X),H0=sort(H0.pval),ONESIDED=ONESIDED)) #,VERBOSE=VERBOSE,MIN.DRAWS=MIN.DRAWS))
    } else {
      return(GK.dist.inf.1s(X=X,ALPHA=alpha,PLOT.FLAG=PLOT.FLAG,tt=tt,FRAC=FRAC,
                            ONESIDED=ONESIDED,PRETEST.FLAG=PRETEST.FLAG,STEPDOWN.FLAG=STEPDOWN.FLAG,
                            H0.invCDF.fn=H0.invCDF.fn,VERBOSE=VERBOSE,MIN.DRAWS=MIN.DRAWS,PARALLEL=PARALLEL))
    }
  } else{ #2-sample
    if (PRETEST.FLAG || STEPDOWN.FLAG) warning("PRETEST.FLAG and STEPDOWN.FLAG for two-sample inference have not been theoretically verified (although they appear to work well in simulations); use cautiously.")
    if (is.null(alpha)) {
      return(GK.dist.inf.2s.pval(X=sort(X),Y=sort(Y),VERBOSE=VERBOSE,MIN.DRAWS=MIN.DRAWS,PARALLEL=PARALLEL,ONESIDED=ONESIDED))
    } else {
      vals <- unique(GK.dist.inf.2s.lookup$ALPHA) #c(0.01,0.05,0.1)
      if (((alpha*(1+(ONESIDED!=0))) %in% vals) || (PRETEST.FLAG || STEPDOWN.FLAG)) {
        return(GK.dist.inf.2s(X=X,Y=Y,alpha=alpha,ONESIDED=ONESIDED,PLOT.FLAG=PLOT.FLAG,ttX=tt,ttY=ttY,FRAC=FRAC,STEPDOWN.FLAG=STEPDOWN.FLAG,PRETEST.FLAG=PRETEST.FLAG,QTEONLY.FLAG=QTEONLY.FLAG))
      } else {
        warning("Value of alpha not found in lookup table, so computing p-value for hypothesis test; may take longer to compute.")
        pval <- GK.dist.inf.2s.pval(X=sort(X),Y=sort(Y),VERBOSE=VERBOSE,MIN.DRAWS=MIN.DRAWS,PARALLEL=PARALLEL,ONESIDED=ONESIDED)
        return(list(pval=pval,rej=(pval<alpha)))
      }
    }
  }
}

#
# One-sample helper functions
#

# Function (main): 1-sample distributional inference.  Construct 
#  1-sample distributional CIs given data and desired confidence level.
#  X: vector of observed data.
#  alpha: desired confidence level is 1-alpha, e.g. alpha=0.05
#  for 95% confidence level or alpha=0.1 for 90% confidence level.
#  tt: optionally give different relative weights to different parts
#  of the distribution; by default, equal weight is given everywhere, so
#  the pointwise type I error is the same at all order statistics.  
#  FRAC: for now, FRAC=1 should not be changed; in the future, possibly FRAC<1
#  will be allowed, to interpolate between order statistics and construct
#  more than n pointwise CIs.  
#  PLOT.FLAG: if TRUE, graph confidence bands
#             (valid to interpret as band for CDF or quantile function).
#  Confidence bands: Connecting the returned CIs via a stair-step
#  pattern will give confidence bands with slightly conservative inference.
#  (The band should go straight up from each upper CI endpoint, 
#  then to the right horizontally to the next upper endpoint; and 
#  straight down from each lower endpoint, then horizontally to the left
#  to connect to the previous lower endpoint.)
#  Testing: reject any null hypothesis distribution that is not
#  contained entirely within the confidence bands.
GK.dist.inf.1s <- 
  function(X=NULL,ALPHA=0.05,PLOT.FLAG=FALSE,tt=NULL,FRAC=1,
           ONESIDED=0,PRETEST.FLAG=FALSE,STEPDOWN.FLAG=FALSE,
           H0.invCDF.fn=NULL,VERBOSE=FALSE,MIN.DRAWS=1e4,PARALLEL=1) {
    if (is.null(X)) { stop("X must contain vector of observed data.") }
    X <- sort(X)
    n <- length(X)
    if (is.null(tt)) { tt <- array(1/2,dim=c(n,2)) } else { stop("tt must be NULL (for now).") }
    if (FRAC!=1) { stop("FRAC must equal 1 (for now).") }
    if (is.null(H0.invCDF.fn) && (PRETEST.FLAG || STEPDOWN.FLAG)) stop("Must provide quantile function in argument H0.invCDF.fn in order to use PRETEST.FLAG=TRUE or STEPDOWN.FLAG=TRUE.")
    if (PRETEST.FLAG || STEPDOWN.FLAG) {
      oldseed <- NULL
      if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
        oldseed <- get(".Random.seed",.GlobalEnv)
      }
      on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
      if (PARALLEL>1) warning("PARALLEL not currently supported for one-sample power improvements.")
    }
    alpha.tilde <- GK.dist.1s.alpha.tilde(n,ifelse(ONESIDED!=0,2*ALPHA,ALPHA))
    CIs <- GK.dist.inf.1s.CIs(n=n,alpha.tilde=alpha.tilde,FRAC=FRAC,tt=tt)
    ret <- NULL #data.frame to be returned
    if (PRETEST.FLAG || STEPDOWN.FLAG) {
      qs <- c(CIs[,1],CIs[,2])
      suppressWarnings(H0.qs <- apply(matrix(qs,ncol=1),1,H0.invCDF.fn)) #warning() if Inf, when ONESIDED!=0
      tmp <- character(2*n);  tmp[1:n] <- ">=";  tmp[-(1:n)] <- "<="
      ret <- data.frame(OSind=c(1:n,1:n),Lstat=c(X,X),
                        reject=c(X<H0.qs[1:n],X>H0.qs[-(1:n)]),
                        quantile=c(CIs[,1],CIs[,2]),H0.dir=tmp,H0.val=H0.qs,
                        pre.rej=F,
                        ONEs=c(rep.int(-1,n),rep.int(1,n)))
      if (ONESIDED>0) ret <- ret[1:n,] else if (ONESIDED<0) ret <- ret[-(1:n),]
      rownames(ret) <- 1:dim(ret)[1]
      qs <- ret$quantile;  H0.qs <- NULL #use ret$H0.val instead
    } else {
      if (ONESIDED>0) CIs[,2] <- 1 else if (ONESIDED<0) CIs[,1] <- 0
      ret <- list(CIs=CIs,X=X)
      if (PLOT.FLAG) GK.dist.inf.plot.1s(ret)
      return(ret)
    }
    cat("Reminder: you may get confidence bands by setting PRETEST.FLAG=FALSE and STEPDOWN.FLAG=FALSE.\n",file="")
    N.BLOCKS <- NULL
    if (PRETEST.FLAG) { #set ret$pre.rej  (Note: already checked ONESIDED!=0)
      # Calibrate pre-test
      # setup
      N.BLOCKS <- 1;  BLOCK.SIZE <- MIN.DRAWS
      fail <- tryCatch(expr={Bs <- quantile.inf.betas(u=1:n/(n+1),n=n,reps=BLOCK.SIZE);F},
                       error=function(err){warning(err$message);TRUE})
      while (fail && BLOCK.SIZE>1) {
        N.BLOCKS <- N.BLOCKS*2
        BLOCK.SIZE <- ceiling(MIN.DRAWS/N.BLOCKS)
        fail <- tryCatch(expr={quantile.inf.betas(u=1:n/(n+1),n=n,reps=BLOCK.SIZE);F},
                         error=function(err)TRUE)
      }
      if (fail && BLOCK.SIZE<=1) {
        warning("Not enough memory for Dirichlet simulation; skipping pretest and returning.")
        return(ret[,-dim(ret)[2]])
      }
      pre.alpha <- ALPHA/log(log(max(15,n)))
      #iterate: guess alpha.tilde, check FWER of pre-test
      pre.FWER <- 0;  PRE.REL.TOL <- 0.1;  TILDE.FACTOR <- 1.5
      pre.OS <- pre.qs.ind <- NULL
      tmp.alpha.tilde <- GK.dist.1s.alpha.tilde(n,2*pre.alpha)/TILDE.FACTOR
      PREV.TOO.HIGH <- FALSE
      while (abs(pre.FWER-pre.alpha)/pre.alpha > PRE.REL.TOL) {
        if (pre.FWER>pre.alpha) {
          tmp.alpha.tilde <- tmp.alpha.tilde/TILDE.FACTOR 
          PREV.TOO.HIGH <- TRUE
        } else {
          if (PREV.TOO.HIGH) TILDE.FACTOR <- 1+(TILDE.FACTOR-1)/2
          PREV.TOO.HIGH <- FALSE
          tmp.alpha.tilde <- tmp.alpha.tilde*TILDE.FACTOR
          if (tmp.alpha.tilde>2*pre.alpha) break
        }
        tmp.CIs <- GK.dist.inf.1s.CIs(n=n,alpha.tilde=tmp.alpha.tilde)[,1+(ONESIDED>0)]
        suppressWarnings(tmp <- ifelse(ONESIDED>0,min(which(qs>=min(tmp.CIs))),max(which(qs<=max(tmp.CIs))))) #ok if +/-Inf, so suppress warning
        if (is.infinite(tmp)) {
          pre.FWER <- 0;  pre.qs.ind <- pre.OS <- NULL
        } else {
          if (ONESIDED>0) pre.qs.ind <- tmp:n else pre.qs.ind <- 1:tmp
          pre.qs <- qs[pre.qs.ind]
          pre.OS <- rep.int(NA,length(pre.qs))
          for (i in 1:length(pre.qs)) {
            if (ONESIDED>0) pre.OS[i] <- max(which(tmp.CIs<=pre.qs[i])) else pre.OS[i] <- min(which(tmp.CIs>=pre.qs[i]))
          }
          if (N.BLOCKS>1) set.seed(112358)
          blk.RPs <- rep.int(NA,N.BLOCKS)
          for (iB in 1:N.BLOCKS) {
            if (N.BLOCKS>1) Bs <- quantile.inf.betas(u=1:n/(n+1),n=n,reps=BLOCK.SIZE)
            n.rej <- 0
            for (i in 1:BLOCK.SIZE) n.rej <- n.rej + max(ONESIDED*Bs[i,pre.OS]>ONESIDED*pre.qs)
            blk.RPs[iB] <- n.rej/BLOCK.SIZE #percentage of rows/draws where *any* rejection
          }
          pre.FWER <- mean(blk.RPs)
        }
      }
      if (N.BLOCKS>1) rm(Bs)
      # Run pre-test using pre.OS and pre.qs.ind
      ret$pre.rej[pre.qs.ind] <- (ONESIDED*X[pre.OS]>ONESIDED*ret$H0.val[pre.qs.ind])
    }
    # Now, step down
    n.newrejs <- sum(ret$reject | ret$pre.rej)
    if (n.newrejs==0) {cat("No initial rejections, so no power improvement possible.\n");suppressWarnings(rm(Bs));return(ret[,-dim(ret)[2]])}
    if (!STEPDOWN.FLAG && sum(ret$pre.rej)==0) {cat("No pre-test rejections, so no power improvement possible.\n");suppressWarnings(rm(Bs));return(ret[,-dim(ret)[2]])}
    if (sum(ret$reject)==n) {cat("n rejections, so no need to step down.\n");suppressWarnings(rm(Bs));return(ret[,-dim(ret)[2]])}
    cat(sprintf("Re-calibrating to improve power. (This may take some time; you may disable all power improvements by setting STEPDOWN.FLAG=FALSE and PRETEST.FLAG=FALSE.)\n"),file="")
    set.seed(112358) #for replicability
    # determine blocksize/#blocks
    if (is.null(N.BLOCKS)) {
      N.BLOCKS <- 1;  BLOCK.SIZE <- MIN.DRAWS
      fail <- tryCatch(expr={Bs <- quantile.inf.betas(u=1:n/(n+1),n=n,reps=BLOCK.SIZE);F},
                       error=function(err){warning(err$message);TRUE})
      while (fail && BLOCK.SIZE>1) {
        N.BLOCKS <- N.BLOCKS*2
        BLOCK.SIZE <- ceiling(MIN.DRAWS/N.BLOCKS)
        fail <- tryCatch(expr={quantile.inf.betas(u=1:n/(n+1),n=n,reps=BLOCK.SIZE);F},
                         error=function(err)TRUE)
      }
      if (fail && BLOCK.SIZE<=1) {
        warning("Not enough memory for Dirichlet simulation; skipping stepdown and returning.")
        return(ret[,-dim(ret)[2]])
      }
    }
    SWITCH.N <- 20
    if (n>SWITCH.N) {
      skip <- tryCatch({RPmat <- matrix(0.1,nrow=n,ncol=sum(!ret$reject));FALSE},
                       error=function(err) { warning(sprintf("Skipping power improvement: %s.",err$message)); TRUE},
                       warning=function(w) { warning(sprintf("Skipping power improvement: %s.",w$message)); TRUE })
      if (skip) return(ret[,-dim(ret)[2]])
    }
    if (n>100 && N.BLOCKS>1) warning("This may take a very long time to compute; try freeing as much RAM as you can and re-running (if you haven't already).")
    rem.ind <- !(ret$reject | ret$pre.rej)
    if (!STEPDOWN.FLAG) rem.ind <- !ret$pre.rej
    rem.ONE <- ret$ONEs[rem.ind]
    rem.qs <- qs[rem.ind]
    RPmat <- NULL #for integer-OS method
    if (n>SWITCH.N) {
      RPmat <- matrix(NA,nrow=n,ncol=length(rem.qs))
      for (i in 1:dim(RPmat)[1]) {
        for (j in 1:dim(RPmat)[2]) {
          RPmat[i,j] <- pbeta(rem.qs[j],i,n+1-i,lower.tail=rem.ONE[j]<0)
        }
      }
    }
    n.rej.init <- sum(ret$reject)
    while (n.newrejs>0 && (!STEPDOWN.FLAG || sum(!(ret$reject | ret$pre.rej))>0)) {
      rem.ind <- !(ret$reject | ret$pre.rej)
      if (!STEPDOWN.FLAG) rem.ind <- !ret$pre.rej #only look at pre-test
      rem.ONE <- ret$ONEs[rem.ind]
      rem.qs <- qs[rem.ind]
      if (n<=SWITCH.N) { #use L-stat
        a.to.CIu.fn <- function(p,n,a,ONESIDED) {
          us <- rep.int(NA,length(p))
          for (i in 1:length(us)) {
            if (ONESIDED[i]<0) { 
              us[i] <- quantile.inf.CIuh(p=p[i],n=n,a=a,APPROX=FALSE)
            } else {
              us[i] <- quantile.inf.CIul(p=p[i],n=n,a=a,APPROX=FALSE)
            }
          }
          return(pmax(1/(n+1)+.Machine$double.eps,pmin(n/(n+1)-.Machine$double.eps,us)))
        }
        cal.fn <- function(a) { #similar to quantile_inf.R::quantile.inf.calibrate()
          us <- a.to.CIu.fn(rem.qs,n,a,rem.ONE)
          blk.RPs <- rep.int(NA,N.BLOCKS)
          set.seed(112358)
          for (iB in 1:N.BLOCKS) {
            Bs <- quantile.inf.betas(u=sort(us),n=n,reps=BLOCK.SIZE)
            Bs[,order(us)] <- Bs
            n.rej <- 0
            for (i in 1:BLOCK.SIZE) {
              n.rej <- n.rej + max(rem.ONE*Bs[i,]>rem.ONE*rem.qs)
            }
            blk.RPs[iB] <- n.rej/BLOCK.SIZE #percentage of rows/draws where *any* rejection
          }
          return(mean(blk.RPs)-ALPHA) #mean is linear fn, so take mean of means
        }
        ALPHAtilde <- tryCatch(uniroot(f=cal.fn, interval=c(.Machine$double.eps,ALPHA), f.lower=-ALPHA, tol=0.0001)$root,error=function(err) { if (err$message=="f() values at end points not of opposite sign") ALPHA else stop(err) })
        end.u <- a.to.CIu.fn(rem.qs,n,ALPHAtilde,rem.ONE)
        ret$OSind[rem.ind] <- end.u*(n+1)
        Lstats <- rep.int(NA,length(end.u))
        for (i in 1:length(end.u)) Lstats[i] <- quantile.inf.interp(X=X,u=end.u[i])
        ret$Lstat[rem.ind] <- Lstats
        ret$reject[rem.ind] <- ((rem.ONE*Lstats)>(rem.ONE*ret$H0.val[rem.ind]))
        n.newrejs <- sum(ret$reject[rem.ind])
      } else { #stick to order statistics
        rem.OSind <- ret$OSind[rem.ind]
        rejs <- rep.int(NA,BLOCK.SIZE*N.BLOCKS)
        if (N.BLOCKS>1) set.seed(112358)
        for (kB in 1:N.BLOCKS) {
          if (N.BLOCKS>1) Bs <- quantile.inf.betas(u=(1:n)/(n+1),n=n,reps=BLOCK.SIZE)
          for (iB in 1:BLOCK.SIZE) {
            rejs[(kB-1)*BLOCK.SIZE+iB] <- max((rem.ONE*Bs[iB,rem.OSind])>(rem.ONE*rem.qs))
          }
        }
        while (mean(rejs)<=ALPHA) { #keep increasing \tilde\alpha's as long as FWER<=ALPHA
          # find smallest aTilde incr from RPmat
          ats <- RPmat[cbind(pmax(1,pmin(dim(RPmat)[1],rem.OSind+rem.ONE)),
                             1:dim(RPmat)[2])] + 
            (rem.OSind==n & rem.ONE>0) + (rem.OSind==1 & rem.ONE<0) #make sure don't go out of bounds
          chg.ind <- which.min(ats)
          if (ats[chg.ind]>1) break
          # compute FWER, mean(rejs)
          tmp1 <- -rem.ONE[chg.ind];  tmpq <- rem.qs[chg.ind]
          newrejs <- rep.int(NA,BLOCK.SIZE*N.BLOCKS)
          if (N.BLOCKS==1) {
            newrejs <- ((tmp1*Bs[,rem.OSind[chg.ind]+rem.ONE[chg.ind]]) < (tmp1*tmpq)) & 
              ((tmp1*tmpq) <= (tmp1*Bs[,rem.OSind[chg.ind]]))
          } else {
            if (N.BLOCKS>1) set.seed(112358)
            for (kB in 1:N.BLOCKS) { #only draw the 2 relevant order statistics
              if (N.BLOCKS>1) Bs <- quantile.inf.betas(u=(1:n)/(n+1),n=n,reps=BLOCK.SIZE)
              Us <- c(rem.OSind[chg.ind]+rem.ONE[chg.ind],rem.OSind[chg.ind])/(n+1)
              Bs <- quantile.inf.betas(u=sort(Us),n=n,reps=BLOCK.SIZE)
              if (rem.ONE[chg.ind]>0) Bs <- Bs[,2:1] #to undo sort(Us)
              newrejs[((kB-1)*BLOCK.SIZE+1):(kB*BLOCK.SIZE)] <- 
                ((tmp1*Bs[,1]) < (tmp1*tmpq)) & ((tmp1*tmpq) <= (tmp1*Bs[,2]))
            }
          }
          rejs <- pmax(rejs,newrejs)
          if (mean(rejs)<=ALPHA) { #update ret
            rem.OSind[chg.ind] <- rem.OSind[chg.ind]+rem.ONE[chg.ind]
            ret$OSind[rem.ind] <- rem.OSind
          } else break
        }
        # update Lstat, reject, n.newrejs, and RPmat for next iteration
        ret$Lstat[rem.ind] <- X[ret$OSind[rem.ind]]
        ret$reject[rem.ind] <- ((rem.ONE*ret$Lstat[rem.ind])>(rem.ONE*ret$H0.val[rem.ind]))
        RPmat <- RPmat[,!ret$reject[rem.ind]]
        if (sum(!ret$reject[rem.ind])<=1) RPmat <- matrix(RPmat,ncol=1)
        n.newrejs <- sum(ret$reject[rem.ind])
      }
      if (!STEPDOWN.FLAG) {
        cat(sprintf("%d new rejection(s).\n",sum(ret$reject)-n.rej.init),file="")
        break #only go through once when PRETEST but not STEPDOWN
      }
      cat(sprintf("%d new rejection(s).\n",n.newrejs),file="")
    }
    rm(Bs); rm(RPmat)
    return(ret[,-dim(ret)[2]])
  }

# compute p-value:
#   X = observed data vector, sorted in increasing order.  (Not actually used, but helpful reminder.)
#   H0 = null hypothesis (CDF evaluated at X) vector, also sorted in increasing order.
#   VERBOSE: if TRUE, then print progress updates during simulation step.
#   MIN.DRAWS: set higher (e.g. 2e5) for increased accuracy.
#   Not currently supported: FRAC, TT.
GK.dist.inf.1s.pval <- function(X=NULL,H0=NULL,ONESIDED=NULL) { #,VERBOSE=FALSE,MIN.DRAWS=1e4) {
  n <- length(H0)
  k <- 1:n
  #Smallest alpha-tilde leading to any pointwise rejection
  a.min.L <- min(2*pbeta(H0[k],k,n+1-k))
  a.min.U <- min(2*(1-pbeta(H0[k],k,n+1-k)))
  if (ONESIDED==1) {
    if (a.min.U>1) return(1) #asymptotically, this is a probability zero event under F()=F0()
    a.min <- a.min.U
  } else if (ONESIDED==-1) {
    if (a.min.L>1) return(1) #asymptotically, this is a probability zero event under F()=F0()
    a.min <- a.min.L
  } else {
    a.min <- min(a.min.L,a.min.U)
  }
  # Compute corresponding alpha for that alpha.tilde
  TOL <- 0.0001
  p.lo <- 0.001;  p.hi <- 0.9
  a.lo <- aval <- GK.dist.1s.alpha.tilde(n=n,alpha=p.lo)
  if (a.min<aval) {
    aval <- GK.dist.1s.alpha.tilde(n=n,alpha=p.lo-TOL)
    if (a.min>aval) { return(p.lo) }
    p.hi <- p.lo-TOL;  a.hi <- aval
    p.lo <- a.lo <- 0
  }
  a.hi <- aval <- GK.dist.1s.alpha.tilde(n=n,alpha=p.hi)
  if (a.min>aval) {
    aval <- GK.dist.1s.alpha.tilde(n=n,alpha=p.hi+TOL)
    if (a.min<aval) { return(p.hi+TOL) }
    p.lo <- 0.9+TOL;  a.lo <- aval
    p.hi <- a.hi <- 1
  }
  while (p.hi-p.lo > TOL) {
    p.cur <- (p.lo+p.hi)/2
    a.cur <- GK.dist.1s.alpha.tilde(n=n,alpha=p.cur)
    if (a.cur>a.min) {
      p.hi <- p.cur;  a.hi <- a.cur
    } else {
      p.lo <- p.cur;  a.lo <- a.cur
    }
  }
  adj1s <- function(x) uniroot(f=function(p)2*p-p^2-x,lower=0,upper=1,f.lower=-x,f.upper=1-x)$root
  if (abs(a.min-a.hi)<abs(a.min-a.lo)) {
    if (ONESIDED==0) return(p.hi) else return(min(1,adj1s(p.hi)))
  } else {
    if (ONESIDED==0) return(p.lo) else return(min(1,adj1s(p.lo)))
  }
}

#
# Two-sample helper functions
#

# Function (main): 2-sample distributional inference.  The returned
#  "CIs" in CIs.x correspond to each order statistic of the X sample,
#  and similarly for CIs.y and Y.  rej=1 if a two-sided test rejects
#  equality of the two distributions, rej=0 otherwise.  alpha.sim is
#  the overall type I error found in simulations, which should be near
#  the true type I error (up to simulation error); this value should be
#  very close to the desired alpha (the argument), except if 1) the
#  combination of sample sizes is not found in the lookup table, in which
#  case the value should be slightly conservative (alpha.sim<alpha), or
#  2) the sample sizes are very small, in which case there is an issue
#  with discontinuity, and again a conservative value is used.  If (1)
#  is the case, you may use GK.dist.inf.2s.calibrate() to add an entry
#  to the lookup table for the sample sizes length(X) and length(Y).
#  X: vector of observed data from first sample.
#  Y: vector of observed data from second sample.  It is assumed that X and Y are
#  sampled independently of each other.
#  alpha: level of the test, e.g. alpha=0.05.
#  tt: optionally give different relative weights to different parts
#  of the distribution; by default, equal weight is given everywhere, so
#  the pointwise type I error is the same at all order statistics.  
#  FRAC: for now, FRAC=1 should not be changed; in the future, possibly FRAC<1
#  will be allowed, to interpolate between order statistics and construct
#  more than n pointwise CIs.  
#  ONESIDED: +1 for H0:Fx()<=Fy(), -1 for H0:Fx()>=Fy(), 0 for two-sided; note: opposite inequalities if writing quantile functions
#  PLOT.FLAG: if TRUE, then draw graph of rejection bands.
#  QTEONLY.FLAG: 
# test simulation code: NREPLIC=100;n=100;set.seed(112358);rfn=rnorm;rejs=matrix(NA,NREPLIC,7);system.time(for(irep in 1:NREPLIC) {tmp=GK.dist.inf.2s(X=rfn(n),Y=rfn(n),alpha=0.1,ONESIDED=1,STEPDOWN.FLAG=F,PLOT.FLAG=F,PRETEST.FLAG=F,QTEONLY.FLAG=TRUE);rejs[irep,]=tmp$prejs}); c(apply(rejs,2,mean),FWER=mean(apply(rejs,1,any)))
# test code output:
#   user  system elapsed 
# 301.75    5.27  924.83  
#                                    FWER 
# 0.03 0.01 0.01 0.03 0.04 0.03 0.01 0.10 
GK.dist.inf.2s <- function(X=NULL,Y=NULL,alpha=0.05,ONESIDED=0,
                           PLOT.FLAG=FALSE,
                           ttX=NULL,ttY=NULL,FRAC=1,
                           STEPDOWN.FLAG=FALSE,PRETEST.FLAG=FALSE,
                           QTEONLY.FLAG=NULL) {
  if (is.null(X)) { stop("X must contain vector of observed data.") }
  if (is.null(Y)) { stop("Y must contain vector of observed data.") }
  if (ONESIDED!=0 && PLOT.FLAG) stop("Can only use PLOT.FLAG if ONESIDED=0")
  n <- c(length(X),length(Y))
  if (is.null(ttX) && is.null(ttY)) { 
    ttX <- array(1/2,dim=c(n[1],2));  ttY <- array(1/2,dim=c(n[2],2))
  } else { stop("tt must be NULL (for now).") }
  if (FRAC!=1) { stop("FRAC must equal 1 (for now).") }
  if (PRETEST.FLAG && ONESIDED==0) stop("Cannot use PRETEST.FLAG=TRUE if ONESIDED=0")
  if (is.null(QTEONLY.FLAG)) QTEONLY.FLAG <- STEPDOWN.FLAG
  U <- list(x=sort(X),y=sort(Y))
  if (!PRETEST.FLAG && !QTEONLY.FLAG) {
    if (ONESIDED!=0) alpha <- 2*alpha
    tmp <- GK.dist.2s.alpha.tilde(min(n),max(n),alpha)
    alpha.tilde <- tmp[1];  alpha.sim <- tmp[2]
    CIs.x <- GK.dist.inf.1s.CIs(n=n[1],alpha.tilde=alpha.tilde,FRAC=FRAC,tt=ttX)
    CIs.y <- GK.dist.inf.1s.CIs(n=n[2],alpha.tilde=alpha.tilde,FRAC=FRAC,tt=ttY)
    if (ONESIDED>0) {
      CIs.x[,2] <- Inf;  CIs.y[,1] <- -Inf
    } else if (ONESIDED<0) {
      CIs.x[,1] <- -Inf;  CIs.y[,2] <- Inf
    }
    rej <- any(GK.dist.inf.rej.2s.r(list(CIs.x=CIs.x, CIs.y=CIs.y, U=U))$reject)
    # rej <- GK.dist.inf.rej.2s(n=n,L=list(xlo=CIs.x[,1],ylo=CIs.y[,1]),U=U,
    #                           H=list(xhi=CIs.x[,2],yhi=CIs.y[,2]))
    ret <- list(CIs.x=CIs.x,CIs.y=CIs.y,rej=rej,alpha.sim=alpha.sim,U=U)
    if (ONESIDED!=0) alpha <- alpha/2 #to undo above
  }
  if (PRETEST.FLAG || QTEONLY.FLAG || (rej==1 && STEPDOWN.FLAG)) {
    # Make sequence of evenly n^(-DELTA)-spaced quantiles
    DELTA <- 2/5
    tmp <- 1.15*round(min(n)^(DELTA))
    ps <- seq(from=1/(tmp+1),to=tmp/(tmp+1),length.out=tmp)
    prejs <- prerejs <- rep.int(FALSE,length(ps))
    # Figure out which quantile(s) already rejected
    ADDL <- 1
    if (PRETEST.FLAG) {
      # Need to flip ONESIDED for quantile (vs. CDF)...then again for PRETEST
      prerejs <- GK.joint.QTE.inf(X=U$x,Y=U$y,p=ps[!prejs],
                                  ALPHA=alpha*min(1,1/log(min(n))),ONESIDED=ONESIDED)
    } else if (STEPDOWN.FLAG && !QTEONLY.FLAG) {
      for (ip in 1:length(ps)) {
        p <- ps[ip]
        if (ONESIDED<=0) {
          OSy1 <- tryCatch(min(which(CIs.y[,1]>p)),warning=function(w)NA)
          OSx2 <- tryCatch(max(which(CIs.x[,2]<p)),warning=function(w)NA)
          if (U$y[OSy1]<U$x[OSx2]) prejs[ip] <- TRUE
        }
        if (ONESIDED>=0) {
          OSx1 <- tryCatch(min(which(CIs.x[,1]>p)),warning=function(w)NA)
          OSy2 <- tryCatch(max(which(CIs.y[,2]<p)),warning=function(w)NA)
          if (U$x[OSx1]<U$y[OSy2]) prejs[ip] <- TRUE
        }
      }
      ADDL <- sum(prejs)
    }
    # Loop GK.joint.QTE.inf() until no additional rejections
    while (ADDL>0 && sum(prejs)<length(ps[!prerejs])) {
      prejs.old <- prejs
      prejs[!prejs & !prerejs] <- 
        GK.joint.QTE.inf(X=U$x,Y=U$y,p=ps[!prejs & !prerejs],ALPHA=alpha,ONESIDED=-ONESIDED) #need to change ONESIDED for Quantile (vs. CDF)
      if (STEPDOWN.FLAG) ADDL <- sum(prejs)-sum(prejs.old) else ADDL <- 0
    }
    # Add ps and prejs (and prerejs) to ret
    if (PRETEST.FLAG) ret <- list(prerejs=prerejs) else if (QTEONLY.FLAG) ret <- list()
    ret$ps <- ps;  ret$prejs <- prejs
  }
  if (PLOT.FLAG && !PRETEST.FLAG && !QTEONLY.FLAG) {
    GK.dist.inf.plot.2s(ret)
  }
  return(ret)
}

#
# Helper function: joint QTE with all gamma=1
# Precondition: X and Y already sorted (ascending)
# p: vector of quantiles
# ONESIDED: +1 for H0:QTE(X-Y)<=0, -1 for >=0, 0 for two-sided; note: opposite of CDF inequality direction
GK.joint.QTE.inf <- function(X,Y,p,ALPHA,ONESIDED=0,BETA.BLOCK.SIZE=1e5,BETA.BLOCKS=5,ALPHAtilde=NULL) {
  # Set RNG seed for replicability (and save current seed to re-seed on exit)
  oldseed <- NULL
  if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
    oldseed <- get(".Random.seed",.GlobalEnv)
  }
  on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
  set.seed(112358) #for replicability
  # 
  n <- list(t=length(X),c=length(Y));  len.p <- length(p)
  if (len.p==1) {
    tmp <- quantile.inf(X=list(t=X,c=Y),p=p,ALPHA=ALPHA,ONESIDED=ONESIDED,METHOD.TYPE='qte',GAMMA=list(t=1,c=1))
    return(tmp$CI$lo>0 || tmp$CI$hi<0)
  }
  # similar to cal.fn(a) in quantile_inf.R
  cal.fn <- function(a) {
    ust <- usc <- rep.int(NA,(1+(ONESIDED==0))*len.p) #lower endpoint first
    for (i in 1:len.p) {
      if (ONESIDED>=0) ust[i] <- quantile.inf.CIul(p=p[i],n=n$t,a=a/(1+(ONESIDED==0))) else ust[i] <- quantile.inf.CIuh(p=p[i],n=n$t,a=a)
      if (ONESIDED<=0) usc[i] <- quantile.inf.CIul(p=p[i],n=n$c,a=a/(1+(ONESIDED==0))) else usc[i] <- quantile.inf.CIuh(p=p[i],n=n$c,a=a)
      if (ONESIDED==0) {
        ust[len.p+i] <- quantile.inf.CIuh(p=p[i],n=n$t,a=a/2)
        usc[len.p+i] <- quantile.inf.CIuh(p=p[i],n=n$c,a=a/2)
      }
    }
    blk.cps <- rep.int(NA,BETA.BLOCKS)
    for (i in 1:BETA.BLOCKS) {
      Bst <- quantile.inf.betas(u=sort(ust),n=n$t, reps=BETA.BLOCK.SIZE)
      Bst[,order(ust)] <- Bst
      Bsc <- quantile.inf.betas(u=sort(usc),n=n$c, reps=BETA.BLOCK.SIZE)
      Bsc[,order(usc)] <- Bsc
      #
      if (ONESIDED==0) {
        cover <- ((Bst[,1+len.p]>Bsc[,1]) & (Bsc[,1+len.p]>Bst[,1]))
        for (j in 2:len.p) {
          cover <- cover * ((Bst[,j+len.p]>Bsc[,j]) & (Bsc[,j+len.p]>Bst[,j]))
        }
      } else if (ONESIDED<0) {
        cover <- (Bst[,1]>Bsc[,1])
        for (j in 2:len.p) cover <- cover * (Bst[,j]>Bsc[,j])
      } else {
        cover <- (Bst[,1]<Bsc[,1])
        for (j in 2:len.p) cover <- cover * (Bst[,j]<Bsc[,j])
      }
      blk.cps[i] <- mean(cover)
    }
    return(mean(blk.cps) - (1-ALPHA))
  }
  if (is.null(ALPHAtilde)) {
    tmp <- cal.fn(0.99)
    if (tmp>0) {
      ALPHAtilde <- 0.99
    } else {
      ALPHAtilde <- uniroot(f=cal.fn, interval=c(0,0.99), tol=max(0.00001,min(0.001,max(n$t,n$c)^(-3/2))), f.upper=tmp, f.lower=ALPHA)$root
    }
  }
  #
  # Use ALPHAtilde to make CIs, rejs
  rejs <- rep.int(FALSE,len.p)
  for (ip in 1:len.p) {
    if (ONESIDED<=0) {
      uth <- quantile.inf.CIuh(p[ip],n$t,ALPHAtilde)
      ucl <- quantile.inf.CIul(p[ip],n$c,ALPHAtilde)
      tmp <- tryCatch({Qth <- quantile.inf.interp(X,uth)
      Qcl <- quantile.inf.interp(Y,ucl)
      (Qcl>Qth)},
      error=function(err) FALSE,
      warning=function(w) FALSE)
      if (tmp) rejs[ip] <- TRUE
    }
    if (ONESIDED>=0) {
      uch <- quantile.inf.CIuh(p[ip],n$c,ALPHAtilde)
      utl <- quantile.inf.CIul(p[ip],n$t,ALPHAtilde)
      tmp <- tryCatch({Qch <- quantile.inf.interp(Y,uch)
      Qtl <- quantile.inf.interp(X,utl)
      (Qtl>Qch)},
      error=function(err) FALSE,
      warning=function(w) FALSE)
      if (tmp) rejs[ip] <- TRUE
    }
  }
  return(rejs)
}

# Compute two-sided p-value:
#   X = observed data, sample #1
#   Y = observed data, sample #2
#   VERBOSE: if TRUE, then print progress updates during simulation step
#   MIN.DRAWS: set higher for increased accuracy
#   Not currently supported: FRAC, TT
GK.dist.inf.2s.pval <- function(X=NULL,Y=NULL,VERBOSE=FALSE,MIN.DRAWS=1e4,PARALLEL=1,ONESIDED=ONESIDED) {
  start.time <- Sys.time()
  if (is.null(X)) { stop("X must contain vector of observed data.") }
  if (is.null(Y)) { stop("Y must contain vector of observed data.") }
  nx <- length(X); ny <- length(Y)
  X <- sort(X);  Y <- sort(Y)
  #   L <- U <- H <- vector("list",2) #L(ow) and H(igh) CI endpoints, plus order statistics
  U <- list(x=X,y=Y)
  FRAC <- 1 
  TT.x <- array(1/2,dim=c(ceiling((nx+1)/FRAC-1),2)) 
  TT.y <- array(1/2,dim=c(ceiling((ny+1)/FRAC-1),2)) 
  # Search for smallest alpha.tilde s.t. reject
  TOL.RATIO <- 1.01 #stop when a.hi/a.lo<TOL.RATIO
  a.lo <- 0.100 #condition: always fail to reject w/ alpha.tilde=a.lo
  a.hi <- 0.900 #condition: always reject w/ alpha.tilde=a.hi
  # Check conditions first--adjust if don't hold (and return if really small)
  rej <- 1
  while (rej==1 && a.lo>1e3*.Machine$double.eps) {
    a.lo <- a.lo/10
    CIs.x <- GK.dist.inf.1s.CIs(n=nx,alpha.tilde=a.lo,FRAC=FRAC,tt=TT.x)
    CIs.y <- GK.dist.inf.1s.CIs(n=ny,alpha.tilde=a.lo,FRAC=FRAC,tt=TT.y)
    if (ONESIDED>0) {
      CIs.x[,2] <- Inf;  CIs.y[,1] <- -Inf
    } else if (ONESIDED<0) {
      CIs.x[,1] <- -Inf;  CIs.y[,2] <- Inf
    }
    rej <- any(GK.dist.inf.rej.2s.r(list(CIs.x=CIs.x, CIs.y=CIs.y, U=U))$reject)
    # rej <- GK.dist.inf.rej.2s(n=c(nx,ny),L=list(xlo=CIs.x[,1],ylo=CIs.y[,1]),U=U,
    #                           H=list(xhi=CIs.x[,2],yhi=CIs.y[,2]),
    #                           ONESIDED=ONESIDED)
  }
  if (a.lo<=1e3*.Machine$double.eps) { 
    return(GK.dist.inf.2s.verify(alpha.tildes=a.lo,nx=nx,ny=ny,VERBOSE=VERBOSE,MIN.DRAWS=MIN.DRAWS,PARALLEL=PARALLEL)$RP0s)
  }
  #
  rej <- 0
  while (rej==0 && a.hi<1-1e3*.Machine$double.eps) {
    a.hi <- 1-(1-a.hi)/10
    CIs.x <- GK.dist.inf.1s.CIs(n=nx,alpha.tilde=a.hi,FRAC=FRAC,tt=TT.x)
    CIs.y <- GK.dist.inf.1s.CIs(n=ny,alpha.tilde=a.hi,FRAC=FRAC,tt=TT.y)
    if (ONESIDED>0) {
      CIs.x[,2] <- Inf;  CIs.y[,1] <- -Inf
    } else if (ONESIDED<0) {
      CIs.x[,1] <- -Inf;  CIs.y[,2] <- Inf
    }
    rej <- any(GK.dist.inf.rej.2s.r(list(CIs.x=CIs.x, CIs.y=CIs.y, U=U))$reject)
    # rej <- GK.dist.inf.rej.2s(n=c(nx,ny),L=list(xlo=CIs.x[,1],ylo=CIs.y[,1]),U=U,
    #                           H=list(xhi=CIs.x[,2],yhi=CIs.y[,2]),
    #                           ONESIDED=ONESIDED)
  }
  if (a.hi>=1-1e3*.Machine$double.eps) { 
    return(GK.dist.inf.2s.verify(alpha.tildes=a.hi,nx=nx,ny=ny,VERBOSE=VERBOSE,MIN.DRAWS=MIN.DRAWS)$RP0s)
  }
  #
  while (a.hi/a.lo>TOL.RATIO) {
    a.cur <- (a.lo+a.hi)/2
    CIs.x <- GK.dist.inf.1s.CIs(n=nx,alpha.tilde=a.cur,FRAC=FRAC,tt=TT.x)
    CIs.y <- GK.dist.inf.1s.CIs(n=ny,alpha.tilde=a.cur,FRAC=FRAC,tt=TT.y)
    if (ONESIDED>0) {
      CIs.x[,2] <- Inf;  CIs.y[,1] <- -Inf
    } else if (ONESIDED<0) {
      CIs.x[,1] <- -Inf;  CIs.y[,2] <- Inf
    }
    rej <- any(GK.dist.inf.rej.2s.r(list(CIs.x=CIs.x, CIs.y=CIs.y, U=U))$reject)
    # rej <- GK.dist.inf.rej.2s(n=c(nx,ny),L=list(xlo=CIs.x[,1],ylo=CIs.y[,1]),U=U,
    #                           H=list(xhi=CIs.x[,2],yhi=CIs.y[,2]),
    #                           ONESIDED=ONESIDED)
    if (rej==0) { a.lo <- a.cur } else { a.hi <- a.cur }
  }
  # Simulate overall type I error for alpha.tilde found above
  tmp <- GK.dist.inf.2s.verify(alpha.tildes=a.hi,nx=nx,ny=ny,VERBOSE=VERBOSE,MIN.DRAWS=MIN.DRAWS)
  if (VERBOSE) {
    cat(sprintf("Overall time elapsed: %s\n",format(Sys.time()-start.time)))
  }
  adj1s <- function(x) uniroot(f=function(p)2*p-p^2-x,lower=0,upper=1,f.lower=-x,f.upper=1-x)$root
  pval <- tmp$RP0s;  if (ONESIDED!=0) pval <- adj1s(pval)
  return(pval)
}



#
# More one-sample helper functions
#

#Function: use prior simulated results to get calibrated alpha.tilde
# given sample size and desired alpha.
GK.dist.1s.alpha.tilde <- function(n,alpha) {
  if (n==1) { return(alpha) }
  ALPHAS <- c(0.001,0.01,0.05,0.1,0.2,0.5,0.7,0.9)
  nA <- length(ALPHAS)
  wgts1 <- 1-(alpha-ALPHAS[-nA])/(ALPHAS[-1]-ALPHAS[-nA])
  wgts1 <- c(ifelse(wgts1<0 | wgts1>1,0,wgts1),0)
  wgts2 <- 1-(ALPHAS[-1]-alpha)/(ALPHAS[-1]-ALPHAS[-nA])
  wgts2 <- c(0,ifelse(wgts2<0 | wgts2>1,0,wgts2))
  wgts <- pmax(wgts1,wgts2)
  if (n<4) {
    if (sum(wgts)==0) { stop(sprintf("Error in GK.dist.1s.alpha.tilde(): extrapolation of alpha-tilde for alpha<%g or alpha>%g is inaccurate for n=2 and n=3.",min(ALPHAS),max(ALPHAS))) }
    at.n2n3 <- rbind(c(0.00050841,0.00512235,0.02642680,0.05458812,0.11359550,0.31450684,0.47263448,0.68690046),
                     c(0.00034078,0.00353723,0.01852887,0.03882469,0.08242757,0.24019615,0.37413242,0.56906730))
    at.wgt <- rbind((n==2)*wgts,(n==3)*wgts)
    return(sum(at.wgt*at.n2n3))
  }
  c.combs <- 
    rbind(c( 4.436013,  2.140632,  0.3830217, -0.4311821, -1.023106, -2.169891, -2.3528158, -2.5381499),
          c( 4.864360,  4.701783,  4.6436968,  4.6331936,  4.422736,  4.294522,  3.9715237,  3.5937777),
          c( 3.144887,  2.923568,  3.2227556,  3.1971612,  3.088053,  2.859663,  2.679227,  2.505629),
          c(-3.462333, -3.095870, -2.9618197, -2.7534682, -2.792809, -2.321693, -2.3602093, -2.2852256))
  coeffs <- rep.int(NA,4)
  if (sum(wgts)==0) {
    warning("Using extrapolation to compute alpha-tilde; may be inaccurate for values of alpha extremely close to zero or one.")
    coeffs[1] <- -2.749199  -1.044047*log(alpha)
    coeffs[2] <- 4.758833   -1.196696*alpha
    coeffs[3] <- exp(1.1526013  -0.2386897*alpha)
    coeffs[4] <- -3.9560887 +1.7150781*alpha^0.1709839
    if (alpha<min(ALPHAS)) {
      ret <- min(exp(-coeffs[1]-coeffs[2]*sqrt(log(log(n)))-coeffs[3]*(log(n))^coeffs[4]),
                 exp(-c.combs[1,1]-c.combs[2,1]*sqrt(log(log(n)))-c.combs[3,1]*(log(n))^c.combs[4,1]))
    } else { #alpha>max(ALPHAS)
      ret <- max(exp(-coeffs[1]-coeffs[2]*sqrt(log(log(n)))-coeffs[3]*(log(n))^coeffs[4]),
                 exp(-c.combs[1,nA]-c.combs[2,nA]*sqrt(log(log(n)))-c.combs[3,nA]*(log(n))^c.combs[4,nA]))
    }
  } else {
    coeffs <- rowSums(c.combs*rbind(wgts,wgts,wgts,wgts))
  }
  return(exp(-coeffs[1]-coeffs[2]*sqrt(log(log(n)))-coeffs[3]*(log(n))^coeffs[4]))
}

# Function: construct CIs given n and \tilde\alpha
GK.dist.inf.1s.CIs <- function(n,alpha.tilde,FRAC=1,tt=array(1/2,dim=c(ceiling((n+1)/FRAC-1),2))) {
  ret <- array(0,dim=c(ceiling((n+1)/FRAC-1),2))
  for (k in seq(from=FRAC,by=FRAC,to=n+1-min(FRAC,1e-6))) { #1:n) {
    ret[k,1] <- qbeta(alpha.tilde*tt[k,1],k,n+1-k)
    ret[k,2] <- qbeta(1-tt[k,2]*alpha.tilde,k,n+1-k)
  }
  return(ret)
}

# Function: given CIs, reject null hypothesis? (1-sample)
GK.dist.inf.rej.1s <- function(X, CIs, H0.fn) {
  U <- sort(X)
  return(max((CIs[,1]>H0.fn(U)) | (H0.fn(U)>CIs[,2])))
}

# Simulate (joint) type I error given some alpha-tilde(s), 1-sample
# PARALLEL: # CPUs available (requires foreach and doParallel packages)
GK.dist.inf.1s.verify <- function(alpha.tildes=NULL,n=NULL,VERBOSE=FALSE,MIN.DRAWS=2e5,PARALLEL=1) {
  if (is.null(alpha.tildes) || is.null(n)) {stop("Error: argument NULL for alpha.tildes and/or n not allowed in function GK.dist.inf.1s.verify()")}
  if (PARALLEL>1) {
    if (!require("parallel") || !require("foreach")) {warning("Install package foreach in order to run PARALLEL."); PARALLEL <- 0}
    if (!require("doParallel")) {warning("Install package doParallel in order to run PARALLEL."); PARALLEL <- 0}
    PARALLEL <- tryCatch({workers <- makeCluster(PARALLEL); registerDoParallel(workers); on.exit(stopCluster(workers),add=TRUE); PARALLEL}, error=function(Err){warning(sprintf("Error creating %d clusters",PARALLEL)); 0})
  } else if (PARALLEL<0) {
    warning("PARALLEL must be a non-negative integer: 0 (or 1) for not parallel, positive for number of parallel CPUs to use.")
  }
  overall.start.time <- Sys.time()
  oldseed <- NULL
  if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
    oldseed <- get(".Random.seed",.GlobalEnv)
  }
  on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
  set.seed(112358) #for replicability
  # Constant parameters
  BLOCKSIZE <- min(max(1,round(1e6/n)),MIN.DRAWS)
  N.DRAWS <- BLOCKSIZE*ceiling(MIN.DRAWS/BLOCKSIZE) #2*10^5 #X*10^5; see paper
  FRAC <- 1 
  TT <- array(1/2,dim=c(ceiling((n+1)/FRAC-1),2)) 
  # Storage variables
  rejs <- matrix(NA,N.DRAWS,length(alpha.tildes))
  CIs <- vector("list",length(alpha.tildes))
  for (i in 1:length(alpha.tildes)) {
    CIs[[i]] <- t(GK.dist.inf.1s.CIs(n=n,alpha.tilde=alpha.tildes[i],FRAC=FRAC,tt=TT))
  }
  if (PARALLEL>1) {
    clusterSetRNGStream(workers,112358)
    rejs <- foreach(irep=1:(N.DRAWS/BLOCKSIZE),.combine=rbind,.inorder=FALSE) %dopar% {
      #       if (VERBOSE) {
      #         cat(sprintf("Start block %d/%d; ",irep,N.DRAWS/BLOCKSIZE),file="",sep='')
      #       }
      unirands <- matrix(runif(BLOCKSIZE*ceiling((n+1)/FRAC-1)),BLOCKSIZE)
      frac.order.stats <- NULL
      if (FRAC==1) {
        frac.order.stats <- t(apply(unirands,1,sort))
      } else {
        # Draw from Dirichlet for (fractional) uniform order statistics
        frac.order.stats <- array(0,dim=c(BLOCKSIZE,ceiling((n+1)/FRAC-1)))
        frac.order.stats[,1] <- qbeta(unirands[,1],FRAC,n+1-FRAC) #rbeta(N.DRAWS,FRAC,n+1-FRAC)
        for (i in 2:ceiling((n+1)/FRAC-1)) {
          frac.order.stats[,i] <- frac.order.stats[,i-1] + 
            qbeta(unirands[,i],FRAC,n+1-i*FRAC)*(1-frac.order.stats[,i-1])
        }
      }
      rejs <- NULL
      for (j in 1:length(alpha.tildes)) {
        #         rejs[(BLOCKSIZE*(irep-1)+1):(BLOCKSIZE*irep),j] <- 1 - 
        rejs <- cbind(rejs, 1 -
                        apply(frac.order.stats,1,function(x) prod(CIs[[j]][1,]<x,na.rm=TRUE))*
                        apply(frac.order.stats,1,function(x) prod(x<CIs[[j]][2,],na.rm=TRUE))
        )
      }
      rejs
    }
  } else {
    for (irep in 1:(N.DRAWS/BLOCKSIZE)) {
      if (VERBOSE) {
        cat(sprintf("Start block %d/%d; ",irep,N.DRAWS/BLOCKSIZE),file="",sep='')
      }
      unirands <- matrix(runif(BLOCKSIZE*ceiling((n+1)/FRAC-1)),BLOCKSIZE)
      frac.order.stats <- NULL
      if (FRAC==1) {
        frac.order.stats <- t(apply(unirands,1,sort))
      } else {
        # Draw from Dirichlet for (fractional) uniform order statistics
        frac.order.stats <- array(0,dim=c(BLOCKSIZE,ceiling((n+1)/FRAC-1)))
        frac.order.stats[,1] <- qbeta(unirands[,1],FRAC,n+1-FRAC) #rbeta(N.DRAWS,FRAC,n+1-FRAC)
        for (i in 2:ceiling((n+1)/FRAC-1)) {
          frac.order.stats[,i] <- frac.order.stats[,i-1] + 
            qbeta(unirands[,i],FRAC,n+1-i*FRAC)*(1-frac.order.stats[,i-1])
        }
      }
      for (j in 1:length(alpha.tildes)) {
        rejs[(BLOCKSIZE*(irep-1)+1):(BLOCKSIZE*irep),j] <- 1 - 
          apply(frac.order.stats,1,function(x) prod(CIs[[j]][1,]<x,na.rm=TRUE))*
          apply(frac.order.stats,1,function(x) prod(x<CIs[[j]][2,],na.rm=TRUE))
      }
    }
  }
  RP0s <- colMeans(rejs)
  return(list(RP0s=RP0s,n=n,N.DRAWS=N.DRAWS,alpha.tildes=alpha.tildes,
              time.elapsed=Sys.time()-overall.start.time,FRAC=FRAC,TT=TT))
}

# Plotting function, 1-sample
# NOTE: you may need to open device (X11, quartz, etc.) before calling this function.
GK.dist.inf.plot.1s <- 
  function(dist.inf.obj=NULL,main="exact confidence band",sub="",
           xlab="X",ylab="F(X)",xlim=NULL,ylim=NULL,xaxs='r',yaxs='r',
           lty="2121",pch=NULL,draw.legend=FALSE,legend="confidence band",
           panel.first=NULL, lwd=2, mgp=c(2.1,0.5,0), xaxp=NULL, yaxp=NULL,
           cex.main=2,cex.axis=2,cex.lab=2,parmar=c(5.0,6.0,6.0,2.1)) {
    if (is.null(dist.inf.obj)) {stop("Argument should be returned list from GK.dist.inf.1s")}
    X <- dist.inf.obj$X;  CIs <- dist.inf.obj$CIs
    U <- sort(dist.inf.obj$X);  n <- length(U)
    plot.band.x.L <- rbind(c(U[1]-100*max(abs(U)),0),
                           cbind(c(rep(U,each=2),U[n]),
                                 c(0,rep(CIs[,1],each=2))),
                           c(U[n]+100*max(abs(U)),CIs[n,1]))
    plot.band.x.H <- rbind(c(U[1]-100*max(abs(U)),CIs[1,2]),
                           cbind(c(U[1],rep(U,each=2)),
                                 c(rep(CIs[,2],each=2),1)),
                           c(U[n]+100*max(abs(U)),1))
    #   x11() #  quartz()
    par(family="serif",mar=parmar)
    plot(x=plot.band.x.L[,1],y=plot.band.x.L[,2],type="n",pch=1,col=1, 
         main=main, sub=sub, xlab=xlab, ylab=ylab, cex.main=cex.main, 
         mgp=mgp, cex.axis=cex.axis, cex.lab=cex.lab, 
         xlim=ifelse(rep.int(is.null(xlim),2),c(U[1],U[n]),xlim), 
         ylim=ifelse(rep.int(is.null(ylim),2),c(0,1),ylim), 
         xaxs=xaxs, yaxs=yaxs, panel.first=panel.first, lwd=lwd,
         xaxp=xaxp, yaxp=yaxp)
    lines(plot.band.x.L, type='l', col=1, pch=pch, lwd=lwd, lty=lty)
    lines(plot.band.x.H, type='l', col=1, pch=pch, lwd=lwd, lty=lty)
    if (draw.legend) {
      legend('bottomright', legend=legend, inset=0.01, col=1, pch=pch, lty=lty, lwd=lwd, cex=1.5, y.intersp=0.9, bg='white')
    }
  }

#
# Draw band on top of existing plot, for visualizing posterior
#
GK.dist.inf.band.lines <- function(dist.inf.obj=NULL,
                                   lwd=2,lty="2121",pch=NULL,col=1) {
  if (is.null(dist.inf.obj)) {stop("Argument should be returned list from GK.dist.inf.1s")}
  X <- dist.inf.obj$X;  CIs <- dist.inf.obj$CIs
  U <- sort(dist.inf.obj$X);  n <- length(U)
  plot.band.x.L <- rbind(c(U[1]-100*max(abs(U)),0),
                         cbind(c(rep(U,each=2),U[n]),
                               c(0,rep(CIs[,1],each=2))),
                         c(U[n]+100*max(abs(U)),CIs[n,1]))
  plot.band.x.H <- rbind(c(U[1]-100*max(abs(U)),CIs[1,2]),
                         cbind(c(U[1],rep(U,each=2)),
                               c(rep(CIs[,2],each=2),1)),
                         c(U[n]+100*max(abs(U)),1))
  lines(plot.band.x.L, type='l', col=col, pch=pch, lwd=lwd, lty=lty)
  lines(plot.band.x.H, type='l', col=col, pch=pch, lwd=lwd, lty=lty)
}

#
# Compute alpha-tilde entry and append to lookup table, 1-sample
#
# NS: matrix where top row is sample sizes for first sample, 2nd (and bottom) row is second sample sizes.  Loops through columns.
# ALPHAS: vector of alpha values (nominal type I error).
# N.DRAWS, CALIB.DEC.ERR, CALIB.ERR.PERC: control accuracy and precision of calibrated alpha-tilde.  It is recommended to not change these without a good understanding of Appendix B from Kaplan and Goldman (2013).  CALIB.DEC.ERR is the overall calibration error as decimal percent, including search tolerance *and* simulation error, while CALIB.ERR.PERC is probability of exceeding CALIB.DEC.ERR.
# OUTFILE: where to output the newly calibrated value(s).  If the file exists, then values will be appended to the end, otherwise a new file will be created.
# LOGFILE: where to output details of computations.
# VERBOSE: will output additional information if TRUE.
GK.dist.inf.1s.calibrate <- function(NS=NULL,ALPHAS=NULL,N.DRAWS=NULL,CALIB.DEC.ERR=NULL,CALIB.ERR.PERC=0.05,OUTFILE="GK_dist_inf_1s_lookup.txt",LOGFILE="GK_dist_inf_1s_lookup_log.txt",VERBOSE=TRUE) {
  FRAC <- 1 #only tested for FRAC=1, but in principle could use FRAC=1/k for some integer k for smaller n, or FRAC=k for larger n
  SIGFMTa <- sprintf("%%%d.%df",SIGFIGa+2,SIGFIGa)
  SIGFMTf <- SIGFMTa
  if (is.null(NS) || is.null(ALPHAS)) {stop("Error in GK.dist.inf.1s.calibrate: NULL value for NS and/or ALPHAS not allowed.")}
  if (is.null(CALIB.DEC.ERR)) {
    if (min(ALPHAS)>0.01) {CALIB.DEC.ERR <- 0.02} else {CALIB.DEC.ERR <- 0.05}
  }
  if (is.null(N.DRAWS)) {
    if (min(ALPHAS)>0.01) {N.DRAWS <- 2e5} else {N.DRAWS <- 1e6}
  }
  oldseed <- NULL
  if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
    oldseed <- get(".Random.seed",.GlobalEnv)
  }
  set.seed(112358) #for replicability  
  
  #Appends to end of specified file now, not overwrite (unless file doesn't exist)
  if (OUTFILE!="") {
    if (!file.exists(OUTFILE)) {  # Write to file: header w/ variable names
      cat(sprintf("Creating new file: %s",OUTFILE), file="",sep="\n")
      cat('"ALPHA","n","alpha.tilde","ALPHA.simulated"',file=OUTFILE,sep="\n",append=FALSE)
    } else { # append to existing file
      cat(sprintf("Appending to existing file: %s",OUTFILE), file="",sep="\n")
    }
  } else {
    cat("No file specified, so output is to console only.", file="",sep="\n")
  }
  
  # Loop over sample sizes
  NS.SORT <- sort(NS,decreasing=FALSE)
  if (FRAC==1) {
    NS.SORT <- sort(NS,decreasing=TRUE)
  }
  prev.alpha.tilde <- ALPHAS  #bad guess
  if (exists("GK.dist.inf.1s.lookup") && !is.null(GK.dist.inf.1s.lookup)) {
    for (k in 1:length(ALPHAS)) {
      tryCatch({
        lka <- GK.dist.inf.1s.lookup[GK.dist.inf.1s.lookup$ALPHA==ALPHAS[k],]
        tmpn <- sort(unique(lka$n))
        tmpind <- which.min(abs(NS.SORT[1]-tmpn))
        prev.alpha.tilde[k] <- lka$alpha.tilde[lka$n==tmpn[tmpind]]
      }, error=function(ERR) {})
    }
  }
  frac.order.stats <- NULL
  for (n in NS.SORT) {
    start.time <- Sys.time()
    TT <- array(1/2,dim=c(ceiling((n+1)/FRAC-1),2)) #equal type I error @ all quantiles, all equal-tailed
    if (FRAC==1) {
      if (n==NS.SORT[1]) {
        gc(verbose=VERBOSE)
        if (VERBOSE) cat(sprintf("Trying prefill..."))
        prefill <- tryCatch({
          frac.order.stats <- matrix(runif(N.DRAWS*n),ncol=N.DRAWS);TRUE},
          error=function(MemErr) {FALSE})
        if (!prefill) {
          tryCatch(frac.order.stats <- matrix(as.numeric(NA),nrow=n,ncol=N.DRAWS), error=function(MemErr) {stop(sprintf("\nNot enough memory for n=%d, N.DRAWS=%d; try reducing N.DRAWS (or suggesting that FRAC>1 be fully implemented)",n,N.DRAWS))} )
          #...so throw error in above line if just not enough room in memory.
          if (VERBOSE) {
            cat(sprintf("Due to memory constraints, will proceed in blocks.\n"),file="",sep='',append=TRUE)
            cat(sprintf("Due to memory constraints, will proceed in blocks.\n"),file=LOGFILE,sep='',append=TRUE)
          }
        }
        gc(verbose=VERBOSE)
        win <- FALSE
        for (NBLOCKS in c(1,2,4,8,16,32,64)) {
          if (!prefill && NBLOCKS<4) { next }
          if (VERBOSE) cat(sprintf("Trying NBLOCKS=%d;",NBLOCKS))
          success <- tryCatch({
            BLOCKSIZE <- floor(N.DRAWS/NBLOCKS)
            for (iblock in 1:NBLOCKS) {
              if (prefill) {
                frac.order.stats[,(((iblock-1)*BLOCKSIZE)+1):(iblock*BLOCKSIZE)] <- 
                  apply(frac.order.stats[,(((iblock-1)*BLOCKSIZE)+1):(iblock*BLOCKSIZE)],2,sort)
              } else {
                frac.order.stats[,(((iblock-1)*BLOCKSIZE)+1):(iblock*BLOCKSIZE)] <- 
                  apply(matrix(runif(n*BLOCKSIZE),nrow=n),2,sort)
              }
              if (VERBOSE) {
                cat(sprintf("Block %d/%d",iblock,NBLOCKS),file="",sep='',append=TRUE)
                cat(sprintf("Block %d/%d",iblock,NBLOCKS),file=LOGFILE,sep='',append=TRUE)
              }
            }
            if (N.DRAWS>BLOCKSIZE*NBLOCKS) {
              for (icol in (1+BLOCKSIZE*NBLOCKS):N.DRAWS) {
                frac.order.stats[,icol] <- sort(frac.order.stats)
              }
            }
            win <- TRUE
            TRUE}
            ,error=function(Err){cat(sprintf("Not enough memory for sorting in %d blocks...",NBLOCKS),file="",sep='',append=TRUE);FALSE})
          if (success) { break }
        }
        if (!win) { stop(sprintf("\nNot enough memory for n=%d, N.DRAWS=%d; try reducing N.DRAWS (or suggesting that FRAC>1 be fully implemented)",n,N.DRAWS)) }
        gc(verbose=VERBOSE)
      } else {
        if (FALSE) { #(prev.n-n==1) {
          gc(verbose=VERBOSE)
          tin <- sample.int(n=prev.n,size=N.DRAWS,replace=TRUE)
          for (icol in 1:N.DRAWS) {
            if (tin[icol]==1) {
              frac.order.stats[,icol] <- frac.order.stats[c(2:(n+1),-(2:(n+1))),icol]
            } else if (tin[icol]<=n) {
              frac.order.stats[,icol] <- frac.order.stats[c(1:(tin[icol]-1),(tin[icol]+1):(n+1),tin[icol],-(1:(n+1))),icol]
            } #else: dropping (n+1)th entry anyway, no need to reorder.
          }
          rm(tin)
          gc(verbose=VERBOSE)
        } else {
          ndiff <- prev.n-n
          for (icol in 1:N.DRAWS) {
            tin <- sample.int(n=prev.n,size=ndiff,replace=FALSE)
            keepind <- (1:prev.n)[-tin]
            frac.order.stats[,icol] <- 
              c(frac.order.stats[keepind,icol],frac.order.stats[-keepind,icol])
          }
        }
      }
    } else {
      if (n==NS.SORT[1]) {
        #unirands: more efficient/consistent to reuse random draws from prev NS
        unirands <- matrix(runif(N.DRAWS*ceiling((NS.SORT[1]+1)/FRAC-1)),ncol=N.DRAWS)
      } else {
        unirands <- rbind(unirands,matrix(runif(N.DRAWS*(ceiling((n+1)/FRAC-1)-dim(unirands)[1])),ncol=N.DRAWS))
      }
      # Draw from Dirichlet for (fractional) uniform order statistics
      frac.order.stats <- array(NA,dim=c(ceiling((n+1)/FRAC-1)),N.DRAWS)
      frac.order.stats[1,] <- qbeta(unirands[1,],FRAC,n+1-FRAC) #rbeta(N.DRAWS,FRAC,n+1-FRAC)
      if (dim(frac.order.stats)[1]>1) {
        for (i in 2:ceiling((n+1)/FRAC-1)) {
          frac.order.stats[i,] <- frac.order.stats[i-1,] + 
            qbeta(unirands[i,],FRAC,n+1-i*FRAC)*(1-frac.order.stats[i-1,])
        }
      }
    }
    if (VERBOSE) {
      cat(sprintf("n=%d Dirichlet draw time elapsed: %s",n,format(Sys.time()-start.time)), 
          sep="\n",append=TRUE,file="")
      cat(sprintf("n=%d Dirichlet draw time elapsed: %s",n,format(Sys.time()-start.time)), 
          sep="\n",append=TRUE,file=LOGFILE)
    }
    
    # Loop over alpha values
    for (ALPHA in ALPHAS) {
      start.time <- Sys.time()
      indleft <- 1:N.DRAWS
      CALIB.ERR <- ALPHA*CALIB.DEC.ERR
      TOL.ALPHA <- CALIB.ERR + qnorm(CALIB.ERR.PERC)*sqrt((ALPHA+CALIB.ERR)*(1-ALPHA-CALIB.ERR)) / sqrt(N.DRAWS) #see paper; rem: qnorm(.05)<0
      # Binary search (algorithm could be improved)
      a.lo <- 0;  f.lo <- 0
      a.hi <- 1;  f.hi <- 1
      rejs <- rejs.known <- 0
      f.prev <- Inf
      a.cur <- a.prev <- prev.alpha.tilde[which(ALPHAS==ALPHA)]
      while (abs(f.prev-ALPHA)>TOL.ALPHA) {
        CIs <- t(GK.dist.inf.1s.CIs(n=n,alpha.tilde=a.cur,FRAC=FRAC,tt=TT))
        rejvec <- rep.int(NA,length(indleft))
        for (indind in 1:length(indleft)) {
          ind <- indleft[indind]
          rejvec[indind] <- 1 - prod(CIs[1,]<frac.order.stats[1:n,ind])*prod(frac.order.stats[1:n,ind]<CIs[2,])
        }
        rejs <- rejs.known + sum(rejvec)
        a.prev <- a.cur;  f.prev <- rejs/N.DRAWS
        if (f.prev>ALPHA) {
          a.hi <- a.cur;  f.hi <- f.prev
          if (a.lo>0) {
            a.cur <- (a.lo+a.cur)/2
          } else {
            a.cur <- a.cur / 1.4
          }
          indleft <- indleft[as.logical(rejvec)]
        } else { #a.cur is too small
          a.lo <- a.cur;  f.lo <- f.prev
          if (a.hi<1) {
            a.cur <- (a.hi+a.cur)/2
          } else {
            a.cur <- min((1+a.cur)/2, a.cur * 1.4)
          }
          rejs.known <- rejs #increased
          indleft <- indleft[!as.logical(rejvec)]
        }
        #
        if (VERBOSE) {
          cat(sprintf(sprintf("f.prev=%s,a.prev=%s",SIGFMTf,sprintf("%%%d.%df",SIGFIGa+3,SIGFIGa+1)),f.prev,a.prev),
              file="",append=T,sep="\n")
          cat(sprintf(sprintf("f.prev=%s,a.prev=%s",SIGFMTf,sprintf("%%%d.%df",SIGFIGa+3,SIGFIGa+1)),f.prev,a.prev),
              file=LOGFILE,append=T,sep="\n")
        }
        if (round(a.hi,digits=SIGFIGa)==round(a.lo,digits=SIGFIGa)) {
          stop(sprintf(sprintf("{a.lo,a.hi,f.lo,f.hi}={%s,%s,%s,%s}",SIGFMTa,SIGFMTa,SIGFMTf,SIGFMTf),
                       a.lo,a.hi,f.lo,f.hi))
        }
      }
      alpha.tilde <- a.prev
      prev.alpha.tilde[which(ALPHAS==ALPHA)] <- alpha.tilde
      
      # Output results
      cat(sprintf(sprintf("%%5.3f,%%d,%s,%s",SIGFMTa,SIGFMTf),
                  ALPHA, n, alpha.tilde, f.prev),
          file=OUTFILE,sep="\n",append=TRUE)
      if (VERBOSE) {
        cat(sprintf(sprintf("alphatilde(alpha=%%5.3f,n=%%d)=%s;f.prev=%s",SIGFMTa,SIGFMTf),
                    ALPHA, n, alpha.tilde, f.prev),
            file="",sep="\n",append=TRUE)
        cat(sprintf("Time elapsed: %s",format(Sys.time()-start.time)), 
            sep="\n",append=TRUE,file="")
        cat(sprintf(sprintf("alphatilde(alpha=%%5.3f,n=%%d)=%s;f.prev=%s",SIGFMTa,SIGFMTf),
                    ALPHA, n, alpha.tilde, f.prev),
            file=LOGFILE,sep="\n",append=TRUE)
        cat(sprintf("Time elapsed: %s",format(Sys.time()-start.time)), 
            sep="\n",append=TRUE,file=LOGFILE)
      }
    } #ALPHAS loop
    prev.n <- n
  } #NS.SORT loop
  
  if (!is.null(oldseed)) {
    assign(".Random.seed", oldseed, .GlobalEnv)
  }
  return(TRUE) #all went smoothly
}



#
# More two-sample helper functions
#

#Function: look up (or approximate) calibrated alpha.tilde from 2-sample table,
# given sample sizes and desired alpha.
GK.dist.2s.alpha.tilde <- function(nx,ny,alpha) {
  exact.match <- any(GK.dist.inf.2s.lookup$ALPHA==alpha & GK.dist.inf.2s.lookup$nx==nx & GK.dist.inf.2s.lookup$ny==ny)
  if (!exact.match && min(nx/ny,ny/nx)<1/10) {
    warning("Values of alpha.tilde for nx/ny<0.1 have not been thoroughly explored; it is recommended that you add your case to the lookup table using GK.dist.inf.2s.calibrate() or at least verify the suggested alpha.tilde using GK.dist.inf.2s.verify().")
  }
  # Apply formula, if applicable; conservative unless nx=ny, less accurate when n<200
  at.formula <- NULL
  if (alpha==0.01) {
    at.formula <- exp(0.02576944-3.37109911*sqrt(log(log(max(nx,ny)))))
  } else if (alpha==0.05) {
    at.formula <- exp(0.6872947-3.1503102*sqrt(log(log(max(nx,ny)))))
  } else if (alpha==0.10) {
    at.formula <- exp(1.081619-3.125037*sqrt(log(log(max(nx,ny)))))
  }
  if (!is.null(at.formula) && nx==ny && nx>=200) { return(c(at.formula,NA)) }
  #
  lookup.ALPHA <- GK.dist.inf.2s.lookup[which(GK.dist.inf.2s.lookup$ALPHA==alpha),]
  if (dim(lookup.ALPHA)[1]==0) {
    stop(sprintf("alpha must be in {%s} to use the lookup table, but the argument's value is %g.  Use argument alpha=NULL instead, to get a p-value, and compare it to alpha; reject H0 if the p-value is below alpha.",
                 paste0(unique(GK.dist.inf.2s.lookup$ALPHA),collapse=","),alpha))
  }
  # exact match or else err on conservative side
  lookup.conservative <- lookup.ALPHA[(lookup.ALPHA$nx>=nx)&(lookup.ALPHA$ny>=ny),] #each alpha-tilde in here is no bigger than exact alpha-tilde, since alpha-tilde decreases with sample size.
  if (dim(lookup.conservative)[1]==0) { #nx,ny too big for table--use formula or throw error
    if (!is.null(at.formula)) {
      return(c(at.formula,NA))
    } else {
      stop(sprintf("No row in the lookup table exists with nx>=%d and ny>=%d for alpha=%g, and the alpha.tilde formula is only calibrated for alpha in (0.01,0.05,0.10).  You may run GK.dist.inf with alpha=NULL to get a p-value (which can also be used for hypothesis testing), or use the function GK.dist.inf.2s.calibrate to add the desired row to the lookup table in order to construct uniform confidence bands, in which case make sure to source('GK_dist_inf.R') again to reload the table.",nx,ny,alpha))
    }
  }
  ind <- which(lookup.conservative$alpha.tilde==max(lookup.conservative$alpha.tilde)) #take least conservative among these, which may match nx and ny exactly (if such an entry in the lookup exists).
  lookup.ret <- lookup.conservative[ind,]
  ret.alpha.tilde <- lookup.ret$alpha.tilde[1] #  In the unlikely event of multiple rows with identical alpha-tilde, the min index is chosen--but it's the same alpha-tilde anyway, so results are unaffected by the choice.
  ret.alpha.sim.lo <- mean(lookup.ret$ALPHA.simulated.lo) #if multiple rows w/ same alpha-tilde, average their simulated alpha values.
  ret.alpha.sim.hi <- mean(lookup.ret$ALPHA.simulated.hi)
  if (ret.alpha.sim.lo!=ret.alpha.sim.hi) {
    ret.alpha.tilde <- ret.alpha.tilde - 0.5*10^(-SIGFIGa) #make a little smaller when there's a discontinuity in alpha as a function of alpha-tilde around the current value
  }
  # Finally...just use the formula if significantly less conservative (and n>=200)
  if (!is.null(at.formula)) {
    if (max(nx,ny)>=200 && at.formula>1.01*ret.alpha.tilde) {
      return(c(at.formula,NA))
    }
  }
  return(c(ret.alpha.tilde, ret.alpha.sim.lo))
}

# Function (subroutine): does 2-sample test reject global null hypothesis?
# pre-condition: U[[1]] and U[[2]] are numeric vectors
#                of all distinct/unique values (fine for simulated data but not all real data samples)
#                sorted in increasing order.
GK.dist.inf.rej.2s <- function(n,L,U,H,ONESIDED=0) {
  ind <- c(1,1)
  xy <- 1 #1=X, 2=Y
  U[[1]] <- c(U[[1]],Inf);  U[[2]] <- c(U[[2]],Inf)
  if (U[[1]][1] < U[[2]][1]) { xy <- 2 }
  #
  while (ind[xy]<=n[xy]) {
    xy <- 3-xy
    while (U[[xy]][ind[xy]+1] < U[[3-xy]][ind[3-xy]]) {
      ind[xy] <- ind[xy]+1
    }
    if (L[[xy]][ind[xy]] > H[[3-xy]][ind[3-xy]]) {
      if (abs(ONESIDED)*abs((3-ONESIDED)/2-xy)<0.001) return(1)
    }
    ind[xy] <- ind[xy]+1
  }
  return(0) #if made it this far, don't reject
}

# Function (subroutine): MTP rejections by sample value.
# argument GK.dist.inf.obj is return value from two-sample call to GK.dist.inf()
# return value: data.frame showing whether H_{0r} is rejected (reject=TRUE) or not (reject=FALSE) for each interval [from,to)
GK.dist.inf.rej.2s.r <- function(GK.dist.inf.obj) {
  obj <- GK.dist.inf.obj
  # consolidate any duplicate Ux or Uy, pull back "next" upper CI endpoint to mimic stairstep band
  dfX <- data.frame(from=c(-Inf, obj$U$x, Inf),
                    CIx.lo=c(0, obj$CIs.x[,1], 1),
                    CIx.hi=c(0, obj$CIs.x[,2], 1) )
  dfY <- data.frame(from=c(-Inf, obj$U$y, Inf),
                    CIy.lo=c(0, obj$CIs.y[,1], 1),
                    CIy.hi=c(0, obj$CIs.y[,2], 1) )
  dfX <- dfX[order(dfX$from, dfX$CIx.lo, decreasing=FALSE) , ]
  dfY <- dfY[order(dfY$from, dfY$CIy.lo, decreasing=FALSE) , ]
  # code for X
  excl.inds <- NULL # indices to exclude (duplicate 'from' value)
  for (i in 1:(nrow(dfX)-1)) {
    dfX$CIx.hi[i] <- dfX$CIx.hi[i+1] # stairstep band: take upper endpoint from *next* CI
    if (dfX$from[i]==dfX$from[i+1]) excl.inds <- c(excl.inds, i)
  }
  dfX <- dfX[-excl.inds , ]
  # code for Y
  excl.inds <- NULL # indices to exclude (duplicate 'from' value)
  for (i in 1:(nrow(dfY)-1)) {
    dfY$CIy.hi[i] <- dfY$CIy.hi[i+1] # stairstep band: take upper endpoint from *next* CI
    if (dfY$from[i]==dfY$from[i+1]) excl.inds <- c(excl.inds, i)
  }
  dfY <- dfY[-excl.inds , ]
  # combine X and Y
  dfX$CIy.lo <- dfX$CIy.hi <- NA
  dfY$CIx.lo <- dfY$CIx.hi <- NA
  df <- rbind(dfX, dfY)
  df <- df[order(df$from, df$CIx.lo) , ] # first duplicate row has non-NA CIx, NA CIy
  df$to <- NA # for the ranges of values on which to assess rejection
  df <- df[,c('from','to','CIx.lo','CIx.hi','CIy.lo','CIy.hi')]
  # deal with duplicates and add 'to' for ranges
  excl.inds <- NULL
  for (i in 2:(nrow(df)-1)) {
    df$to[i] <- df$from[i+1]
    if (is.na(df$CIx.lo[i])) {
      df$CIx.lo[i] <- df$CIx.lo[i-1]
      df$CIx.hi[i] <- df$CIx.hi[i-1]
    } else {
      df$CIy.lo[i] <- df$CIy.lo[i-1]
      df$CIy.hi[i] <- df$CIy.hi[i-1]
    }
    if (df$from[i]==df$from[i-1]) excl.inds <- c(excl.inds, i-1)
  }
  df <- df[-excl.inds,]
  df <- df[!(df$from==Inf) , ]
  df$reject <- (df$CIx.lo > df$CIy.hi) | (df$CIx.hi < df$CIy.lo)
  return(df[,c('from','to','reject','CIx.lo','CIx.hi','CIy.lo','CIy.hi')])
}
# If you need to compare: results from before 24oct2024 used the following
OLD.GK.dist.inf.rej.2s.r <- function(GK.dist.inf.obj) {
  obj <- GK.dist.inf.obj
  Ux <- obj$U$x;  Uy <- obj$U$y
  CIx <- rbind(c(-Inf,NA), obj$CIs.x, c(NA,Inf))
  CIy <- rbind(c(-Inf,NA), obj$CIs.y, c(NA,Inf))
  r <- data.frame( values=c(Ux, Uy), src=I(c(rep('x',length(Ux)), rep('y',length(Uy)))))
  tmp <- order(r[,1])
  r <- r[tmp,]
  r$xind <- cumsum(r[,2]=='x');  r$yind <- cumsum(r[,2]=='y')
  rej <- rep(NA,dim(r)[1])
  for (i in 1:length(rej)) {
    #"extra" +1 b/c insert c(-Inf,NA) above
    tmp1 <- CIx[r[i,3]+1,1]   > CIy[r[i,4]+1+1,2]
    tmp2 <- CIx[r[i,3]+1+1,2] < CIy[r[i,4]+1,1]
    rej[i] <- tmp1 || tmp2
  }
  return(data.frame(from=c(-Inf,r[,1]), to=c(r[,1],Inf), rej=c(F,rej)))
}

# Simulate (joint) type I error given alpha-tilde(s), 2-sample
#PARALLEL: # CPUs available (requires foreach and doParallel packages)
GK.dist.inf.2s.verify <- function(alpha.tildes=NULL,nx=NULL,ny=NULL,VERBOSE=FALSE,MIN.DRAWS=2e5,PARALLEL=1) {
  if (PARALLEL>1) {
    if (!require("parallel") || !require("foreach")) {warning("Install package foreach in order to run PARALLEL."); PARALLEL <- 0}
    if (!require("doParallel")) {warning("Install package doParallel in order to run PARALLEL."); PARALLEL <- 0}
    PARALLEL <- tryCatch({workers <- makeCluster(PARALLEL); registerDoParallel(workers); on.exit(stopCluster(workers),add=TRUE); PARALLEL}, error=function(Err){warning(sprintf("Error creating %d clusters",PARALLEL)); 0})
  } else if (PARALLEL<0) {
    warning("PARALLEL must be a non-negative integer: 0 (or 1) for not parallel, positive for number of parallel CPUs to use.")
  }
  overall.start.time <- Sys.time()
  oldseed <- NULL
  if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
    oldseed <- get(".Random.seed",.GlobalEnv)
  }
  on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
  set.seed(112358) #for replicability
  # Constant parameters
  BLOCKSIZE <- max(1,round(1e6/(nx+ny)))
  N.DRAWS <- BLOCKSIZE*ceiling(MIN.DRAWS/BLOCKSIZE) #2*10^5 #X*10^5; see paper
  FRAC <- 1 
  TT.x <- array(1/2,dim=c(ceiling((nx+1)/FRAC-1),2)) 
  TT.y <- array(1/2,dim=c(ceiling((ny+1)/FRAC-1),2)) 
  # Storage variables
  rejs <- matrix(NA,N.DRAWS,length(alpha.tildes))
  CIs.x <- vector("list",length(alpha.tildes))
  CIs.y <- vector("list",length(alpha.tildes))
  for (i in 1:length(alpha.tildes)) {
    CIs.x[[i]] <- GK.dist.inf.1s.CIs(n=nx,alpha.tilde=alpha.tildes[i],FRAC=FRAC,tt=TT.x)
    CIs.y[[i]] <- GK.dist.inf.1s.CIs(n=ny,alpha.tilde=alpha.tildes[i],FRAC=FRAC,tt=TT.y)
  }
  if (PARALLEL>1) {
    clusterSetRNGStream(workers,112358)
    rejs <- foreach(irep=1:(N.DRAWS/BLOCKSIZE),.combine=rbind,.inorder=FALSE,.export=c("GK.dist.inf.rej.2s")) %dopar% {
      unirands.x <- matrix(runif(BLOCKSIZE*ceiling((nx+1)/FRAC-1)),BLOCKSIZE)
      unirands.y <- matrix(runif(BLOCKSIZE*ceiling((ny+1)/FRAC-1)),BLOCKSIZE)
      frac.order.stats.x <- frac.order.stats.y <- NULL
      if (FRAC==1) {
        frac.order.stats.x <- t(apply(unirands.x,1,sort))
        frac.order.stats.y <- t(apply(unirands.y,1,sort))
      } else {
        # Draw from Dirichlet for (fractional) uniform order statistics
        frac.order.stats.x <- array(0,dim=c(BLOCKSIZE,ceiling((nx+1)/FRAC-1)))
        frac.order.stats.x[,1] <- qbeta(unirands.x[,1],FRAC,nx+1-FRAC) 
        for (i in 2:ceiling((nx+1)/FRAC-1)) {
          frac.order.stats.x[,i] <- frac.order.stats.x[,i-1] + 
            qbeta(unirands.x[,i],FRAC,nx+1-i*FRAC)*(1-frac.order.stats.x[,i-1])
        }
        frac.order.stats.y <- array(0,dim=c(BLOCKSIZE,ceiling((ny+1)/FRAC-1)))
        frac.order.stats.y[,1] <- qbeta(unirands.y[,1],FRAC,ny+1-FRAC) 
        for (i in 2:ceiling((ny+1)/FRAC-1)) {
          frac.order.stats.y[,i] <- frac.order.stats.y[,i-1] + 
            qbeta(unirands.y[,i],FRAC,ny+1-i*FRAC)*(1-frac.order.stats.y[,i-1])
        }
      }
      #
      reprejs <- matrix(NA,BLOCKSIZE,length(alpha.tildes))
      L <- H <- vector("list",length(alpha.tildes))
      for (aind in 1:length(alpha.tildes)) {
        L[[aind]] <- list(xlo=CIs.x[[aind]][,1], ylo=CIs.y[[aind]][,1])
        H[[aind]] <- list(xhi=CIs.x[[aind]][,2], yhi=CIs.y[[aind]][,2])
      }
      for (aind in 1:length(alpha.tildes)) {
        for (j in 1:BLOCKSIZE) {
          #           rejs[BLOCKSIZE*(irep-1)+j,aind] <- 
          reprejs[j,aind] <- 
            GK.dist.inf.rej.2s(n=c(nx,ny), L=L[[aind]], U=list(Ux=frac.order.stats.x[j,],Uy=frac.order.stats.y[j,]), H=H[[aind]])
        }
      }
      reprejs
    }
  } else {
    for (irep in 1:(N.DRAWS/BLOCKSIZE)) {
      if (VERBOSE) {
        cat(sprintf("Start block %d/%d; ",irep,N.DRAWS/BLOCKSIZE),file="",sep=';')
      }
      unirands.x <- matrix(runif(BLOCKSIZE*ceiling((nx+1)/FRAC-1)),BLOCKSIZE)
      unirands.y <- matrix(runif(BLOCKSIZE*ceiling((ny+1)/FRAC-1)),BLOCKSIZE)
      frac.order.stats.x <- array(0,dim=c(BLOCKSIZE,ceiling((nx+1)/FRAC-1)))
      frac.order.stats.y <- array(0,dim=c(BLOCKSIZE,ceiling((ny+1)/FRAC-1)))
      # Draw from Dirichlet for (fractional) uniform order statistics
      start.time <- Sys.time()
      frac.order.stats.x[,1] <- qbeta(unirands.x[,1],FRAC,nx+1-FRAC) 
      for (i in 2:ceiling((nx+1)/FRAC-1)) {
        frac.order.stats.x[,i] <- frac.order.stats.x[,i-1] + 
          qbeta(unirands.x[,i],FRAC,nx+1-i*FRAC)*(1-frac.order.stats.x[,i-1])
      }
      if (VERBOSE) {
        cat(sprintf("nx=%d Dirichlet draw time elapsed: %s; ",nx,format(Sys.time()-start.time)), 
            sep=";",append=TRUE,file="")
      }
      #
      start.time <- Sys.time()
      frac.order.stats.y[,1] <- qbeta(unirands.y[,1],FRAC,ny+1-FRAC) 
      for (i in 2:ceiling((ny+1)/FRAC-1)) {
        frac.order.stats.y[,i] <- frac.order.stats.y[,i-1] + 
          qbeta(unirands.y[,i],FRAC,ny+1-i*FRAC)*(1-frac.order.stats.y[,i-1])
      }
      if (VERBOSE) {
        cat(sprintf("ny=%d Dirichlet draw time elapsed: %s; ",ny,format(Sys.time()-start.time)), 
            sep=";",append=TRUE,file="")
      }
      #
      L <- H <- vector("list",length(alpha.tildes))
      for (aind in 1:length(alpha.tildes)) {
        L[[aind]] <- list(xlo=CIs.x[[aind]][,1], ylo=CIs.y[[aind]][,1])
        H[[aind]] <- list(xhi=CIs.x[[aind]][,2], yhi=CIs.y[[aind]][,2])
      }
      for (aind in 1:length(alpha.tildes)) {
        for (j in 1:BLOCKSIZE) {
          rejs[BLOCKSIZE*(irep-1)+j,aind] <- 
            GK.dist.inf.rej.2s(n=c(nx,ny), L=L[[aind]], U=list(Ux=frac.order.stats.x[j,],Uy=frac.order.stats.y[j,]), H=H[[aind]])
        }
      }
    }
  }
  RP0s <- colMeans(rejs)
  return(list(RP0s=RP0s,n=c(nx,ny),N.DRAWS=N.DRAWS,alpha.tildes=alpha.tildes,
              time.elapsed=format(Sys.time()-overall.start.time),FRAC=FRAC,TT.x=TT.x,TT.y=TT.y))
}

# Plotting function for two-sample bands.  dist.inf.obj is the returned list from GK.dist.inf.2s
# NOTE: you may need to open device (X11, quartz, etc.) before calling this function.
GK.dist.inf.plot.2s <- function(dist.inf.obj=NULL,main="Bands determining rejection\n(NOT 'confidence bands')",sub="",xlab="X, Y",ylab="F(X), F(Y)",
                                pchX=3,pchY=4,ltyX="31",ltyY="22",lwdX=2,lwdY=2,
                                draw.legend=TRUE,legend=c("X","Y")) {
  if (is.null(dist.inf.obj)) {stop("Argument should be returned list from GK.dist.inf.2s")}
  U <- dist.inf.obj$U
  n <- dist.inf.obj$n
  CIs.x <- dist.inf.obj$CIs.x
  CIs.y <- dist.inf.obj$CIs.y
  n <- c(length(U[[1]]),length(U[[2]]))
  tmp.big <- 100*max(abs(unlist(U)))
  plot.band.x.L <- rbind(c(U[[1]][1]-tmp.big,0),
                         cbind(c(rep(U[[1]],each=2),U[[1]][n[1]]),
                               c(0,rep(CIs.x[,1],each=2))),
                         c(U[[1]][n[1]]+tmp.big,CIs.x[n[1],1]))
  plot.band.y.L <- rbind(c(U[[2]][1]-tmp.big,0),
                         cbind(c(rep(U[[2]],each=2),U[[2]][n[2]]),
                               c(0,rep(CIs.y[,1],each=2))),
                         c(U[[2]][n[2]]+tmp.big,CIs.y[n[2],1]))
  plot.band.x.H <- rbind(c(U[[1]][1]-tmp.big,CIs.x[1,2]),
                         cbind(c(U[[1]][1],rep(U[[1]],each=2)),
                               c(rep(CIs.x[,2],each=2),1)),
                         c(U[[1]][n[1]]+tmp.big,1))
  plot.band.y.H <- rbind(c(U[[2]][1]-tmp.big,CIs.y[1,2]),
                         cbind(c(U[[2]][1],rep(U[[2]],each=2)),
                               c(rep(CIs.y[,2],each=2),1)),
                         c(U[[2]][n[2]]+tmp.big,1))
  par(family="serif",mar=c(5.0,6.0,6.0,2.1))
  plot(x=plot.band.x.L[,1],y=plot.band.x.L[,2],type="n",pch=1,col=1, main=main, sub=sub, xlab=xlab, ylab=ylab, cex.main=2.0, mgp=c(3,1,0), cex.axis=2, cex.lab=2, xlim=c(min(unlist(U)),max(unlist(U))), ylim=c(0,1))
  if (!all(CIs.x[,1]<=0)) {
    lines(plot.band.x.L, type='l', col=1, lwd=lwdX, lty=ltyX)
    points(U[[1]],CIs.x[,1], type='p', col=1, pch=pchX, lwd=lwdX)
  }
  if (!all(CIs.x[,2]>=1)) {
    lines(plot.band.x.H, type='l', col=1, lwd=lwdX, lty=ltyX)
    points(U[[1]],CIs.x[,2], type='p', col=1, pch=pchX, lwd=lwdX)
  }
  if (!all(CIs.y[,1]<=0)) {
    lines(plot.band.y.L, type='l', col=1, lwd=lwdY, lty=ltyY)
    points(U[[2]],CIs.y[,1], type='p', col=1, pch=pchY, lwd=lwdY)
  }
  if (!all(CIs.y[,2]>=1)) {
    lines(plot.band.y.H, type='l', col=1, lwd=lwdY, lty=ltyY)
    points(U[[2]],CIs.y[,2], type='p', col=1, pch=pchY, lwd=lwdY)
  }
  if (draw.legend) {
    legend('bottomright', legend=legend, inset=0.01, col=1, pch=c(pchX,pchY), lty=c(ltyX,ltyY), lwd=c(lwdX,lwdY), cex=1.8, y.intersp=0.9, bty='n') #bg='white'
  }
  return(list(band.XL=plot.band.x.L,band.XH=plot.band.x.H,
              band.YL=plot.band.y.L,band.YH=plot.band.y.H))
}

#
# Compute alpha-tilde entry and append to lookup table, 2-sample
#
# NS: first sample's size is in first row, while second row is for second sample's size; iterates through columns.
# prev.alpha.tilde: vector of same length as ALPHAS, containing best guess of alpha.tilde for first iteration; NULL is fine.
# PARALLEL: number of CPUs available--requires foreach and doParallel packages.  For now, only distributes different ALPHAS to different CPUs, so no effect if length(ALPHAS)=1.
# other args: see comments before GK.dist.inf.1s.calibrate
GK.dist.inf.2s.calibrate <- function(NS=NULL,ALPHAS=NULL,N.DRAWS=NULL,CALIB.DEC.ERR=NULL,CALIB.ERR.PERC=0.05,OUTFILE="GK_dist_inf_2s_lookup.txt",LOGFILE="GK_dist_inf_2s_lookup_log.txt",VERBOSE=TRUE,prev.alpha.tilde=NULL,PARALLEL=1) {
  FRAC <- 1 #only implemented for FRAC=1, but in principle could use FRAC=1/k for some integer k for smaller n, or FRAC=k for larger n
  SIGFMTa <- sprintf("%%%d.%df",SIGFIGa+2,SIGFIGa)
  SIGFMTf <- SIGFMTa
  if (is.null(NS) || is.null(ALPHAS)) {stop("Error in GK.dist.inf.2s.calibrate: NULL value for NS and/or ALPHAS not allowed.")}
  if (is.null(CALIB.DEC.ERR)) {
    if (min(ALPHAS)>0.01) {CALIB.DEC.ERR <- 0.02} else {CALIB.DEC.ERR <- 0.05}
  }
  if (is.null(N.DRAWS)) {
    if (min(ALPHAS)>0.01) {N.DRAWS <- 2e5} else {N.DRAWS <- 1e6}
  }
  if (PARALLEL>1) {
    if (!require("parallel") || !require("foreach")) {warning("Install package foreach in order to run PARALLEL."); PARALLEL <- 1}
    if (!require("doParallel")) {warning("Install package doParallel in order to run PARALLEL."); PARALLEL <- 1}
    PARALLEL <- tryCatch({workers <- makeCluster(PARALLEL); registerDoParallel(workers); on.exit(stopCluster(workers),add=TRUE); PARALLEL}, 
                         error=function(Err){warning(sprintf("Error creating %d clusters",PARALLEL)); 1})
  } else if (PARALLEL<1) {
    warning("PARALLEL must be a positive integer: 1 for not parallel, positive for number of parallel CPUs to use.")
  }
  if (PARALLEL>1) {
    iblkcol.perm <- function(dat, cols, ...) {
      i <- 1
      it <- idiv(length(cols), ...) #idiv(ncol(dat), ...)
      nextEl <- function() {
        n <- nextElem(it)
        r <- seq(i, length=n)
        i <<- i + n
        dat[,cols[r],drop=FALSE] #dat[,r, drop=FALSE]
      }
      obj <- list(nextElem=nextEl)
      class(obj) <- c('abstractiter', 'iter')
      obj
    }
  }
  oldseed <- NULL
  if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
    oldseed <- get(".Random.seed",.GlobalEnv)
  }
  on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
  
  if (VERBOSE) {
    cat(sprintf("N.DRAWS=%d; CALIB.DEC.ERR=%g; CALIB.ERR.PERC=%g; ALPHAS={%s}; NXS={%s}; NYS={%s}; OUTFILE=%s",
                N.DRAWS, CALIB.DEC.ERR, CALIB.ERR.PERC, paste0(sprintf("%g",ALPHAS),collapse=","),
                paste0(sprintf("%g",NS[1,]),collapse=","), paste0(sprintf("%g",NS[2,]),collapse=","),
                OUTFILE), file="",sep="\n")
  }
  
  #Appends to end of specified file now, not overwrite (unless file doesn't exist)
  if (OUTFILE!="") {
    if (!file.exists(OUTFILE)) {  # Write to file: header w/ variable names
      cat(sprintf("Creating new file: %s",OUTFILE), file="",sep="\n")
      cat('"ALPHA","nx","ny","alpha.tilde","ALPHA.simulated.lo","ALPHA.simulated.hi"',file=OUTFILE,sep="\n",append=FALSE)
    } else { # append to existing file
      cat(sprintf("Appending to existing file: %s",OUTFILE), file="",sep="\n")
    }
  } else {
    cat("No file specified, so output is to console only.", file="",sep="\n")
  }
  
  NS <- apply(X=NS,MARGIN=2,FUN=sort) #w/in each column, top row n smaller than bottom.
  if (is.null(prev.alpha.tilde)) {
    prev.alpha.tilde <- ALPHAS  #bad guess
    if (!is.null(GK.dist.inf.2s.lookup)) {
      for (k in 1:length(ALPHAS)) {
        tryCatch({
          lka <- GK.dist.inf.2s.lookup[GK.dist.inf.2s.lookup$ALPHA==ALPHAS[k],]
          tmpy <- sort(unique(lka$ny))
          tmpind <- which.min(abs(NS[2,1]-tmpy))
          lkay <- lka[lka$ny==tmpy[tmpind],]
          tmpx <- sort(unique(lkay$nx))
          tmpind <- which.min(abs(NS[1,1]-tmpx))
          prev.alpha.tilde[k] <- lkay[lkay$nx==tmpx[tmpind],]$alpha.tilde
        }, error=function(ERR) {})
      }
    }
  }
  for (nind in 1:dim(NS)[2]) {
    nx <- NS[1,nind];  ny <- NS[2,nind]
    TT.x <- array(1/2,dim=c(ceiling((nx+1)/FRAC-1),2)) #equal type I error @ all quantiles, all equal-tailed
    TT.y <- array(1/2,dim=c(ceiling((ny+1)/FRAC-1),2))
    permleft <- 1:N.DRAWS #could be significant memory--do before below
    rej.ind <- rep.int(as.integer(NA),N.DRAWS) #could be significant memory
    start.time <- Sys.time()
    REDRAW <- FALSE
    PREFILL <- FALSE
    if (PARALLEL>1) {
      NBLOCKS <- max(100,8*PARALLEL)
      BLOCKSIZE <- ceiling(N.DRAWS/NBLOCKS)
      PREFILL <- tryCatch({
        clusterSetRNGStream(workers,112358)
        N.DRAWS <- NBLOCKS*BLOCKSIZE
        permdata <- foreach(i=1:NBLOCKS,.combine=cbind,.inorder=FALSE) %dopar% {
          replicate(n=BLOCKSIZE,expr=sample.int(nx+ny),simplify='array')  
        }
        TRUE}, error=function(Err) FALSE )
    }
    if (!PREFILL) {
      set.seed(112358) #for replicability
      gc(verbose=FALSE)
      if (0.95*memory.limit() < 4*N.DRAWS*(nx+ny)/2^20 + memory.size()) {
        problem <- TRUE
      } else if (0.4*memory.limit() < 4*N.DRAWS*(nx+ny)/2^20 + memory.size()) {
        problem <- tryCatch({
          permdata <- matrix(as.integer(NA),(nx+ny),N.DRAWS)
          for (j in 1:N.DRAWS) { permdata[,j] <- sample.int(nx+ny) }
          FALSE}, 
          error=function(Errr)TRUE)
      } else {
        problem <- tryCatch({permdata <- replicate(n=N.DRAWS,expr=sample.int(nx+ny),simplify='array');FALSE}, error=function(Errr)TRUE)
      }
      if (problem) { #need to re-draw sim'd data each time
        REDRAW <- TRUE; warning("Not enough memory to store all simulated data; will be slower."); if (VERBOSE) {cat("Not enough memory to store all simulated data; will be slower.",sep='\n',append=TRUE,file="");cat("Not enough memory to store all simulated data; will be slower.",sep='\n',append=TRUE,file=LOGFILE)} 
      } else {
        if (VERBOSE) {
          cat(sprintf("(nx,ny)=(%d,%d) permutation draw time elapsed: %s",nx,ny,format(Sys.time()-start.time)), 
              sep="\n",append=TRUE,file="")
          cat(sprintf("(nx,ny)=(%d,%d) permutation draw time elapsed: %s",nx,ny,format(Sys.time()-start.time)), 
              sep="\n",append=TRUE,file=LOGFILE)
        }
      }
    }
    
    # Loop over alpha values
    for (ALPHA in ALPHAS) {
      permleft <- 1:N.DRAWS
      rej.ind <- rep.int(as.integer(NA),N.DRAWS) #if we're out of memory, let's know now
      CALIB.ERR <- ALPHA*CALIB.DEC.ERR
      TOL.ALPHA <- CALIB.ERR + qnorm(CALIB.ERR.PERC)*sqrt((ALPHA+CALIB.ERR)*(1-ALPHA-CALIB.ERR)) / sqrt(N.DRAWS) #see paper; rem: qnorm(.05)<0
      start.time <- Sys.time()
      # Binary search
      L <- U <- H <- vector("list",2) #L(ow) and H(igh) CI endpoints; U(niform) order stats
      a.cur <- a.prev <- prev.alpha.tilde[which(ALPHAS==ALPHA)]
      a.lo <- 0;  a.hi <- 1;  f.lo <- 0;  f.hi <- 1
      rejs <- rejs.known <- 0
      rejs.TARGET <- ALPHA*N.DRAWS;  rejs.TOL <- TOL.ALPHA*N.DRAWS
      while ((abs(rejs-rejs.TARGET)>rejs.TOL) && (round(a.hi,digits=SIGFIGa)!=round(a.lo,digits=SIGFIGa))) {
        START.VERIFY.TIME <- Sys.time()
        if (REDRAW) {
          f.prev <- GK.dist.inf.2s.verify(alpha.tildes=a.cur,nx=nx,ny=ny,VERBOSE=VERBOSE,MIN.DRAWS=N.DRAWS,PARALLEL=PARALLEL)$RP0s
          rejs <- N.DRAWS*f.prev
        } else {
          CIs.x <- GK.dist.inf.1s.CIs(n=nx,alpha.tilde=a.cur,FRAC=FRAC,tt=TT.x)
          CIs.y <- GK.dist.inf.1s.CIs(n=ny,alpha.tilde=a.cur,FRAC=FRAC,tt=TT.y)
          L[[1]] <- CIs.x[,1];  L[[2]] <- CIs.y[,1]
          H[[1]] <- CIs.x[,2];  H[[2]] <- CIs.y[,2]
          # compute rejections given (fixed) CIs and (simulated) data
          PAR.SUCCESS <- FALSE
          if (PARALLEL>1 && length(permleft)>500) {
            PAR.SUCCESS <- tryCatch({
              rej.ind <- foreach(samps=iblkcol.perm(dat=permdata,cols=permleft,chunkSize=max(1,ceiling(1e4/(nx+ny)))), .combine=c, .inorder=TRUE, .export=c("GK.dist.inf.rej.2s")) %dopar% {
                itret <- rep.int(as.integer(NA),dim(samps)[2])
                for (j in 1:dim(samps)[2]) {
                  itret[j] <- GK.dist.inf.rej.2s(n=c(nx,ny),L=L,U=list(x=sort(samps[1:nx,j]),y=sort(samps[(nx+1):(nx+ny),j])),H=H)
                }
                as.integer(itret)
              }
              TRUE}, error=function(Err) FALSE)
          }
          if (!PAR.SUCCESS) {
            rej.ind <- rep.int(as.integer(NA),length(permleft))
            for (j in 1:length(permleft)) {
              permind <- permleft[j]
              U[[1]] <- sort(permdata[1:nx,permind])
              U[[2]] <- sort(permdata[(nx+1):(nx+ny),permind])
              rej.ind[j] <- GK.dist.inf.rej.2s(n=c(nx,ny),L=L,U=U,H=H)
            }
          }
          rejs <- rejs.known + sum(rej.ind)
          f.prev <- rejs/N.DRAWS
        }
        a.prev <- a.cur
        if (rejs>rejs.TARGET) {
          a.hi <- a.cur;  f.hi <- f.prev
          if (a.lo>0) {
            a.cur <- (a.lo+a.cur)/2
          } else {
            a.cur <- a.cur / 1.1
          }
          if (!REDRAW) {
            permleft <- permleft[as.logical(rej.ind)]
          }
        } else { #a.cur is too small
          a.lo <- a.cur;  f.lo <- f.prev
          if (a.hi<1) {
            a.cur <- (a.hi+a.cur)/2
          } else {
            a.cur <- a.cur * 1.4
          }
          if (!REDRAW) {
            rejs.known <- rejs #increased
            permleft <- permleft[!as.logical(rej.ind)]
          }
        }
        #
        if (VERBOSE) {
          cat(sprintf(sprintf("f.cur=%s,a.cur=%s; verify time: %s",SIGFMTf,sprintf("%%%d.%df",SIGFIGa+3,SIGFIGa+1),format(Sys.time()-START.VERIFY.TIME)),rejs/N.DRAWS,a.prev),
              file="",append=T,sep="\n")
          cat(sprintf(sprintf("f.cur=%s,a.cur=%s; verify time: %s",SIGFMTf,sprintf("%%%d.%df",SIGFIGa+3,SIGFIGa+1),format(Sys.time()-START.VERIFY.TIME)),rejs/N.DRAWS,a.prev),
              file=LOGFILE,append=T,sep="\n")
          cat(sprintf("It's %s",format(Sys.time(), "%X, %A, %d %b %Y")),sep='\n',append=TRUE,file="")
        }
      }
      alpha.tilde <- a.prev
      if (round(a.hi,digits=SIGFIGa)==round(a.lo,digits=SIGFIGa)) {
        #alpha.tilde=a.lo=a.hi, but f.lo<f.hi, so print both
      } else {
        f.lo <- f.hi <- f.prev #we hit f.prev close to ALPHA exactly
      }
      prev.alpha.tilde[which(ALPHAS==ALPHA)] <- alpha.tilde
      
      # Output results
      cat(sprintf(sprintf("%%5.3f,%%d,%%d,%s,%s,%s",SIGFMTa,SIGFMTf,SIGFMTf),
                  ALPHA, nx, ny, alpha.tilde, f.lo, f.hi),
          file=OUTFILE,sep="\n",append=TRUE)
      if (VERBOSE) {
        cat(sprintf(sprintf("alphatilde(alpha=%%5.3f,nx=%%d,ny=%%d)=%s;f.lo=%s;f.hi=%s",SIGFMTa,SIGFMTf,SIGFMTf),
                    ALPHA, nx, ny, alpha.tilde, f.lo, f.hi),
            file="",sep="\n",append=TRUE)
        cat(sprintf("Time elapsed: %s",format(Sys.time()-start.time)), 
            sep="\n",append=TRUE,file="")
        cat(sprintf(sprintf("alphatilde(alpha=%%5.3f,nx=%%d,ny=%%d)=%s;f.lo=%s;f.hi=%s",SIGFMTa,SIGFMTf,SIGFMTf),
                    ALPHA, nx, ny, alpha.tilde, f.lo, f.hi),
            file=LOGFILE,sep="\n",append=TRUE)
        cat(sprintf("Time elapsed: %s",format(Sys.time()-start.time)), 
            sep="\n",append=TRUE,file=LOGFILE)
      }
    } #ALPHAS loop
  } #nind loop
  return(TRUE) #all went smoothly
}

#EOF

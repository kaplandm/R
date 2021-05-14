# Feedback: KaplanDM@Missouri.edu
# Author (original): David M. Kaplan
# 18 Sept 2013
# Implementation of Kaplan (2015) confidence intervals based on fixed-smoothing asymptotics and Edgeworth expansion.
# Note: for empirical research, quantile_inf.R is recommended instead; it is a hybrid of methods including Goldman and Kaplan (2013), which allows for more types of object of interest and performs better when it is computable.
# References:
# 'Improved Quantile Inference via Fixed-smoothing Asymptotics and Edgeworth Expansion' by David M. Kaplan (2015), https://doi.org/10.1016/j.jeconom.2014.08.011
# See also https://doi.org/10.1016/j.jeconom.2016.09.015  and  https://doi.org/10.1111/ectj.12095


###################
#                 #
# Main function   #
#                 #
###################
# Function
# # Given sample data and quantile of interest, this function produces a two-sided confidence interval (CI) for a quantile of a population distribution, or for quantile treatment effects (given treatment and control samples).  See below (and examples in quantile_inf_examples.R) for details.
# Arguments
# # X: numeric vector of data, or for METHOD.TYPE 'qte', X$t is the numeric vector of treatment group data and X$c is the numeric vector of control group data.
# # p: quantile of interest, a numeric vector (or scalar) with values strictly between zero and one; e.g., the 0.5-quantile is the median.
# # ALPHA: the nominal coverage probability of the returned confidence interval (CI) is 1-ALPHA; e.g., ALPHA=0.05 produces a 95% CI.
# # ONESIDED: ONESIDED=0 produces a two-sided CI; ONESIDED=1 produces an upper one-sided CI; ONESIDED=-1 produces a lower one-sided CI.
# # METHOD.TYPE: 'single' produces a CI for a single quantile using one sample of data; 'qte' produces a CI for the quantile treatment effect (or, under fewer assumptions, simply the difference between the p-quantiles of two populations).  The two numeric vectors of data from the two (i.e., treatment and control) populations should be passed in X$t and X$c, respectively.
#
# Return value: a vector with the lower and upper endpoints of the CI, respectively
#
quantile.inf.K11 <- function(X,p,ALPHA=0.05,ONESIDED=0,METHOD.TYPE='single',PDFratio=NULL) {

  # Check arguments: errors
  if (missing(X)) stop("Argument X is missing.")
  if (missing(p)) stop("Argument p is missing.")
  quantile.inf.K11.validate.inputs(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED,METHOD.TYPE=METHOD.TYPE,PDFratio)

  # Check arguments: warnings
  if (!is.null(ALPHA) && ALPHA>0.2) warning('ALPHA is usually 0.01, 0.05, or 0.10; e.g., 0.05 gives a 95% confidence interval.')

  # Constants
  PDFMIN <- 1e-4 #for 2-sample
  UNIROOT.TOL <- 0.0001
  UNIROOT.INT <- c(0.000000001,1-0.000000001)
  KERNEL.TYPE <- 'gaussian'
  BW.TYPE <- 'SJ'

  # Prep data
  if (METHOD.TYPE=='qte') { X$t <- sort(X$t);  X$c <- sort(X$c) } else { X <- sort(X) } #X[i] is now order statistic i
  
  # Return if METHOD.TYPE 'single'
  if (METHOD.TYPE=='single') {
    ret <- quantile.inf.single.K11(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED)
    return(ret)
  } else if (METHOD.TYPE=='qte') {
    ret <- quantile.inf.qte.K11(X=X,p=p,ALPHA=ALPHA,ONESIDED=ONESIDED,PDFratio=PDFratio)
    return(ret)
  } #else: should already have been caught

} #end of function quantile.inf.K11()



###############################################################
# Pre-calculated simulated critical values from Kaplan (2011) #
###############################################################
sub.simcv <- function() {
simalphas <- c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,
             0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,
             0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,
             0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,
             0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,
             0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,
             0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,
             0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,
             0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,
             0.97,0.98,0.99)
simCVs.bi <- matrix(0,2,length(simalphas))
simCVs.bi[1,] <-
  c(10.9435,8.8905,7.8981,7.2232,6.7377,6.351,6.0229,5.7584,
    5.5418,5.3417,4.1981,3.6015,3.2189,2.9354,2.714,2.535,2.3849,2.2609,
    2.149,2.0503,1.9611,1.8804,1.8057,1.7389,1.6765,1.6177,1.5631,
    1.5123,1.4647,1.4196,1.3779,1.3379,1.2998,1.2629,1.2281,1.1953,
    1.1624,1.1308,1.1012,1.0723,1.0442,1.0169,0.99113,0.96603,0.94104,
    0.91779,0.8943,0.8723,0.85036,0.82907,0.80833,0.7879,0.76856,
    0.74943,0.73068,0.7125,0.69448,0.67636,0.65898,0.64175,0.62518,
    0.60846,0.59187,0.57614,0.56043,0.54512,0.53009,0.51467,0.49973,
    0.48494,0.47025,0.45579,0.44124,0.42726,0.4134,0.3995,0.38613,
    0.37253,0.35938,0.34629,0.33318,0.32009,0.30725,0.29457,0.28195,
    0.26964,0.25698,0.24456,0.23229,0.22018,0.20801,0.19615,0.18394,
    0.17234,0.16065,0.14903,0.13742,0.12555,0.11394,0.10242,0.090718,
    0.079071,0.067377,0.055898,0.04451,0.033114,0.021501,0.010147)
simCVs.bi[2,] <-
  c(5.9696,5.2541,4.8539,4.5537,4.3365,4.1601,4.0172,3.8932,3.7904,
    3.6945,3.1109,2.7768,2.5495,2.3795,2.2387,2.1231,2.0241,1.9387,
    1.8618,1.7919,1.7286,1.6696,1.6155,1.5656,1.5192,1.4756,1.4325,
    1.3925,1.3547,1.3189,1.2847,1.2522,1.2219,1.1917,1.1633,1.1353,
    1.1082,1.0819,1.0566,1.032,1.0086,0.98543,0.96225,0.94023,0.91818,
    0.89684,0.87671,0.85666,0.83658,0.81765,0.79923,0.78084,0.76276,
    0.74445,0.72695,0.70947,0.69251,0.67591,0.65949,0.64369,0.62773,
    0.6121,0.59668,0.58137,0.56623,0.55174,0.53661,0.52165,0.50672,
    0.49258,0.47876,0.46438,0.45022,0.43646,0.42261,0.40875,0.39521,
    0.38141,0.36817,0.35498,0.34185,0.3288,0.31576,0.30296,0.29003,
    0.27721,0.26461,0.25188,0.23945,0.2274,0.21519,0.20266,0.19044,
    0.17843,0.16643,0.15419,0.14254,0.13036,0.11814,0.10607,0.094184,
    0.08209,0.069939,0.057994,0.046069,0.034169,0.0224,0.010453)
simCVs.uni <- matrix(0,2,length(simalphas))
simCVs.uni[1,] <-
c(44.8943,31.1802,25.2111,21.6883,19.2587,17.4917,16.1646,14.8876,
 13.9785,13.1943,8.939,7.0955,5.9917,5.226,4.6791,4.2355,3.8938,3.6012,
 3.3534,3.1435,2.9647,2.8014,2.6569,2.5254,2.4059,2.2988,2.1986,2.1045,
 2.0205,1.9421,1.8691,1.8022,1.7373,1.6756,1.617,1.5641,1.5122,1.4644,
 1.4184,1.3736,1.3299,1.289,1.249,1.2112,1.1749,1.1398,1.107,1.0753,
 1.044,1.0146,0.98538,0.95693,0.92908,0.90212,0.87608,0.85127,0.82678,
 0.80332,0.78023,0.75758,0.7353,0.71397,0.69282,0.67244,0.65224,0.63245,
 0.61313,0.5938,0.57494,0.55698,0.5391,0.52137,0.50417,0.4875,0.47138,
 0.45488,0.43891,0.42298,0.40778,0.39269,0.37756,0.3624,0.34747,0.33319,
 0.3183,0.30394,0.2894,0.27537,0.26169,0.24795,0.23425,0.22045,0.20687,
 0.19345,0.18045,0.16714,0.15436,0.14132,0.12855,0.11578,0.1028,
 0.090287,0.077914,0.064961,0.052481,0.040066,0.027524,0.014898)
simCVs.uni[2,] <-
c(11.8392,9.6368,8.5465,7.7846,7.2393,6.8183,6.502,6.2144,5.9611,5.7361,
4.531,3.8996,3.4738,3.1766,2.942,2.7496,2.5918,2.4522,2.3315,2.2232,
2.1273,2.0418,1.9634,1.8916,1.8254,1.7632,1.7051,1.6499,1.5986,1.5497,
1.5024,1.4595,1.4175,1.3791,1.3411,1.3046,1.27,1.2356,1.2038,1.1726,
1.1425,1.1137,1.0863,1.059,1.0322,1.0064,0.98103,0.9573,0.93325,0.91049,
0.88797,0.86565,0.84382,0.82256,0.80216,0.78183,0.76214,0.74232,0.72384,
0.705,0.6862,0.66792,0.64977,0.63206,0.61481,0.598,0.58103,0.56456,
0.54839,0.53288,0.51723,0.50117,0.4856,0.47051,0.45574,0.44051,0.42602,
0.41153,0.39682,0.38232,0.36861,0.35478,0.34084,0.32665,0.313,0.29897,
0.28536,0.27191,0.25866,0.24548,0.23193,0.21882,0.20566,0.19261,0.17974,
0.16659,0.15379,0.14096,0.12814,0.11529,0.10286,0.090379,0.077713,
0.064897,0.052336,0.040043,0.027465,0.014858)
return(list(simalphas=simalphas,simCVs.bi=simCVs.bi,simCVs.uni=simCVs.uni))
}


#########################
#                       #
# Check input args      #
#                       #
#########################
# Check for valid input data
quantile.inf.K11.validate.inputs <- function(X,p,ALPHA,ONESIDED,METHOD.TYPE,PDFratio) {
  for (varname in c('X','p','ALPHA','ONESIDED','METHOD.TYPE')) {
    if (is.null(get(varname))) stop(sprintf("%s cannot be NULL.",varname))
  }
  if (!is.character(METHOD.TYPE)) { stop("Need is.character(METHOD.TYPE)=TRUE") }
  if (!(METHOD.TYPE %in% c('single','qte'))) stop("Argument METHOD.TYPE must be either 'single' or 'qte', which correspond respectively to one-sample inference on a single quantile and two-sample quantile treatment effect inference.")
  if (is.list(X)) {
    if (METHOD.TYPE!='qte') stop(sprintf("X should be a numeric vector for METHOD.TYPE '%s'",METHOD.TYPE))
    chk.num.vec.fn(X$t);  chk.num.vec.fn(X$c)
  } else {
    if (METHOD.TYPE=='qte') stop("X should be a list with X$t and X$c containing the two data samples for METHOD.TYPE 'qte'")
    chk.num.vec.fn(X)
  }
  if (missing(p)) stop("Need to supply quantile of interest, e.g. p=0.5 for median")
  chk.num.vec.fn(p)
  if (min(p)<0 || max(p)>1) stop('Argument p must be between 0 and 1; e.g., 0.5 for the median.')
  if (length(p)!=1) stop('p should be a scalar.')
  tmp <- length(X)
  if (METHOD.TYPE=='qte') tmp <- min(length(X$t),length(X$c))
  if (min(p)<1/tmp || max(p)>=1-1/tmp) stop(paste0('For this function, 1/n<=p<1-1/n is required.  ',
           'For values nearer to 0 or 1, ',
           'use methods based on extreme value theory.'))
  if (!(ONESIDED %in% c(-1,0,1))) stop(paste0('ONESIDED must be one of three values: 0 for two-sided inference;',
           ' -1 for one-sided with alternative hypothesis H1:xi<xi0',
           ' where confidence interval is of form (-Inf,b];',
           ' and 1 for one-sided with alternative hypothesis H1:xi>xi0',
           ' where confidence interval is of form [a,Inf).'))
  if (ONESIDED!=0 && !is.null(ALPHA) && ALPHA>=0.5) stop('ALPHA must be less than 0.5 for one-sided inference.')
  if (METHOD.TYPE=='qte' && !is.null(PDFratio) && !is.numeric(PDFratio)) { stop("PDFratio must be either NULL or numeric.") }
}


#########################
#                       #
# Check numeric vector  #
#                       #
#########################
chk.num.vec.fn <- function(X) {
    varname <- deparse(substitute(X))
    if (!is.vector(X) || is.list(X)) stop(sprintf("%1$s must be a vector (and not a list).  If %1$s is a one-dimensional array or matrix, you can use as.vector(%1$s) in order to cast %1$s as a vector.",varname))
    if (!is.numeric(X)) stop(sprintf("%1$s must be numeric; i.e., is.numeric(%1$s) must be TRUE.",varname))
}


#########################
#                       #
# list[a,b,...] <- fn() #
#                       #
#########################
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


#########################
#                       #
# Kaplan single CI      #
#                       #
#########################
quantile.inf.single.K11 <- function(X,p,ALPHA,ONESIDED) {
  #Load critical values simulated ahead of time for small values of m.
  simcvs <- sub.simcv()
  simalphas <- simcvs$simalphas
  simCVs.uni <- simcvs$simCVs.uni
  #Initialize depending on ONESIDED.
  if (ONESIDED==-1) {
    CI.lo <- -Inf;  zz <- qnorm(1-ALPHA,0,1)
  } else if (ONESIDED==1) {
    CI.hi <- Inf;  zz <- qnorm(ALPHA,0,1)
  } else if (ONESIDED==0) {
    zz <- qnorm(1-ALPHA/2,0,1)
  } else {
    stop('Argument ONESIDED can only be -1, 0, or 1.')
  }
  nx <- length(X);  sqrtnx <- sqrt(nx);  sqrtp1p <- sqrt(p*(1-p))
  rx <- floor(nx*p)+1  #order statistic index for the sample quantile
  Xnr <- X[rx]  #order statistic, and estimator of quantile
  mmax <- min(rx-1,nx-rx)      #smoothing parameter m can't be larger than this
  if (mmax<1) stop(sprintf("Can't have p as close to zero or one as p=%4.2f is for n=%d",p,nx))
  #
  localpwr <- 0.5 #1st-order power of alternative against which to max power
  #OPTIMIZE IS FASTER THAN OPTIM, AND MORE ROBUST THAN UNIROOT
  CC <- optimize(f=function(x) (pnorm(x+zz)-pnorm(x-zz)-(1-localpwr))^2,
                 interval=c(.1,6), maximum=FALSE, tol=0.00001)$minimum
  m <- nx^(2/3)*(CC*zz)^(1/3)*(3/4)^(1/3) *
       (((dnorm(qnorm(p)))^2)/(2*(qnorm(p))^2+1))^(1/3) *
       ((dnorm(zz-CC)-dnorm(zz+CC))/(dnorm(zz-CC)+dnorm(zz+CC)))^(1/3)
  m <- max(1,min(mmax,floor(m))) #make sure not <1 or >mmax
  Smn <- (nx/(2*m)) * (X[rx+m]-X[rx-m])  #SBG sparsity estimator
  ccv <- abs(zz + (zz^3)/(4*m))  #+ (zz^5+8*zz^3)/(96*m^2); #2nd-order works fine, no need for 3rd-order
  #use absolute value; eqn holds for 1- or 2-sided (ONESIDED=-1,0,1)
  if (m<=dim(simCVs.uni)[1]) { #For small m, ccv approx is inaccurate=>use simulated ccv
    if (ONESIDED!=0) ALPHA <- 2*ALPHA
    if (sum(simalphas==ALPHA)>0) {
      ccv <- simCVs.uni[m,simalphas==ALPHA]
    } else {
      #INTERPOLATE BTWN CONSECUTIVE ALPHAS
      tmp <- which(simalphas>=ALPHA)[1]
      if (ALPHA>max(simalphas)) {
        ccv <- simCVs.uni[m,end]
      } else if (ALPHA<min(simalphas)) {
        ccv <- simCVs.uni[m,1]
      } else {
        tmp2 <- (simalphas[tmp]-ALPHA)/(simalphas[tmp]-simalphas[tmp-1])
        ccv <- tmp2*simCVs.uni[m,tmp-1] + (1-tmp2)*simCVs.uni[m,tmp]
      }
    }
  }
  #Usual calculation of CI endpoints given critical value & std error.
  if (ONESIDED>=0) CI.lo <- Xnr - ccv*Smn*sqrtp1p/sqrtnx
  if (ONESIDED<=0) CI.hi <- Xnr + ccv*Smn*sqrtp1p/sqrtnx
  return(c(CI.lo,CI.hi))
}


#########################
#                       #
# Kaplan QTE CI         #
#                       #
#########################
#Pre-condition: X$t and X$c already sorted (ascending)
#SIM.CV=FALSE recommended since doesn't account for theta as well as ccv
quantile.inf.qte.K11 <- function(X,p,ALPHA,ONESIDED,PDFMIN=0.0001,KERNEL.TYPE='epanechnikov',BW.TYPE='SJ',PDFratio=NULL,SIM.CV=FALSE) {
  Y <- X$c;  X <- X$t
  #Load simulated critical values.
  if (SIM.CV) {
    simCVs.uni <- simCVs.bi <- NULL
    tryCatch({
      simcvs <- sub.simcv()
      simalphas<- simcvs$simalphas
      simCVs.uni <- simcvs$simCVs.uni
      simCVs.bi  <- simcvs$simCVs.bi
    }, error=function(ERR) {}) #warning("Could not load simulated fixed-m critical values inside function quantile.inf.qte.K(); using approximated CVs instead.")
  }
  nx <- length(X);  ny <- length(Y)
  rx <- floor(nx*p)+1; ry <- floor(ny*p)+1
  Xnr <- X[rx];  Ynr <- Y[ry] #sample quantiles for X, Y
  sqrtp1p <- sqrt(p*(1-p));  sqrtnx <- sqrt(nx)  #sqrtny=sqrt(ny);
  zz <- qnorm(1-ALPHA/2,0,1)  #standard normal critical value, 2-sided
  if (ONESIDED!=0) zz <- qnorm(1-ALPHA,0,1)
  localpwr <- 0.5 #1st-order power of alternative against which to max power
  CC <- optimize(f=function(x) (pnorm(x+zz)-pnorm(x-zz)-(1-localpwr))^2,
                 interval=c(.1,6), maximum=FALSE, tol=0.00001)$minimum
  #
  #Pilot estimate of theta to plug into m
  xeval <- X[rx];  yeval<-Y[ry] #where to evaluate PDFs when calculating theta
  if (is.null(PDFratio)) {
    fxest<-max(PDFMIN,
               density(X,from=xeval,to=xeval,n=1,kernel=KERNEL.TYPE,bw=BW.TYPE)$y)
    fyest<-max(PDFMIN,
               density(Y,from=yeval,to=yeval,n=1,kernel=KERNEL.TYPE,bw=BW.TYPE)$y)
    delta <- fxest/fyest
  } else {
    delta <- PDFratio
  }
  #     theta.est <- (fxest^(-2)+(nx/ny)*fyest^(-2))^(-2) * (fxest^(-4)+(nx^2/ny^2)*fyest^(-4))
  theta.est <- (1+(nx/ny)*delta^2)^(-2) * (1+(nx^2/ny^2)*delta^4)
  #
  mmax <- min(min(rx-1,nx-rx),min(ry-1,ny-ry))
  m <- max(nx,ny)^(2/3)*(3/4)^(1/3) *
    ((1-theta.est)+theta.est*(CC*zz)*(dnorm(zz-CC)-dnorm(zz+CC)) /
       (dnorm(zz-CC)+dnorm(zz+CC)))^(1/3) *
    (((dnorm(qnorm(p)))^2)/(2*(qnorm(p))^2+1))^(1/3)
  m <- max(1,min(mmax,floor(m))) #make sure not <1 or >mmax
  Sx <- (nx/(2*m)*(X[rx+m]-X[rx-m]))
  Sy <- (ny/(2*m)*(Y[ry+m]-Y[ry-m]))
  Smn <- sqrt(Sx^2+(nx/ny)*Sy^2)
  if (is.null(PDFratio)) {
    delta <- Sy/Sx #rem: S is sparsity, not density
  } else {
    delta <- PDFratio
  }
  theta.est <- (1+(nx/ny)*delta^2)^(-2) * (1+(nx^2/ny^2)*delta^4)
  #         ((nx/(2*m)*(X[rx+m]-X[rx-m]))^4+(nx^2/ny^2)*(ny/(2*m)*(Y[ry+m]-Y[ry-m]))^4) /  (Smn^4)  #shrink a little toward 1 by adding (.1+...)/(.1+...)?
  c1 <- (1/4)*(theta.est*zz^3-zz*(1-theta.est))
  ccv1 <- zz + c1/m
  c2 <- (1/32)*(zz*(17*theta.est^2-30*theta.est+13) +zz^3*(16*theta.est^2-20*theta.est+(20/3)) +zz^5*(7*theta.est^2-10*theta.est+(10/3)))
  ccv2 <- zz + c1/m + c2/m^2
  ccv <- ccv2 #for 2-sample, sometimes ccv1>ccv2, like for theta=3/4; unlike 1s, can be >simCV, too
  #Using simulated CVs is somewhat ad hoc when nx!=ny
  if (SIM.CV && m<=dim(simCVs.uni)[1]) {
    if (ONESIDED!=0) ALPHA <- 2*ALPHA
    tmpw1 <- 2*theta.est-1; tmpw2 <- 2-2*theta.est #weight btwn uni and "bi" sim'd ccvs
    #rem: theta is between 1/2 and 1, hence not just tmpw1=theta
    if (sum(simalphas==ALPHA)>0)
      ccv <- tmpw1*simCVs.uni[m,simalphas==ALPHA] +
      tmpw2*simCVs.bi[m,simalphas==ALPHA]  #now vector, nreplic x 1
    else {
      #INTERPOLATE BTWN CONSECUTIVE ALPHAS
      tmp <- which(simalphas>=ALPHA)[1]
      if (ALPHA>max(simalphas))
        ccv <- tmpw1*simCVs.uni[m,end] +
        tmpw2*simCVs.bi[m,end]
      else if (ALPHA<min(simalphas))
        ccv <- tmpw1*simCVs.uni[m,1] +
        tmpw2*simCVs.bi[m,1]
      else {
        tmp2<-(simalphas[tmp]-ALPHA)/(simalphas[tmp]-simalphas[tmp-1])
        ccv <- tmpw1*(tmp2*simCVs.uni[m,tmp-1] + (1-tmp2)*simCVs.uni[m,tmp]) +
          tmpw2*(tmp2*simCVs.bi[m,tmp-1] + (1-tmp2)*simCVs.bi[m,tmp])
      }
    }
  }
  #Usual CI endpoint calculation based on critical value, std error.
  CI.lo <- -Inf; CI.hi <- Inf
  if (ONESIDED>=0) CI.lo <- Xnr-Ynr - ccv*Smn*sqrtp1p/sqrtnx
  if (ONESIDED<=0) CI.hi <- Xnr-Ynr + ccv*Smn*sqrtp1p/sqrtnx
  return(c(CI.lo,CI.hi))
}


## EOF
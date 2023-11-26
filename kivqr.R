# k-class IVQR estimation
# by Dave Kaplan & Xin Liu
# Questions? Comments? kaplandm@missouri.edu

# Load code
prob.loaded <- exists("gmmq")
success <-
  tryCatch({source('https://raw.githubusercontent.com/kaplandm/R/main/gmmq.R');TRUE},
           error=function(w) FALSE)
if (!success) {
  success <- tryCatch({source("gmmq.R");TRUE}, error=function(w) FALSE)
  if (!success) {
    if (prob.loaded) {
      warning("Couldn't load gmmq.R, but it seems like you already did.")
    } else {
      stop("Failed to source() gmmq.R from web or local file.  You may download and source() it yourself, or at least make sure it's in your getwd().  Currently available at https://github.com/kaplandm/R")
    }
  }
}

# tau: quantile index, 0<tau<1 (0.5=median, etc.) [required]
# Y: vector of scalar outcome for n observations [required]
# X.endog: matrix of endogenous regressors (n rows) [required]
# X.exog: matrix of exogenous regressors (n rows) [optional]
# Z.excl: matrix of excluded instruments (n rows) [required]
# h: smoothing bandwidth (can set h=0 for minimal smoothing; if unspecified, uses Kaplan & Sun (2017) plug-in)
# k: value of k for k-class; default (NULL) solves for "LIML" k that reliably reduces median bias in simulations
# b.init: initial coefficient vector estimate (default 0,...,0 or if available rq(Y~X.endog+X.exog)$coefficients[2:dB,1] )
# ! an intercept will be included automatically;
#   do *not* include one in X.exog or Z.excl
# ! careful of the order of coefficients: first endogenous,
#   then exogenous, then intercept
# LIMLadj: 
kivqr <- function(tau, Y, X.endog=NULL, X.exog=NULL, Z.excl=NULL, h=0, k=NULL, b.init=NULL, LIMLadj=TRUE) {
  # Parse/set additional arguments
  if (is.null(Y) || is.null(X.endog) || is.null(Z.excl)) stop("Must specify Y (vector of outcome observations), X.endog (matrix of endogenous regressor observations, n rows), and Z.excl (matrix of excluded instruments for X.endog, also n rows).")
  # make any vectors into n-by-1 matrix (and any data.frame into matrix type)
  Y <- as.matrix(Y); X.endog <- as.matrix(X.endog); Z.excl <- as.matrix(Z.excl)
  if (!is.null(X.exog)) X.exog <- as.matrix(X.exog)
  if (is.null(k)) { #solve for LIML k if no k specified by user
    ENDO <- cbind(Y,X.endog)
    if (is.null(X.exog)) {
      detfn <- function(k) det(t(ENDO)%*%lm(ENDO~1)$residuals -
                                 k*t(ENDO)%*%lm(ENDO~Z.excl)$residuals)
    } else {
      detfn <- function(k) det(t(ENDO)%*%lm(ENDO~X.exog)$residuals -
                                 k*t(ENDO)%*%lm(ENDO~X.exog+Z.excl)$residuals)
    }
    FTOL <- 0.0001;  LIMLITER <- 20
    lo <- 0.9;  detlo <- detfn(lo)
    for (r in -LIMLITER:-1) {
      kLIML <- tryCatch(uniroot(f=detfn, tol=FTOL, maxiter=LIMLITER,
                                interval=c(lo,1+2^r), extendInt='no',
                                f.lower=detlo)$root,
                        error=function(err)NA)
      if (!is.na(kLIML)) break
    }
    kLIML2 <- tryCatch(optimize(f=function(x)(detfn(x))^2, tol=FTOL,
                                interval=c(lo, 1+3*(max(KS[KS>0])-1)) )$minimum,
                       error=function(err)NA)
    if (!is.na(kLIML2) && (is.na(kLIML) || kLIML2<kLIML)) kLIML <- kLIML2
    # 
    # sanity check: set to 1 if clearly too low or too high
    tmp <- ncol(cbind(Z.excl))
    if (is.na(kLIML)) {
      warning('LIML k NA; using k=1 instead.')
      kLIML <- 1
    } else if (LIMLadj && kLIML<0.99) {
      warning('LIML k NA or too low; using k=1 instead.')
      kLIML <- 1
    } else if (LIMLadj && kLIML>1+3*(max(c(n/(n-tmp+2), 1+tmp/(n-tmp)))-1)) { #max(bias-adjusted, JIVE)
      warning('LIML k too big; using k=1 instead.')
      kLIML <- 1
    }
    k <- kLIML
  }
  if (is.null(b.init)) {
    if (require(quantreg)) {
      if (is.null(X.exog)) {
        tmp <- coef(rq(Y~X.endog))
      } else if (is.null(X.endog)) {
        tmp <- coef(rq(Y~X.exog))
      } else {
        tmp <- coef(rq(Y~X.endog+X.exog))
      }
      b.init <- tmp[c(2:length(tmp),1)]
    } else {
      b.init <- 0
    }
  }
  # Construct k-class instrument matrix
  if (is.null(X.exog)) {
    Ztilde <- cbind( (1-k)*cbind(X.endog) + 
                       k*cbind(lm(X.endog~Z.excl)$fitted.values) , 1 )
  } else {
    Ztilde <-
      cbind((1-k)*cbind(X.endog, X.exog) + 
            k*cbind(lm(X.endog~X.exog+Z.excl)$fitted.values, X.exog) , 1 )
  }
  # Call gmmq()
  ret <- gmmq(tau=tau, Y=cbind(Y, X.endog), X=cbind(X.exog, rep(1,length(Y))),
              Z.excl=Ztilde, h=h, b.init=b.init)
  ret$b <- c(ret$b)
  if (is.null(X.exog)) {
    names(ret$b) <- c(sprintf('endog.%d',1:ncol(cbind(X.endog))),'(Intercept)')
  } else {
    names(ret$b) <- c(sprintf('endog.%d',1:ncol(cbind(X.endog))),
                      sprintf('exog.%d',1:ncol(cbind(X.exog))),'(Intercept)')
  }
  ret$k <- k
  return(ret)
}

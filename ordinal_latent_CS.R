# Computes inner CS for latent quantile differences given ordinal data
# Paper: "Comparing latent inequality with ordinal data" (Kaplan, Zhao, Zhuo)
# Original author: Wei Zhao
# Updates: David M. Kaplan   kaplandm.github.io


# Arguments:
# X, Y: numeric vectors with ordinal data; for example, 1 1 3 1 2 2 3
# alpha: the confidence level is 1-alpha; default alpha=0.1 yields 90% confidence level
# N: number of draws for simulating critical values
# ordinal.values: manually specify values at which to compare CDFs,
#                 otherwise determined automatically
# weightX, weightY: weights to apply to X and Y (when estimating CDFs)
# FXhat, FYhat, nX, nY: instead of specifying X and Y (raw data), can provide
#                       estimated CDFs and samples sizes.
# 
# Return value: list with "tilde alpha" and "tilde beta" (see paper, Method 1),
#    ordinal values used for analysis (just argument ordinal.values if provided),
#    confidence limits for X and Y,
#    and the inner confidence set itself,
#        which is the union of all the [lower,upper] ranges
ordinal.latent.CS <- function(X, Y, alpha=0.1, N=1e4, ordinal.values, 
                              weightX, weightY, FXhat, FYhat, nX, nY){ 
  # Set/save seeds, for replicability
  oldseed <- NULL
  if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
    oldseed <- get(".Random.seed",.GlobalEnv)
  }
  on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
  set.seed(112358)

  # Find ordinal points for analysis
  if (missing(ordinal.values)) {
    VV <- NULL
  } else if (is.null(ordinal.values) || !is.numeric(ordinal.values) || any(is.na(ordinal.values))) {
    warning("Argument ordinal.values must be numeric and not contain NA; ignoring and determining automatically.")
    VV <- NULL
  } else {
    VV <- ordinal.values
  }
  if (is.null(VV)) {
    uX <- sort(unique(X))
    uY <- sort(unique(Y))
    VV <- sort(intersect(unique(Y), unique(X)))
    tmp <- which(uX==min(VV))
    if (tmp>1) VV <- c(uX[tmp-1], VV)
    tmp <- which(uY==min(VV))
    if (tmp>1) VV <- c(uY[tmp-1], VV)
    tmp <- which(uX==max(VV))
    if (tmp<length(uX)) VV <- c(VV, uX[tmp+1])
    tmp <- which(uY==max(VV))
    if (tmp<length(uY)) VV <- c(VV, uY[tmp+1])
  }
  
  # J=number of categories, labeled 1,2,...,J
  J <- length(VV)
  
  # Estimated CDFs, evaluated at the points for analysis
  if (!missing(FXhat)) {
    if (missing(FYhat) || missing(nX) || missing(nY)) {
      stop("Arguments FXhat, FYhat, nX, nY must all be provided to be used (or else none provided, to use X and Y arguments).")
    } else if (!is.numeric(FXhat) || !is.numeric(FYhat) || !is.numeric(nX) || !is.numeric(nY)) {
      stop("Arguments FXhat, FYhat, nX, nY must be all numeric or all missing.")
    } else if (length(FXhat)!=length(FYhat) || length(FXhat)!=J) {
      stop("Arguments FXhat and FYhat must be numeric vectors with same length as ordinal.values.")
    } else if (length(nX)!=1 || length(nY)!=1) {
      stop("Arguments nX and nY must be scalars.")
    }
    Fhat_X <- FXhat
    Fhat_Y <- FYhat
  } else {
    nX <- length(X) # sample size
    nY<- length(Y) # sample size
    if (missing(weightX)) {
      Fhat_X <- ecdf(X)(VV) # estimated CDF
    } else {
      weightsumX <- sum(weightX)
      Fhat_X <- rep(NA, J)
      Fhat_X[J] <- 1
      for (j in 1:(J-1)) {
        Fhat_X[j] <- sum(weightX[X<=j]) / weightsumX
      }
    }
    if (missing(weightY)) {
      Fhat_Y <- ecdf(Y)(VV) # estimated CDF
    } else {
      weightsumY <- sum(weightY)
      Fhat_Y <- rep(NA, J)
      Fhat_Y[J] <- 1
      for (j in 1:(J-1)) {
        Fhat_Y[j] <- sum(weightY[Y<=j]) / weightsumY
      }
    }
  }
 
  # To be a little conservative: replace 0 with 1/n, replace 1 with 1-1/n
  Fhat_X[Fhat_X==0] <- 1/nX
  Fhat_Y[Fhat_Y==0] <- 1/nY
  Fhat_X[Fhat_X==1] <- (1-1/nX)
  Fhat_Y[Fhat_Y==1] <- (1-1/nY)
  
  # Estimate the covariance matrices
  SigmaXhat <- SigmaYhat <- matrix(data=NA, nrow=J-1, ncol=J-1)
  for (i in 1:(J-1)) {
    for (k in i:(J-1)) {
      SigmaXhat[i,k] <- SigmaXhat[k,i] <- Fhat_X[i]*(1-Fhat_X[k])
      SigmaYhat[i,k] <- SigmaYhat[k,i] <- Fhat_Y[i]*(1-Fhat_Y[k])
    }
  }

  # Estimate Sigma_tX and Sigma_tY from paper (t-statistic covariance matrices)
  var_covX <- diag(1/sqrt(diag(SigmaXhat)) )%*% SigmaXhat %*% t(diag(1/sqrt(diag(SigmaXhat))))
  var_covY <- diag(1/sqrt(diag(SigmaYhat)) )%*% SigmaYhat %*% t(diag(1/sqrt(diag(SigmaYhat))))
  
  # Critical value simulation
  # Simulate N normal vectors (for X t-statistics)
  chX <- t(chol(var_covX))
  mvnX <- t(chX %*% matrix(rnorm(N*(J-1)), nrow=J-1, ncol=N))
  Tmin <- apply(X=mvnX, MARGIN=1, FUN=min) # take minimum
  rm(mvnX)  # free up memory from the big normal matrix
  # Compute \tilde\alpha from paper (1-\tilde\alpha = pointwise coverage prob)
  talpha <- quantile(pnorm(Tmin), probs=1-sqrt(1-alpha), type=6)
  ConfX <- rep(NA, J-1) # confidence limits from Method 1, eqn (15)
  for (j in 1:(J-1)) ConfX[j]<- Fhat_X[j] + qnorm(1-talpha)*sqrt(SigmaXhat[j,j]/nX)

  # Simulate N normal vectors (for Y t-statistics)
  chY <- t(chol(var_covY))
  mvnY <- t(chY %*% matrix(rnorm(N*(J-1)), nrow=J-1, ncol=N))
  Tmax <- apply(X=mvnY, MARGIN=1, FUN=max) # take maximum
  rm(mvnY)  # free up memory from the big normal matrix
  # Compute \tilde\beta from paper (1-\tilde\beta = pointwise coverage prob)
  tbeta <- 1 - quantile(pnorm(Tmax), probs=sqrt(1-alpha), type=6)
  ConfY<- rep(NA,J-1) # confidence limits from Method 1, eqn (15)
  for (j in 1:(J-1)) ConfY[j]<- Fhat_Y[j] - qnorm(1-tbeta)*sqrt(SigmaYhat[j,j]/nY)
  
  # Aggregate quantile index intervals into overall inner confidence set
  innerCSranges <- NULL
  for (j in 1:(J-1)) {
    if (ConfX[j]<ConfY[j]) {
      innerCSranges <- rbind(innerCSranges, c(lower=ConfX[j], upper=ConfY[j]))
    }
  }
  return(list(talpha=talpha, tbeta=tbeta, ordinal.values=VV, 
              ConfX=ConfX, ConfY=ConfY, innerCSranges=innerCSranges))
}
# EOF
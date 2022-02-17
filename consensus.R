# Code implementing (shifted) CRRA utility version of
# "Inference on Consensus Ranking of Distributions"
# by David M. Kaplan (me)
# Questions? kaplandm@missouri.edu

# Future extensions to do: 1) stepdown, 2) pretest

# Input: vectors Y0 and Y1 (and optionally other arguments).
#   The utility function has
#       risk aversion parameter theta varying across THETA.GRID
#       and shift parameter s varying across SHIFT.GRID,
#       applying the CRRA utility function to y-s.
# Output: confidence band for expected utility difference,
#   and inner and outer confidence sets,
#   in this case sets of (theta,s).
#   The return value is a data.frame whose first two columns show
#   s and theta values, followed by the lower and upper functions of the
#   uniform confidence band (Lband and Uband), followed by the inner and
#   outer confidence sets (TRUE if value is included, FALSE if not).
#   The 2-sided band has level (1-ALPHA)*100%, 
#   which is also P(innerCS \subseteq TrueSet \subseteq outerCS).
#   If SYMMETRIC=FALSE, then everything is "equal-tailed" and
#   the 1-sided bands and CSs each have confidence level (1-ALPHA/2)*100%
#   (For example, use ALPHA=0.1 and SYMMETRIC=FALSE to get a one-sided 95% inner CS.)
# Note: if THETA.GRID is not sufficiently rich, then
#   the outer CS may not truly be an outer CS, and
#   the inner CS may be conservative.
consensus <- function(Y0, Y1, ALPHA=0.1, SHIFT.GRID=0, SYMMETRIC=TRUE, THETA.GRID=seq(from=0, to=3, by=0.2), BREP=199L, na.rm=FALSE) {

  ### VALIDATE ARGUMENTS
  
  if (!is.numeric(ALPHA) || length(ALPHA)!=1 || is.na(ALPHA) || ALPHA<=0 || ALPHA>=1) {
    stop("ALPHA must be a number between 0 and 1 (exclusive)")
  }
  if (!is.numeric(SHIFT.GRID)) stop("SHIFT.GRID must be a numeric vector.")
  SHIFT.GRID <- sort(c(SHIFT.GRID))
  if (!is.vector(SHIFT.GRID)) stop("SHIFT.GRID must be a numeric vector.")
  if ((max(THETA.GRID)>=1) && (min(Y0,Y1)-max(SHIFT.GRID)<=0)) stop("Must have strictly positive min(Y0,Y1)-max(SHIFT.GRID)")
  if (!is.numeric(BREP) || length(BREP)!=1 || is.na(BREP) || BREP<=0 || abs(BREP-round(BREP))>0.001) {
    stop("BREP must be a positive integer (number of bootstrap replications)")
  } else { BREP <- as.integer(round(BREP)) }
  if (!is.numeric(THETA.GRID) || length(THETA.GRID)<1 || any(is.na(THETA.GRID)) || any(THETA.GRID<0)) {
    stop("THETA.GRID must be a vector of non-negative numbers")
  }
  if (anyDuplicated(THETA.GRID)>0) {
    warning("Removing duplicate values from THETA.GRID")
    THETA.GRID <- unique(THETA.GRID)
  }
  if (is.unsorted(THETA.GRID)) {
    warning("sorting THETA.GRID (was unsorted)")
    THETA.GRID <- sort(THETA.GRID)
  }
  if (missing(Y0) || missing(Y1)) stop("Must supply Y0 and Y1")
  if (!is.numeric(Y0) || length(Y0)<1) {
    stop("Y0 must be a vector of numbers")
  }
  if (!is.numeric(Y1) || length(Y1)<1) {
    stop("Y1 must be a vector of numbers")
  }
  if (is.na(na.rm)) {
    stop("na.rm must be TRUE or FALSE")
  }
  # Deal with infinite, missing (NA), and NaN in Y0 and Y1
  if (any(is.infinite(Y0)) || any(is.infinite(Y1))) {
    stop("Y0 and Y1 cannot contain infinite values")
  }
  if (any(is.na(Y0)) || any(is.na(Y1)) || any(is.nan(Y0)) || any(is.nan(Y1))) {
    if (na.rm) {
      Y0 <- Y0[!is.na(Y0) & !is.nan(Y0)]
      Y1 <- Y1[!is.na(Y1) & !is.nan(Y1)]
    } else {
      stop("missing values (NA) and NaN's not allowed if 'na.rm' is FALSE")
    }
  }
  
  ### ACTUAL METHOD

  # Sample sizes
  n0 <- length(Y0);  n1 <- length(Y1);  LAMBDA <- sqrt(n1/n0)

  # Define utility function over grid using CRRA Uvecfn
  if (any(THETA.GRID==1)) {
    Uvfn <- function(y) {
      ret <- NULL
      for (s in SHIFT.GRID) {
        ret <- c(ret, Uvecfn(scalar.y=y-s, GRID=THETA.GRID, T1IND=which(THETA.GRID==1)))
      }
      return(ret)
    }
  } else {
    Uvfn <- function(y) {
      ret <- NULL
      for (s in SHIFT.GRID) {
        ret <- c(ret, Uvecfn(scalar.y=y-s, GRID=THETA.GRID))
      }
      return(ret)
    }
  }

  # Compute sample values
  f0 <- apply(X=as.array(Y0), MARGIN=1, FUN=Uvfn)
  f1 <- apply(X=as.array(Y1), MARGIN=1, FUN=Uvfn)
  Pn0 <- rowMeans(x=f0);  Pn1 <- rowMeans(x=f1)
  Pn <- Pn1 - Pn0 # sample expected utility difference
  Vn0 <- apply(X=f0, MARGIN=1, FUN=var)
  Vn1 <- apply(X=f1, MARGIN=1, FUN=var)
  SDn <- sqrt(Vn1+LAMBDA^2*Vn0) # sample std dev for EU diff
  
  # Save/set RNG seed
  oldseed <- NULL
  if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
    oldseed <- get(".Random.seed",.GlobalEnv)
  }
  on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
  set.seed(112358)

  # Exchangeable bootstrap
  Tnstars <- matrix(data=NA, nrow=BREP, ncol=length(Pn))
  for (b in 1:BREP) {
    # Bayesian bootstrap weights (Dirichlet)
    w0 <- rgamma(n=n0, shape=1, scale=1);  w0 <- w0/sum(w0)
    w1 <- rgamma(n=n1, shape=1, scale=1);  w1 <- w1/sum(w1)
    w0mean <- function(Y) sum(w0*Y)
    w1mean <- function(Y) sum(w1*Y)
    Pn0star <- apply(X=f0, MARGIN=1, FUN=w0mean)
    Pn1star <- apply(X=f1, MARGIN=1, FUN=w1mean)
    Pnstar <- Pn1star - Pn0star
    Tnstars[b,] <- sqrt(n1)*(Pnstar - Pn)
  }
  ZIQR <- qnorm(0.75)-qnorm(0.25) #N(0,1) interquartile range
  # if (mean(is.na(Tnstars)+is.nan(Tnstars))>ALPHA) stop("Too many NaN or NA in bootstrap t statistics.")
  SDstar <- apply(X=Tnstars, MARGIN=2, FUN=IQR, type=6, na.rm=TRUE) / ZIQR
  problemind <- max(0,c(which(SDstar==0),which(is.na(SDstar)),which(is.nan(SDstar))))
  if (problemind>0) {
    problemshiftind <- ceiling(problemind/length(THETA.GRID))
    problemthetaind <- problemind - length(THETA.GRID)*(problemshiftind-1)
    stop(sprintf("Insufficient bootstrap variation with SHIFT=%g and THETA=%g",
                 SHIFT.GRID[problemshiftind],
                 THETA.GRID[problemthetaind]))
  }
  if (SYMMETRIC) {
    Tmaxstars <- apply(X=Tnstars, MARGIN=1, FUN=function(row)max(abs(row)/SDstar))
    if (mean(is.na(Tmaxstars)+is.nan(Tmaxstars))>ALPHA) stop("Too many NaN or NA in bootstrap max t statistics.")
    cv1 <- cv2 <- quantile(x=Tmaxstars, probs=1-ALPHA, type=6, na.rm=TRUE)
  } else {
    Tmaxstars <- apply(X=Tnstars, MARGIN=1, FUN=function(row)max(row/SDstar))
    Tminstars <- apply(X=Tnstars, MARGIN=1, FUN=function(row)min(row/SDstar))
    if (mean(is.na(Tmaxstars)+is.nan(Tmaxstars))>ALPHA/2) stop("Too many NaN or NA in bootstrap max t statistics.")
    if (mean(is.na(Tminstars)+is.nan(Tminstars))>ALPHA/2) stop("Too many NaN or NA in bootstrap min t statistics.")
    cv1 <-  quantile(x=Tmaxstars, probs=1-ALPHA/2, type=6, na.rm=TRUE)
    cv2 <- -quantile(x=Tminstars, probs=  ALPHA/2, type=6, na.rm=TRUE)
  }
  
  # Uniform confidence band for expected utility difference
  Lband <- Pn - cv1*pmax(SDn,SDstar)/sqrt(n1)
  Uband <- Pn + cv2*pmax(SDn,SDstar)/sqrt(n1)
  
  # Confidence sets for utility fns w/ E[f(Y0)]>=E[f(Y1)]
  innerCSind <- Lband>=0 #as.integer(Lband>=0)
  outerCSind <- Uband>=0 #as.integer(Uband>=0)
  
  return(data.frame(SHIFT=rep(x=SHIFT.GRID, each=length(THETA.GRID)),
                    THETA.GRID=rep(x=THETA.GRID, times=length(length(SHIFT.GRID))),
                    Lband=Lband, Uband=Uband, 
                    innerCSind=innerCSind, outerCSind=outerCSind))
}

# CRRA utility function
# Input: scalar y; grid of theta values; index of theta=1
# Output: vector of u(y;theta) over the GRID of theta
Uvecfn <- function(scalar.y, GRID, T1IND=NA) {
  scalar.y <- scalar.y
  tmp <- (scalar.y^(1-GRID)-1)/(1-GRID)
  if (!is.na(T1IND)) tmp[T1IND] <- log(scalar.y)
  return(tmp)
}

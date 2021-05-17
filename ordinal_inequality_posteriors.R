# Author: David M. Kaplan  https://kaplandm.github.io
# Code to compute Bayesian posteriors in
#   "Comparing Latent Inequality with Ordinal Data"
#   by David M. Kaplan and Longhao Zhuo
# Computes posteriors for certain relationships
# between two population ordinal distributions,
# assuming iid sampling of both (independent) samples.
# To interpret these ordinal relationships in terms
# of latent variables, see the paper.
# The prior is an uninformative uniform Dirichlet(1,...,1)
# prior over the parameter space (category probabilities)
# for each distributions separately.  This can be adjusted
# optionally; see adj.prior below.

# ARGUMENTS:
# *either* Xcounts and Ycounts *or* pXdraws and pYdraws must be specified.
#   With Xcounts and Ycounts, iid sampling is assumed; 
#   with complex sampling designs, draws from the posterior can be 
#   generated separately and passed in the B-by-J matrices pXdraws and
#   pYdraws, where B is the number of posterior draws 
#   and J is the number of categories, and each row of each matrix 
#   contains the probability for each category (i.e., the PMF).
# Xcounts: vector of counts (integers) for each category; e.g., 
#   for SRHS with 5 categories 
#   poor, fair, good, very good, excellent,
#   Xcounts=c(2,3,9,8,6) indicates having observed a
#   data sample with 
#   2 "poor" observations, 3 fair, 9 good, 8 very good, 
#   6 excellent.
# Ycounts: same as Xcounts but for 2nd data sample.
#   The goal is to learn about the ordinal population 
#   distributions from which Xcounts and Ycounts were sampled.
# adj.prior: FALSE (default) to use a uniform Dirichlet(1,...,1)
#   prior over the parameter vector for the category
#   probabilities for X (i.e., the PMF values),
#   and similarly for Y. TRUE to adjust the prior
#   to have probability 1/2 on each "relationship"
#   examined, like P(X SD1 Y)=1/2 (so there is a separate
#   posterior for each relationship, whereas there is
#   a single posterior when adj.prior=FALSE).
# DRAWS: how many draws (default 1000) from the posterior to take. Ignored if pXdraws and pYdraws are used.
# 
# Return VALUE: a vector with the following named components:
# X.SD1.Y: posterior probability that X SD1 Y
# Y.SD1.X: posterior probability that Y SD1 X
# X.SC.Y: posterior probability that the CDF of X
#   crosses the CDF of Y from below once (single crossing).
# Y.SC.X: posterior probability that the CDF of Y
#   crosses the CDF of X from below once.
# multiple.crossings: posterior probability that the
#   ordinal CDFs cross more than once (not SD1 or SC).
# X.MPS.Y: posterior probability that X is a median-
#   preserving spread of Y; see Allison and Foster (2004).
# Y.MPS.X: posterior probability that Y is a median-
#   preserving spread of X; see Allison and Foster (2004).
# 
# EXAMPLES
# ordinal.inequality.posteriors(Xcounts=c(60,7,7,7,7), Ycounts=c(15,0,0,0,15))
# ordinal.inequality.posteriors(Xcounts=c(1,3,4,19,5), Ycounts=c(5,7,8,9,3))
ordinal.inequality.posteriors <- function(Xcounts=NULL, Ycounts=NULL, pXdraws=NULL, pYdraws=NULL, adj.prior=FALSE, DRAWS=1000L) {
  if (is.na(adj.prior) || !is.logical(adj.prior)) stop("Argument adj.prior must be either TRUE or FALSE.")
  if (!is.numeric(DRAWS)) stop("Argument DRAWS must be an integer, like 100L.")
  if (DRAWS<1) stop("Argument DRAWS must be a positive integer.")
  if (!is.integer(DRAWS)) DRAWS <- as.integer(DRAWS)

  if (is.null(pXdraws) && is.null(pYdraws)) {
    if (is.null(Xcounts) || is.null(Ycounts)) stop("Arguments Xcounts and Ycounts must both be provided if pXdraws and pYdraws are not.  Xcounts and Ycounts should be numerical vectors; see code comments for details.")
    if (!is.numeric(Xcounts) || !is.numeric(Ycounts)) stop("Arguments Xcounts and Ycounts both must be vectors of integers, so is.numeric(Xcounts) and is.numeric(Ycounts) should both be true.")
  } else {
    if (!is.null(Xcounts) || !is.null(Ycounts)) {
      warning("Ignoring Xcounts and Ycounts since pXdraws and pYdraws are provided.")
      Xcounts <- Ycounts <- NULL
    }
    if (!is.numeric(pXdraws) || !is.numeric(pYdraws)) stop("Arguments pXdraws and pYdraws both must be integer matrices, so is.numeric(pXdraws) and is.numeric(pYdraws) should both be true.")
    if (!is.matrix(pXdraws) || !is.matrix(pYdraws)) stop("When provided, must have is.matrix(pXdraws)=TRUE and is.matrix(pYdraws)=TRUE.")
    if (any(dim(pXdraws)!=dim(pYdraws))) stop("Must have dim(pXdraws)=dim(pYdraws).")
    DRAWS <- dim(pXdraws)[1]
    J <- dim(pXdraws)[2]
  }

  relationships <- c('X.SD1.Y','Y.SD1.X', 'X.SC.Y','Y.SC.X', 'multiple.crossings','X.MPS.Y','Y.MPS.X')
  ret <- rep(0L, length(relationships))
  names(ret) <- relationships
  
  # If counts provided, then generate pXdraws and pYdraws with Dirichlet posterior assuming iid sampling
  if (!is.null(Xcounts) && !is.null(Ycounts)) {
    J <- length(Xcounts)
    if (length(Ycounts) != J) stop("Xcounts and Ycounts must have the same number of categories.")
    if (J<=1) stop("Should have at least length(Xcounts)=length(Ycounts)=2 categories, usually more.")
    
    # Set RNG seed for replicability (and save current seed to re-seed on exit)
    oldseed <- NULL
    if (exists(".Random.seed", .GlobalEnv)) {  #.Random.seed #restore state at end
      oldseed <- get(".Random.seed", .GlobalEnv)
    }
    on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
    set.seed(112358) #for replicability
    
    # Posteriors for (p1,...,pJ) based on Dirichlet-multinomial model
    #   with uninformative (flat/uniform) Dir(1,...,1) prior.
    pXdraws <- Dir.draws(n=DRAWS, paramvec=rep(1,J)+Xcounts)
    pYdraws <- Dir.draws(n=DRAWS, paramvec=rep(1,J)+Ycounts)
  }
  
  
  if (J==2) {
    # Only SD1 possible, and no adjustment either way
    ret['X.SD1.Y'] <- mean(pXdraws[,1]<=pYdraws[,1])
    ret['Y.SD1.X'] <- mean(pYdraws[,1]<=pXdraws[,1])
    return(ret)
  }
  
  for (d in 1:DRAWS) {
    FX <- cumsum(pXdraws[d,]);  FY <- cumsum(pYdraws[d,])
    FX[J] <- FY[J] <- 1 #avoid rounding issues
    FXleFY <- (FX<=FY);  FYleFX <- (FY<=FX)
    medX <- min(which(FX>=0.5));  medY <- min(which(FY>=0.5))
    # 
    if (all(FXleFY) && all(FYleFX)) {
      # This is probability zero event, so doesn't really matter...
      ret['X.SD1.Y'] <- ret['X.SD1.Y']+1
      ret['Y.SD1.X'] <- ret['Y.SD1.X']+1
    } else if (all(FXleFY)) {
      ret['X.SD1.Y'] <- ret['X.SD1.Y']+1
    } else if (all(FYleFX)) {
      ret['Y.SD1.X'] <- ret['Y.SD1.X']+1
    } else if (suppressWarnings(max(which(!FYleFX))<min(which(!FXleFY)))) {
      ret['X.SC.Y'] <- ret['X.SC.Y']+1
      if (medY==medX) {
        if (all(ifelse(1:J<medX, FX<=FY, FY<=FX))) {
          ret['Y.MPS.X'] <- ret['Y.MPS.X']+1
        }
      }
    } else if (suppressWarnings(max(which(!FXleFY))<min(which(!FYleFX)))) {
      ret['Y.SC.X'] <- ret['Y.SC.X']+1
      if (medY==medX) {
        if (all(ifelse(1:J<medX, FYleFX, FXleFY))) {
          ret['X.MPS.Y'] <- ret['X.MPS.Y']+1
        }
      }
    } else {
      ret['multiple.crossings'] <- ret['multiple.crossings']+1
    }
  }
  
  ret <- ret / DRAWS
  
  if (adj.prior) {
    # Pasted from print.adj.priors() output
    pi0.over.pi1 <- rbind(
      c(2, 1.0000 , 1.0000 , 0.0000000000 , 0.0000000000 , 0.0000000000 , 0.0000000000  , 0.0000000000 ),
      c(3, 0.5107123300 , 0.5107123300 , 0.1932754879 , 0.1932754879 , 0.0000000000 , 0.0657572326 , 0.0657572326  ),
      c(4, 0.3270781819 , 0.3270781819 , 0.2542417376 , 0.2542417376 , 0.1132138484 , 0.0494866696 , 0.0494866696  ),
      c(5, 0.2520478240 , 0.2520478240 , 0.2482109988 , 0.2482109988 , 0.2495314257 , 0.0317789624 , 0.0317789624  ),
      c(6, 0.2029616386 , 0.2029616386 , 0.2366393933 , 0.2366393933 , 0.3886960144 , 0.0231755143 , 0.0231755143  ),
      c(7, 0.1650940431 , 0.1650940431 , 0.2168417292 , 0.2168417292 , 0.5629884339 , 0.0151768972 , 0.0151768972  ),
      c(8, 0.1416203839 , 0.1416203839 , 0.1992135437 , 0.1992135437 , 0.7232465966 , 0.0119411617 , 0.0119411617  ),
      c(9, 0.1254312971 , 0.1254312971 , 0.1813421789 , 0.1813421789 , 0.8871485186 , 0.0080138345 , 0.0080138345  ),
      c(10, 0.1127218979 , 0.1127218979 , 0.1672096843 , 0.1672096843 , 1.0445716622 , 0.0054299688 , 0.0054299688  ),
      c(11, 0.1029616159 , 0.1029616159 , 0.1642492382 , 0.1642492382 , 1.1331058020 , 0.0058846300 , 0.0058846300  ),
      c(12, 0.0928370085 , 0.0928370085 , 0.1508899509 , 0.1508899509 , 1.3142791021 , 0.0039675086 , 0.0039675086  ),
      c(13, 0.0836640429 , 0.0836640429 , 0.1421448921 , 0.1421448921 , 1.4795437639 , 0.0035627106 , 0.0035627106  ),
      c(14, 0.0777044540 , 0.0777044540 , 0.1309679279 , 0.1309679279 , 1.6609898882 , 0.0028081141 , 0.0028081141  ),
      c(15, 0.0738935768 , 0.0738935768 , 0.1230298957 , 0.1230298957 , 1.8034763106 , 0.0023053425 , 0.0023053425  ),
      c(16, 0.0661549881 , 0.0661549881 , 0.1173842110 , 0.1173842110 , 1.9922202274 , 0.0022048507 , 0.0022048507  ),
      c(17, 0.0616279804 , 0.0616279804 , 0.1114908873 , 0.1114908873 , 2.1575623619 , 0.0019538326 , 0.0019538326  ),
      c(18, 0.0584459799 , 0.0584459799 , 0.1051555540 , 0.1051555540 , 2.3255736615 , 0.0016527295 , 0.0016527295  ),
      c(19, 0.0545205940 , 0.0545205940 , 0.0981770394 , 0.0981770394 , 2.5435861091 , 0.0016025742 , 0.0016025742  ),
      c(20, 0.0546315141 , 0.0546315141 , 0.0943921919 , 0.0943921919 , 2.6218761318 , 0.0011513265 , 0.0011513265  ),
      c(21, 0.0497669889 , 0.0497669889 , 0.0905127484 , 0.0905127484 , 2.8343558282 , 0.0011012214 , 0.0011012214  ),
      c(22, 0.0478344352 , 0.0478344352 , 0.0883259702 , 0.0883259702 , 2.9432176656 , 0.0010511062 , 0.0010511062  ),
      c(23, 0.0457517741 , 0.0457517741 , 0.0841291763 , 0.0841291763 , 3.1203131438 , 0.0011012214 , 0.0011012214  ),
      c(24, 0.0420492980 , 0.0420492980 , 0.0816703339 , 0.0816703339 , 3.3159257661 , 0.0005503252 , 0.0005503252  ),
      c(25, 0.0406915711 , 0.0406915711 , 0.0774168120 , 0.0774168120 , 3.5065344750 , 0.0009509059 , 0.0009509059  ),
      c(26, 0.0389132623 , 0.0389132623 , 0.0761955856 , 0.0761955856 , 3.6189376443 , 0.0006003702 , 0.0006003702  ),
      c(27, 0.0383144948 , 0.0383144948 , 0.0709508401 , 0.0709508401 , 3.8473097431 , 0.0007505855 , 0.0007505855  ),
      c(28, 0.0357890906 , 0.0357890906 , 0.0690045921 , 0.0690045921 , 4.0454086781 , 0.0003501250 , 0.0003501250  ),
      c(29, 0.0340726981 , 0.0340726981 , 0.0657573295 , 0.0657573295 , 4.2826201796 , 0.0004001601 , 0.0004001601  ),
      c(30, 0.0343940544 , 0.0343940544 , 0.0633253223 , 0.0633253223 , 4.3879310345 , 0.0005002601 , 0.0005002601  )
    )
    dimnames(pi0.over.pi1) <- list(NULL, c('J', 'X.SD1.Y', 'Y.SD1.X', 'X.SC.Y', 'Y.SC.X', 'multiple.crossings', 'X.MPS.Y', 'Y.MPS.X'))
    tmp <- (pi0.over.pi1[,'J']==J)
    if (!any(tmp)) {
      adj <- adj.prior.sim(J=J)
    } else {
      adj <- pi0.over.pi1[tmp, ]
    }
    for (rel in relationships) {
      ret[rel] <- ret[rel] / (ret[rel] + adj[rel]*(1-ret[rel]))
    }
    if (J==3) ret['multiple.crossings'] <- 0 #impossible...
  }

  return(ret)
}


# Draw n vectors from a Dirichlet(paramvec) distribution
Dir.draws <- function(n=1, paramvec) {
  J <- length(paramvec)
  Gmat <- matrix(NA, nrow=n, ncol=J)
  for (j in 1:J) Gmat[,j] <- rgamma(n=n, shape=paramvec[j], scale=1)
  return(t(apply(X=Gmat, MARGIN=1, FUN=function(row)row/sum(row))))
}


# Simulate posterior adjustments to get prior P(H0)=1/2
adj.prior.sim <- function(J=2:10, DRAWS=10000L) {
  ret <- NULL
  for (ncat in J) {
    Xcounts <- Ycounts <- rep(0L, ncat)
    tmp <- ordinal.inequality.posteriors(Xcounts, Ycounts, adj.prior=FALSE, DRAWS=DRAWS)
    ret <- rbind(ret, tmp)
  }
  pi0.over.pi1 <- ret / (1-ret)
  pi0.over.pi1[,'X.SD1.Y'] <- pi0.over.pi1[,'Y.SD1.X'] <- 
    (pi0.over.pi1[,'X.SD1.Y'] + pi0.over.pi1[,'Y.SD1.X'])/2
  pi0.over.pi1[,'X.SC.Y'] <- pi0.over.pi1[,'Y.SC.X'] <- 
    (pi0.over.pi1[,'X.SC.Y'] + pi0.over.pi1[,'Y.SC.X'])/2
  pi0.over.pi1[,'X.MPS.Y'] <- pi0.over.pi1[,'Y.MPS.X'] <- 
    (pi0.over.pi1[,'X.MPS.Y'] + pi0.over.pi1[,'Y.MPS.X'])/2
  pi0.over.pi1 <- cbind(J=J, pi0.over.pi1)
  return(pi0.over.pi1)
}


# Print posterior adjustments for J=2..10
print.adj.priors <- function(J=2:30, DRAWS=10000L) {
  pi0.over.pi1 <- adj.prior.sim(J=J, DRAWS=DRAWS)
  cat("pi0.over.pi1 <- rbind(\n")
  cat("c(")
  for (j in 1:dim(pi0.over.pi1)[1]) {
    cat(sprintf("%d", pi0.over.pi1[j,1]))
    cat(sprintf(", %12.10f", pi0.over.pi1[j,-1]))
    if (j==dim(pi0.over.pi1)[1]) cat(")\n") else cat("),\nc(")
  }
  cat(")\nnames(pi0.over.pi1) <- c(")
  cat(paste0(sprintf("'%s'", dimnames(pi0.over.pi1)[[2]]), collapse=', '))
  cat(")\n")
}

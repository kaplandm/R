# Functions for finite-sample inference on bid distributions from auction data
#  - confidence interval for specific quantile
#  - uniform confidence band for CDF or quantile function
# Help with initial code: Ronald D. Flores
# Current code: David M. Kaplan & Xin Liu
# Feedback/support: https://kaplandm.github.io/


prob.loaded <- (exists("GK.dist.1s.alpha.tilde") && exists("quantile.inf"))
success <-
  tryCatch({source('https://raw.githubusercontent.com/kaplandm/R/main/GK_dist_inf.R');
    source('https://raw.githubusercontent.com/kaplandm/R/main/quantile_inf.R');
    TRUE},  error=function(w) FALSE)
if (!success) {
  success <-
    tryCatch({source("GK_dist_inf.R");source("quantile_inf.R");TRUE},
             error=function(w) FALSE)
  if (!success) {
    if (prob.loaded) {
      warning("Couldn't load GK_dist_inf.R and/or quantile_inf.R, but it seems like you already did, so continuing on.")
    } else {
      stop("Failed to source() GK_dist_inf.R and/or quantile_inf.R from web or local file.  You may download and source() it yourself, or at least make sure it's in your getwd().  Files currently available at https://github.com/kaplandm/R")
    }
  }
}


# Arguments for quantile confidence interval (CI)
# wins:    vector of winning bids, one from each auction
# price:   1 if first-price auction (or effectively so, like descending/Dutch)
#          2 if second-price (or effectively so, like ascending/English)
#          3 if third-price, etc.
# n:       number of bidders per auction (scalar, or if differs by auction then vector with same length as "wins")
# conf.level: confidence level as decimal; e.g., 0.90 for 90% CI
# tau:     quantile index of interest; e.g., 0.5 for median (of bid distribution)

# Return value:
# matrix with lower/upper/two-sided CIs for tau-quantile of bid distribution

bid.quantile.CI <- function(wins, price, n, conf.level, tau, NREP=1e3) {
  # Process arguments
  ALPHA <- 1 - conf.level
  wins <- sort(wins) # so wins[K] is the Kth order statistic
  if (!is.numeric(n) || !is.numeric(price) || !is.numeric(wins) || !is.numeric(conf.level) || !is.numeric(tau) || !is.numeric(NREP)) stop("All arguments must be numeric")
  J <- length(wins)
  if (length(n)>1 && length(n)!=J) stop("Argument n must be scalar or else vector with length J")
  if (tau<=0 || tau>=1) stop("Argument tau must satisfy 0<tau<1, like tau=0.5 for the median")
  if (conf.level<=0 || conf.level>=1) stop("Argument conf.level must satisfy 0<conf.level<1, like conf.level=0.90 for 90% confidence level")

  # Set up confidence interval matrix
  CIs <- matrix(data=NA, nrow=3, ncol=2,
                dimnames=list(c('Near-exact lower CI','Near-exact upper CI',
                                'Near-exact two-sided CI'),
                              c('Lower','Upper')))
  CIs[1,1] <- -Inf
  CIs[2,2] <- +Inf

  if (length(n)==1) { # Use Goldman and Kaplan (2017)
    k <- n + 1 - price  #order statistic index of winning bid
    CIs[1,2] <- quantile.inf(X=wins, p=pbeta(tau, k, n+1-k), ALPHA=ALPHA, ONESIDED=-1)$CI$hi
    CIs[2,1] <- quantile.inf(X=wins, p=pbeta(tau, k, n+1-k), ALPHA=ALPHA, ONESIDED= 1)$CI$lo
    CIs[3, ] <- unlist(quantile.inf(X=wins, p=pbeta(tau, k, n+1-k), ALPHA=ALPHA, ONESIDED= 0)$CI)
  } else { # Simulate inid Beta rv's
    rownames(CIs) <- c('Exact (simulated) conservative lower CI',
                       'Exact (simulated) conservative upper CI',
                       'Exact (simulated) conservative two-sided CI')
    # Save/set RNG seed
    oldseed <- NULL
    if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
      oldseed <- get(".Random.seed",.GlobalEnv)
    }
    on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
    set.seed(112358)
    # Simulate inid beta order statistics
    beta.sims <- NULL
    for (j in 1:J) {
      beta.sims <- cbind(beta.sims, rbeta(n=NREP, shape1=n[j]+1-price, shape2=price))
    }
    beta.sims <- t(apply(X=beta.sims, MARGIN=1, FUN=sort))
    # Compute probabilities and choose order statistics accordingly (erring on conservative)
    beta.lt.probs <- apply(X=beta.sims, MARGIN=2, FUN=function(col)mean(col<=tau))
    beta.gt.probs <- apply(X=beta.sims, MARGIN=2, FUN=function(col)mean(col>=tau))
    r <- suppressWarnings(min(which(beta.gt.probs>=conf.level)))
    CIs[1,2] <- ifelse(r<Inf, wins[r], Inf)
    r <- suppressWarnings(max(which(beta.lt.probs>=conf.level)))
    CIs[2,1] <- ifelse(r>-Inf, wins[r], -Inf)
    if (CIs[1,2]<Inf && CIs[2,1]>-Inf) {
      r1a <- max(which(beta.lt.probs>=1-ALPHA/2))
      beta1a <- beta.sims[,r1a]
      beta.2a.probs <- apply(X=beta.sims[,-(1:r1a)], MARGIN=2, FUN=function(col)mean(beta1a<=tau & col>=tau))
      r2a <- r1a + min(which(beta.2a.probs>=conf.level))
      r2b <- min(which(beta.gt.probs>=1-ALPHA/2))
      beta2b <- beta.sims[,r2b]
      beta.2b.probs <- apply(X=beta.sims[,1:(r2b-1)], MARGIN=2, FUN=function(col)mean(beta2b>=tau & col<=tau))
      r1b <- max(which(beta.2b.probs>=conf.level))
      if (r1a>-Inf && r1b<Inf) {
        if (r2a>-Inf && r2b<Inf) {
          if (r2a-r1a < r2b-r1b) {
            CIs[3,] <- c(wins[r1a], wins[r2a])
          } else {
            CIs[3,] <- c(wins[r1b], wins[r2b])
          }
        } else {
          CIs[3,1] <- ifelse(r1a>-Inf, wins[r1a], -Inf)
          CIs[3,2] <- ifelse(r2a< Inf, wins[r2a],  Inf)
        }
      } else {
        CIs[3,1] <- ifelse(r1b>-Inf, wins[r1b], -Inf)
        CIs[3,2] <- ifelse(r2b< Inf, wins[r2b],  Inf)
      }
    } else {
      CIs[3,] <- c(-Inf,Inf)
    }
  }
  return(CIs)
}


# # Examples
# set.seed(112358)
# J <- 100; n <- 10; price <- 2
# fake.bids <- matrix(data=rnorm(J*n, mean=100, sd=10), nrow=J)
# fake.wins <- apply(X=fake.bids, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
# bid.quantile.CI(wins=fake.wins, price=price, n=n, conf.level=0.90, tau=0.75)
# bid.quantile.CI(wins=fake.wins, price=price, n=n, conf.level=0.90, tau=0.75)[3,] #select one single type of CI
# bid.quantile.CI(wins=fake.wins, price=price, n=n, conf.level=0.90, tau=0.50) #not enough data for lower endpoints for median, but can still get upper endpoints
# bid.quantile.CI(wins=fake.wins, price=price, n=rep(10,J), conf.level=0.90, tau=0.75) #vector n, but all n=10
# bid.quantile.CI(wins=fake.wins, price=price, n=sample(x=9:11,size=J,replace=T), conf.level=0.90, tau=0.75) #vector n, varying n[j]
# #
# (st <- Sys.time())
# set.seed(112358)
# NREP <- 1000
# J1 <- J2 <- 20
# n1 <- 4;  n2 <- 10
# price <- 2
# ALPHA <- 0.1
# FB <- punif;  QB <- qunif
# TAU <- 0.75
# covers <- NULL
# for (irep in 1:NREP) {
#   fake.bids1 <- matrix(data=runif(J1*n1), nrow=J1)
#   fake.wins1 <- apply(X=fake.bids1, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
#   fake.bids2 <- matrix(data=runif(J2*n2), nrow=J2)
#   fake.wins2 <- apply(X=fake.bids2, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
#   ret <- bid.quantile.CI(wins=c(fake.wins1,fake.wins2), price=price, n=c(rep(n1,J1),rep(n2,J2)), conf.level=0.90, tau=TAU, NREP=5e3)
#   covers <- rbind(covers, c(QB(TAU)<=ret[1,2], ret[2,1]<=QB(TAU), ret[3,1]<=QB(TAU) && QB(TAU)<=ret[3,2]))
# }
# colMeans(covers) # should be around 1-ALPHA (some simulation error; less if increase NREP)
# Sys.time() - st # Around 4 minutes
# Sys.time()


# Function to compute uniform confidence band for CDF of bid distribution
# wins, price, n: same as for bid.quantile.CI above
# conf.level: now the uniform confidence level, like 0.90 for 90%
# plot.type: NULL for no plot, '1s' for one-sided bands, '2s' for two-sided bands
# plot.quantile: TRUE to plot quantile function bands, FALSE for CDF bands (if non-NULL plot.type)
# NREP: number of replications for simulation (higher = more accurate); only applies to vector n
# pt.beta.flag: use pointwise (average) beta distributions, or pointwise simulated distributions?
# ALPHA.pts: if the pointwise alpha values have been simulated previously (or you just want to try other values), can specify as a vector here
# return.ALPHA.pts: (for vector n only) add another column to the return matrix, whose first three entries are the pointwise alpha values
bid.band <- function(wins, price, n, conf.level, plot.type='1s', plot.quantile=TRUE, NREP=1e3, pt.beta.flag=TRUE, ALPHA.pts=NULL, return.ALPHA.pts=FALSE) {
  # Process arguments
  ALPHA.joint <- 1 - conf.level
  wins <- sort(wins) # so wins[K] is the Kth order statistic
  if (!is.numeric(n) || !is.numeric(price) || !is.numeric(wins) || !is.numeric(conf.level) || !is.numeric(NREP)) stop("All arguments besides plot.type must be numeric")
  J <- length(wins)
  if (length(n)>1 && length(n)!=J) stop("Argument n must be scalar or else vector with length J")
  if (conf.level<=0 || conf.level>=1) stop("Argument conf.level must satisfy 0<conf.level<1, like conf.level=0.90 for 90% confidence level")
  if (!is.null(ALPHA.pts) && (!is.numeric(ALPHA.pts) || length(ALPHA.pts)<2 || length(ALPHA.pts)>3)) stop("Argument ALPHA.pts must be either NULL or a numeric vector of length 2 (for scalar n) or 3 (for vector n)")

  if (length(n)==1) { # Use Goldman and Kaplan (2018)
    k <- n + 1 - price  #order statistic index of winning bid
    if (is.null(ALPHA.pts)) {
      ALPHA.pt.1s <- GK.dist.1s.alpha.tilde(n=J, alpha=2*ALPHA.joint-ALPHA.joint^2) / 2
      ALPHA.pt.2s <- GK.dist.1s.alpha.tilde(n=J, alpha=ALPHA.joint) / 2
    } else {
      ALPHA.pt.1s <- ALPHA.pts[1]
      ALPHA.pt.2s <- ALPHA.pts[2]
    }
    y.pts.lo1 <- qbeta( qbeta(  ALPHA.pt.1s, 1:J, J:1),  k, n+1-k)
    y.pts.up1 <- qbeta( qbeta(1-ALPHA.pt.1s, 1:J, J:1),  k, n+1-k)
    y.pts.lo2 <- qbeta( qbeta(  ALPHA.pt.2s, 1:J, J:1),  k, n+1-k)
    y.pts.up2 <- qbeta( qbeta(1-ALPHA.pt.2s, 1:J, J:1),  k, n+1-k)
    ret <- cbind(winning.bids=wins, UCB.lower.1s=y.pts.lo1, UCB.upper.1s=y.pts.up1,
                 UCB.lower.2s=y.pts.lo2, UCB.upper.2s=y.pts.up2)
  } else { # Simulate inid Beta rv's
    if (pt.beta.flag) {
      n.avg <- mean(n) # Only for parameterization; does not affect uniform coverage probability
      k.avg <- n.avg + 1 - price
      # Alternative: "average CDF" like in David & Nagaraja Thm 5.2.1
      counts <- table(n)
      vals <- as.numeric(names(counts))
      Fbar <- function(x) sum(counts*pbeta(x, shape1=vals+1-price, shape2=price)) / sum(counts)
    }
    if (!pt.beta.flag || is.null(ALPHA.pts)) {
      # Save/set RNG seed
      oldseed <- NULL
      if (exists(".Random.seed",.GlobalEnv)) {  #.Random.seed #restore state at end
        oldseed <- get(".Random.seed",.GlobalEnv)
      }
      on.exit(if (!is.null(oldseed)) { assign(".Random.seed", oldseed, .GlobalEnv) }, add=TRUE)
      set.seed(112358)
      # 
      # Simulate inid beta order statistics
      beta.sims <- NULL
      for (j in 1:J) {
        beta.sims <- cbind(beta.sims, rbeta(n=NREP, shape1=n[j]+1-price, shape2=price))
      }
      beta.sims <- t(apply(X=beta.sims, MARGIN=1, FUN=sort))
    }
    if (is.null(ALPHA.pts)) {
      # Solve for pointwise ALPHA.pt that achieves joint probability equivalent to 1-ALPHA.joint
      ALPHA.check.fn.lo1 <- function(a) {
        if (pt.beta.flag) {
          # Use average k and n (not particularly good)
          # beta.OS.quantiles <- qbeta( qbeta(  a, 1:J, J:1),  k.avg, n.avg+1-k.avg)
          # Use average CDF, and invert to get quantile
          beta.OS.quantiles <- rep(NA, J)
          for (j in 1:J) beta.OS.quantiles[j] <- uniroot(f=function(q)Fbar(q)-qbeta(  a, j, J+1-j), interval=0:1)$root
        } else {
          beta.OS.quantiles <- apply(X=beta.sims, MARGIN=2, FUN=function(col)quantile(x=col, probs=a, type=6))
        }
        pr.joint <- mean(apply(X=beta.sims, MARGIN=1,
                               FUN=function(b) all(b>=beta.OS.quantiles) ))
        return(pr.joint - (1-ALPHA.joint))
      }
      ALPHA.check.fn.up1 <- function(a) {
        if (pt.beta.flag) {
          beta.OS.quantiles <- rep(NA, J)
          for (j in 1:J) beta.OS.quantiles[j] <- uniroot(f=function(q)Fbar(q)-qbeta(1-a, j, J+1-j), interval=0:1)$root
        } else {
          beta.OS.quantiles <- apply(X=beta.sims, MARGIN=2, FUN=function(col)quantile(x=col, probs=1-a, type=6))
        }
        pr.joint <- mean(apply(X=beta.sims, MARGIN=1,
                               FUN=function(b) all(b<=beta.OS.quantiles) ))
        return(pr.joint - (1-ALPHA.joint))
      }
      ALPHA.check.fn.2s <- function(a) {
        if (pt.beta.flag) {
          beta.OS.quantiles.L <- beta.OS.quantiles.U <- rep(NA, J)
          for (j in 1:J) {
            beta.OS.quantiles.L[j] <- uniroot(f=function(q)Fbar(q)-qbeta(  a, j, J+1-j), interval=0:1)$root
            beta.OS.quantiles.U[j] <- uniroot(f=function(q)Fbar(q)-qbeta(1-a, j, J+1-j), interval=0:1)$root
          }
        } else {
          beta.OS.quantiles.L <- apply(X=beta.sims, MARGIN=2, FUN=function(col)quantile(x=col, probs= a , type=6))
          beta.OS.quantiles.U <- apply(X=beta.sims, MARGIN=2, FUN=function(col)quantile(x=col, probs=1-a, type=6))
        }
        all2s <- function(b) all(beta.OS.quantiles.L<=b) && all(b<=beta.OS.quantiles.U)
        pr.joint <- mean( apply(X=beta.sims, MARGIN=1, FUN=all2s) )
        return(pr.joint - (1-ALPHA.joint))
      }
      # If uniroot fails to get close enough (note: uniroot tol is for the root, not the function value)
      alpha.backup.search <- function(fn, init.raw, ftol, lo=0, hi=1, mult=1e3) {
        cur <- init.raw   # / (mult*0.8)
        cur.f <- fn(cur)
        while (abs(cur.f)>ftol && cur>1e-33) {
          if (cur.f>0) lo <- cur else hi <- cur
          if (lo==0) {
            cur <- cur / mult
            if (cur<1e-12) ftol <- ALPHA.joint / 20
          } else if (hi==1) {
            cur <- cur <- sqrt(lo*hi)   # cur * mult
          } else cur <- sqrt(lo*hi)
          cur.f <- fn(cur)
        }
        if (abs(cur.f)>ftol) warning(sprintf("Numerical failure: coverage probability may have error of %g",abs(cur.f)))
        return(cur)
      }
      ftol <- ALPHA.joint / 10;  UNITOL <- 1e-5;  UNITOL2 <- 1e-9
      tmp <- uniroot(ALPHA.check.fn.lo1, interval=0:1, tol=UNITOL)
      if (abs(tmp$f.root)<ftol) {
        ALPHA.pt.lo1 <- tmp$root
      } else {
        tmp <- uniroot(ALPHA.check.fn.lo1, interval=0:1, tol=UNITOL2)
        if (abs(tmp$f.root)<ftol) {
          ALPHA.pt.lo1 <- tmp$root
        } else {
          ALPHA.pt.lo1 <- tmp$root - ifelse(tmp$f.root<0, tmp$estim.prec, 0)
        }
      }
      tmp <- uniroot(ALPHA.check.fn.up1, interval=0:1, tol=UNITOL)
      if (abs(tmp$f.root)<ftol) {
        ALPHA.pt.up1 <- tmp$root
      } else {
        tmp <- uniroot(ALPHA.check.fn.up1, interval=0:1, tol=UNITOL2)
        if (abs(tmp$f.root)<ftol) {
          ALPHA.pt.up1 <- tmp$root
        } else {
          ALPHA.pt.up1 <- tmp$root - ifelse(tmp$f.root<0, tmp$estim.prec, 0)
        }
      }
      tmp <- uniroot(ALPHA.check.fn.2s, interval=0:1, tol=UNITOL)
      if (abs(tmp$f.root)<ftol) {
        ALPHA.pt.2s <- tmp$root
      } else {
        tmp <- uniroot(ALPHA.check.fn.2s, interval=0:1, tol=UNITOL2)
        if (abs(tmp$f.root)<ftol) {
          ALPHA.pt.2s <- tmp$root
        } else {
          ALPHA.pt.2s <- tmp$root - ifelse(tmp$f.root<0, tmp$estim.prec, 0)
        }
      }
    } else {
      ALPHA.pt.lo1 <- ALPHA.pts[1]
      ALPHA.pt.up1 <- ALPHA.pts[2]
      ALPHA.pt.2s  <- ALPHA.pts[3]
    }
    # print(c(ALPHA.pt.lo1, ALPHA.pt.up1, ALPHA.pt.2s))
    if (pt.beta.flag) {
      y.pts.lo1 <- rep(NA, J)
      for (j in 1:J) y.pts.lo1[j] <- uniroot(f=function(q)Fbar(q)-qbeta(  ALPHA.pt.lo1, j, J+1-j), interval=0:1)$root
      y.pts.up1 <- rep(NA, J)
      for (j in 1:J) y.pts.up1[j] <- uniroot(f=function(q)Fbar(q)-qbeta(1-ALPHA.pt.lo1, j, J+1-j), interval=0:1)$root
      y.pts.lo2 <- y.pts.up2 <- rep(NA, J)
      for (j in 1:J) {
        y.pts.lo2[j] <- uniroot(f=function(q)Fbar(q)-qbeta(  ALPHA.pt.2s, j, J+1-j), interval=0:1)$root
        y.pts.up2[j] <- uniroot(f=function(q)Fbar(q)-qbeta(1-ALPHA.pt.2s, j, J+1-j), interval=0:1)$root
      }
    } else {
      y.pts.lo1 <- apply(X=beta.sims, MARGIN=2, FUN=function(col)quantile(x=col, probs=  ALPHA.pt.lo1, type=6))
      y.pts.up1 <- apply(X=beta.sims, MARGIN=2, FUN=function(col)quantile(x=col, probs=1-ALPHA.pt.up1, type=6))
      y.pts.lo2 <- apply(X=beta.sims, MARGIN=2, FUN=function(col)quantile(x=col, probs=  ALPHA.pt.2s , type=6))
      y.pts.up2 <- apply(X=beta.sims, MARGIN=2, FUN=function(col)quantile(x=col, probs=1-ALPHA.pt.2s , type=6))
    }
    ret <- cbind(winning.bids=wins, UCB.lower.1s=y.pts.lo1, UCB.upper.1s=y.pts.up1,
                 UCB.lower.2s=y.pts.lo2, UCB.upper.2s=y.pts.up2)
    if (return.ALPHA.pts) ret <- cbind(ret, ALPHA.pts=c(ALPHA.pt.lo1, ALPHA.pt.up1, ALPHA.pt.2s, rep(NA,length(wins)-3)))
  }

  # Plot if desired
  if (!is.null(plot.type) && !is.na(plot.type)) {
    if (is.null(plot.quantile) || !is.logical(plot.quantile) || is.na(plot.quantile)) {
      warning('Argument plot.quantile should be TRUE or FALSE; setting to TRUE and proceeding')
      plot.quantile <- TRUE
    }
    if (plot.type=='1s') {
      plot.bid.band(band=ret, conf.level=conf.level, two.sided=FALSE, quantile.band=plot.quantile)
    } else if (plot.type=='2s') {
      plot.bid.CDF.band(band=ret, conf.level=conf.level, two.sided=TRUE, quantile.band=plot.quantile)
    } else if (!is.null(plot.type) && !is.na(plot.type) && plot.type!='') warning("Argument plot.type is neither '1s' nor '2s' so is ignored (no plot).")
  }
  return(ret)
}


# Plot uniform confidence band for either quantile function or CDF of bid distribution
# Feel free to customize (or email to request customization)
# band: the return value from calling bid.CDF.band()
# conf.level needs to be provided separately (e.g., 0.90 for 90% confidence level)
# two.sided: FALSE for one-sided, TRUE for two-sided (slightly wider but usually pretty similar)
# main: title for plot
# plot.legend: include a legend? (TRUE/FALSE)
# mar: margins argument for par()
# ...: any other arguments to pass along to plot(), like xlim, yaxs, etc.
plot.bid.band <- function(band, conf.level=NA, two.sided=FALSE, quantile.band=TRUE,
                          main=NULL, plot.legend=TRUE, col=1, add=FALSE, mar=c(5.1,4.1,6.1,2.1), ...) {
  par(mar=mar)
  if (is.null(main)) {
    main <- sprintf('%s uniform confidence band%s\n',
                    ifelse(is.na(conf.level),'',sprintf('%g%%',100*conf.level)),
                    ifelse(two.sided,'','s'))
  }
  wins <- band[,'winning.bids']
  tau.lo1 <- band[,'UCB.lower.1s']
  tau.lo2 <- band[,'UCB.lower.2s']
  tau.up1 <- band[,'UCB.upper.1s']
  tau.up2 <- band[,'UCB.upper.2s']
  if (!add) { # create new plot
    if (quantile.band) {
      plot(x=range(tau.up1,tau.lo1,tau.up2,tau.lo2), y=range(wins), type='n',
           xlab=expression(Quantile~(tau)), ylab="Bid value", xlim=0:1, main=main, xaxt='n', ...)
      axis(side=1, at=0:4/4)
    } else {
      plot(x=range(wins), y=range(tau.up1,tau.lo1,tau.up2,tau.lo2), type='n',
           xlab="Bid value", ylab="Bid CDF", ylim=0:1, main=main, yaxt='n', ...)
      axis(side=2, at=0:4/4)
    }
  }
  bignum <- 100*(max(wins)-min(wins))
  if (two.sided) {
    x <- c(wins[1]-bignum,wins[1],wins,max(wins)+bignum)
    y <- c(0,0,tau.lo2,max(tau.lo2))
    if (quantile.band) { z <- y; y <- x; x <- z }
    lines(x=x, y=y, type='s', lwd=2, lty=2, col=col)
    x <- c(wins[1]-bignum,wins,max(wins),max(wins)+bignum)
    y <- c(min(tau.up2),tau.up2,1,1)
    if (quantile.band) { z <- y; y <- x; x <- z }
    lines(x=x, y=y, type='S', lwd=2, lty=2, col=col)
    if (plot.legend) legend(ifelse(quantile.band,'topleft','bottomright'), legend='Two-sided', lty=2, lwd=2, col=col)
  } else {
    x <- c(wins[1]-bignum,wins[1],wins,max(wins)+bignum)
    y <- c(0,0,tau.lo1,max(tau.lo1))
    if (quantile.band) { z <- y; y <- x; x <- z }
    lines(x=x, y=y, type='s', lwd=2, lty=2, col=col)
    x <- c(wins[1]-bignum,wins,max(wins),max(wins)+bignum)
    y <- c(min(tau.up1),tau.up1,1,1)
    if (quantile.band) { z <- y; y <- x; x <- z }
    lines(x=x, y=y, type='S', lwd=2, lty=2, col=col)
    if (plot.legend) {
      tmp <- legend('top', plot=FALSE, lty=2, lwd=2, col=col, ncol=2,
                    legend=c('Lower 1-sided','Upper 1-sided'))
      newxy <- c(tmp$rect$left, tmp$rect$top + 1.01*tmp$rect$h)
      legend(x=newxy[1], y=newxy[2], legend=c('Lower 1-sided','Upper 1-sided'),
             lty=2, lwd=2, col=col, ncol=2, xpd=NA)
    }
  }
}


# # Examples
# set.seed(112358)
# J <- 100; n <- 10; price <- 2
# fake.bids <- matrix(data=rnorm(J*n, mean=100, sd=10), nrow=J)
# fake.wins <- apply(X=fake.bids, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
# ret <- bid.CDF.band(wins=fake.wins, price=price, n=n, conf.level=0.90)
# #
# (st <- Sys.time())
# set.seed(112358)
# J <- 100; n <- 10; price <- 2
# NREP <- 1000 #simulation replications
# ALPHA <- 0.1
# FB <- punif
# covers <- covers2 <- NULL
# for (irep in 1:NREP) {
#   fake.bids <- matrix(data=runif(J*n), nrow=J)
#   fake.wins <- apply(X=fake.bids, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
#   ret <- bid.CDF.band(wins=fake.wins, price=price, n=n, conf.level=1-ALPHA, plot.type=NULL)
#   cover1L <- all(ret[,2]<=FB(ret[,1]))
#   cover1U <- all(FB(ret[,1])<=ret[,3])
#   cover2  <- all(ret[,4]<=FB(ret[,1]), FB(ret[,1])<=ret[,5])
#   covers  <- rbind(covers, c(cover1L, cover1U, cover2))
#   ret <- bid.CDF.band(wins=fake.wins, price=price, n=rep(n,J), conf.level=1-ALPHA, plot.type=NULL, NREP=1e3)
#   cover1L <- all(ret[,2]<=FB(ret[,1]))
#   cover1U <- all(FB(ret[,1])<=ret[,3])
#   cover2  <- all(ret[,4]<=FB(ret[,1]), FB(ret[,1])<=ret[,5])
#   covers2  <- rbind(covers2, c(cover1L, cover1U, cover2))
# }
# rbind( colMeans(covers), colMeans(covers2) ) # Output: 0.903 0.922 0.910 \\ 0.904 0.918 0.914
# Sys.time() - st # Around 3-4 minutes
# Sys.time()
# #
# print(st <- Sys.time())
# pt.beta.flag <- TRUE
# set.seed(112358)
# NREP <- 100 #simulation replications
# beta.sim.rep <- 1e4 #draws from beta order statistic distribution
# J1 <- J2 <- 20
# n1 <- 4;  n2 <- 10
# price <- 2
# ALPHA <- 0.1
# FB <- punif
# covers <- NULL
# A.pts <- NULL
# for (irep in 1:NREP) {
#   fake.bids1 <- matrix(data=runif(J1*n1), nrow=J1)
#   fake.wins1 <- apply(X=fake.bids1, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
#   fake.bids2 <- matrix(data=runif(J2*n2), nrow=J2)
#   fake.wins2 <- apply(X=fake.bids2, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
#   if (irep==1) {
#     ret <- bid.CDF.band(wins=c(fake.wins1,fake.wins2), price=price, n=c(rep(n1,J1),rep(n2,J2)), conf.level=1-ALPHA, plot.type=NULL, NREP=beta.sim.rep, pt.beta.flag=pt.beta.flag, return.ALPHA.pts=TRUE)
#     A.pts <- ret[1:3,'ALPHA.pts']
#     print(sprintf("A.pts=%g",A.pts)) # original with pt.beta.flag=TRUE: [1] "A.pts=1.28295e-06" "A.pts=0.0108473"   "A.pts=1.28338e-06"
#   } else {
#     ret <- bid.CDF.band(wins=c(fake.wins1,fake.wins2), price=price, n=c(rep(n1,J1),rep(n2,J2)), conf.level=1-ALPHA, plot.type=NULL, NREP=beta.sim.rep, pt.beta.flag=pt.beta.flag, ALPHA.pts=A.pts)
#   }
#   cover1L <- all(ret[,2]<=FB(ret[,1]))
#   cover1U <- all(FB(ret[,1])<=ret[,3])
#   cover2  <- all(ret[,4]<=FB(ret[,1]), FB(ret[,1])<=ret[,5])
#   covers  <- rbind(covers, c(cover1L, cover1U, cover2))
# }
# print(colMeans(covers))  # should be all near 0.90 (esp. w/ larger NREP)
# print(Sys.time() - st)
# print(Sys.time())
#
# set.seed(112358)
# J <- 10; n <- 4; price <- 2
# fake.bids <- matrix(data=rnorm(J*n, mean=100, sd=10), nrow=J)
# fake.wins <- apply(X=fake.bids, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
# ret <- bid.CDF.band(wins=fake.wins, price=price, n=n, conf.level=0.90)

#EOF

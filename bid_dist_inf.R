# Functions for finite-sample inference on bid distributions from auction data
#  - confidence interval for specific quantile
#  - uniform confidence band for CDF
# Original: Ronald D. Flores (2020-2022)
# Updates: David M. Kaplan (July 2022)
# Feedback/support: https://kaplandm.github.io/ has my contact info


# Arguments for quantile confidence interval (CI)
  # J:       number of auctions
  # n:       number of bidders per auction
  # price:   1 if first-price auction (or effectively so, like descending/Dutch)
  #          2 if second-price (or effectively so, like ascending/English)
  #          3 if third-price, etc.
  # wins:    vector of winning bids from J auctions
  # conf.level: confidence level as decimal; e.g., 0.90 for 90% CI
  # tau:     quantile index of interest; e.g., 0.5 for median (of bid distribution)

# Return value:
  # matrix with rows for exact and conservative CIs, lower/upper/two-sided,
  # with associated coverage probabilities. 

bid.quantile.CI <- function(J, n, price, wins, conf.level, tau) {
  # Process arguments
  k <- n + 1 - price  #order statistic index of winning bid
  ALPHA <- 1 - conf.level
  wins <- sort(wins) #so wins[k] is the k'th order statistic

  # Solve for fractional order statistic indices r; 1=one-sided, 2=two-sided; lo=lower endpoint; up=upper endpoint
  r.lo1 <- tryCatch(expr=uniroot(function(x){pbeta(pbeta(tau, k, n+1-k), x, J+1-x) - (1-ALPHA)  }, interval=c(1, J))$root,
                    error=function(w)NA)
  r.lo2 <- tryCatch(uniroot(function(x){pbeta(pbeta(tau, k, n+1-k), x, J+1-x) - (1-ALPHA/2)}, interval=c(1, J))$root,
                    error=function(w)NA)
  r.up1 <- tryCatch(uniroot(function(x){pbeta(pbeta(tau, k, n+1-k), x, J+1-x) -    ALPHA   }, interval=c(1, J))$root,
                    error=function(w)NA)
  r.up2 <- tryCatch(uniroot(function(x){pbeta(pbeta(tau, k, n+1-k), x, J+1-x) -    ALPHA/2 }, interval=c(1, J))$root,
                    error=function(w)NA)

  # Set up confidence interval matrix
  CIs <- matrix(data=NA, nrow=6, ncol=2,
                dimnames=list(c('Exact lower CI','Conservative lower CI',
                                'Exact upper CI','Conservative upper CI',
                                'Exact two-sided CI','Conservative two-sided CI'),
                              c('Lower','Upper')))
  CIs[1:2,1] <- -Inf
  CIs[3:4,2] <- +Inf

  # Given indices, compute corresponding fractional order statistics from data,
  #   as well as certain fundamental probabilities used to compute coverage probabilities,
  #   and compute confidence intervals
  pr.tmp <- pbeta(pbeta(tau, k, n+1-k), 1:J, J:1)
  pr.r.up1.frac <- pr.r.up1.ceil <- pr.r.lo1.frac <- pr.r.lo1.floor <- 
    pr.r.up2.frac <- pr.r.lo2.frac <- pr.r.up2.ceil <- pr.r.lo2.floor <- NA
  if (!is.na(r.lo1)) { #one-sided lower endpoint (for upper CI)
    win.r.lo1.floor <- wins[floor(r.lo1)]
    win.r.lo1.ceil  <- wins[ceiling(r.lo1)]
    win.r.lo1.frac  <- win.r.lo1.floor + (r.lo1 - floor(r.lo1))*(win.r.lo1.ceil - win.r.lo1.floor)
    pr.r.lo1.floor <- pr.tmp[floor(r.lo1)]
    pr.r.lo1.frac  <- pbeta(pbeta(tau, k, n+1-k), r.lo1, J+1-r.lo1)
    CIs[3,1] <- win.r.lo1.frac
    CIs[4,1] <- win.r.lo1.floor
  }
  if (!is.na(r.lo2)) { #two-sided lower endpoint
    win.r.lo2.floor <- wins[floor(r.lo2)]
    win.r.lo2.ceil  <- wins[ceiling(r.lo2)]
    win.r.lo2.frac  <- win.r.lo2.floor + (r.lo2 - floor(r.lo2))*(win.r.lo2.ceil - win.r.lo2.floor)
    pr.r.lo2.floor <- pr.tmp[floor(r.lo2)]
    pr.r.lo2.frac  <- pbeta(pbeta(tau, k, n+1-k), r.lo2, J+1-r.lo2)
    CIs[5,1] <- win.r.lo2.frac
    CIs[6,1] <- win.r.lo2.floor
  }
  if (!is.na(r.up1)) { #one-sided upper endpoint (for lower CI)
    win.r.up1.floor <- wins[floor(r.up1)]
    win.r.up1.ceil  <- wins[ceiling(r.up1)]
    win.r.up1.frac  <- win.r.up1.floor + (r.up1 - floor(r.up1))*(win.r.up1.ceil - win.r.up1.floor)
    pr.r.up1.ceil  <- pr.tmp[ceiling(r.up1)]
    pr.r.up1.frac  <- pbeta(pbeta(tau, k, n+1-k), r.up1, J+1-r.up1)
    CIs[1,2] <- win.r.up1.frac
    CIs[2,2] <- win.r.up1.ceil
  }
  if (!is.na(r.up2)) { #two-sided upper endpoint
    win.r.up2.floor <- wins[floor(r.up2)]
    win.r.up2.ceil  <- wins[ceiling(r.up2)]
    win.r.up2.frac  <- win.r.up2.floor + (r.up2 - floor(r.up2))*(win.r.up2.ceil - win.r.up2.floor)
    pr.r.up2.ceil  <- pr.tmp[ceiling(r.up2)]
    pr.r.up2.frac  <- pbeta(pbeta(tau, k, n+1-k), r.up2, J+1-r.up2)
    CIs[5,2] <- win.r.up2.frac
    CIs[6,2] <- win.r.up2.ceil
  }

  # Compute corresponding CP
  CPs <- c(1 - pr.r.up1.frac, 1 - pr.r.up1.ceil, pr.r.lo1.frac, pr.r.lo1.floor,
           pr.r.lo2.frac - pr.r.up2.frac, pr.r.lo2.floor - pr.r.up2.ceil)

  # Return values
  ret <- cbind(CIs,CP=CPs)
  for (i in 1:nrow(ret)) {
    if (any(is.na(ret[i,]))) ret[i,] <- NA
  }
  return(ret)
}  


# # Examples
# set.seed(112358)
# J <- 100; n <- 10; price <- 2
# fake.bids <- matrix(data=rnorm(J*n, mean=100, sd=10), nrow=J)
# fake.wins <- apply(X=fake.bids, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
# bid.quantile.CI(J=J, n=n, price=price, wins=fake.wins, conf.level=0.90, tau=0.75)
# bid.quantile.CI(J=J, n=n, price=price, wins=fake.wins, conf.level=0.90, tau=0.75)[3,] #select one single type of CI
# bid.quantile.CI(J=J, n=n, price=price, wins=fake.wins, conf.level=0.90, tau=0.50) #not enough data for lower endpoints for median, but can still get upper endpoints



# Function to compute uniform confidence band for CDF of bid distribution
# J, n, price, wins: same as for bid.quantile.CI above
# conf.level: now the uniform confidence level, like 0.90 for 90%
# upper: TRUE for upper band, FALSE for lower
# NREP: number of replications for simulation
# plot.type: NULL for no plot, '1s' for one-sided bands, '2s' for two-sided bands
bid.CDF.band <- function(J, n, price, wins, conf.level, NREP=1e4, plot.type='1s') {
  # Process arguments
  k <- n + 1 - price  #order statistic index of winning bid
  ALPHA.joint <- 1 - conf.level
  wins <- sort(wins) #so wins[k] is the k'th order statistic
  
  # Simulate beta rv's and order statistics
  beta.J.r.matrix <- matrix(data=NA, nrow=NREP, ncol=J)
  for (i in 1:NREP) beta.J.r.matrix[i,] <- sort(rbeta(J, k, n+1-k))
  
  # Solve for pointwise ALPHA.pt that achieves joint probability equivalent to 1-ALPHA.joint
  ALPHA.check.fn.lo1 <- function(a) {
    beta.OS.quantiles <- qbeta( qbeta(  a, 1:J, J:1),  k, n+1-k)
    pr.joint <- mean(apply(beta.J.r.matrix, MARGIN=1,
                           FUN=function(b) all(b>=beta.OS.quantiles) ))
    return(pr.joint - (1-ALPHA.joint))
  }
  ALPHA.check.fn.up1 <- function(a) {
    beta.OS.quantiles <- qbeta( qbeta(1-a, 1:J, J:1),  k, n+1-k)
    pr.joint <- mean(apply(beta.J.r.matrix, MARGIN=1,
                           FUN=function(b) all(b<=beta.OS.quantiles) ))
    return(pr.joint - (1-ALPHA.joint))
  }
  ALPHA.check.fn.2s <- function(a) {
    beta.OS.quantiles.L <- qbeta( qbeta(  a, 1:J, J:1),  k, n+1-k)
    beta.OS.quantiles.U <- qbeta( qbeta(1-a, 1:J, J:1),  k, n+1-k)
    all2s <- function(b) all(beta.OS.quantiles.L<=b) && all(b<=beta.OS.quantiles.U)
    pr.joint <- mean( apply(beta.J.r.matrix, MARGIN=1, FUN=all2s) )
    return(pr.joint - (1-ALPHA.joint))
  }
  ALPHA.pt.lo1 <- uniroot(ALPHA.check.fn.lo1, interval=c(0, 1))$root
  ALPHA.pt.up1 <- uniroot(ALPHA.check.fn.up1, interval=c(0, 1))$root
  ALPHA.pt.2s  <- uniroot(ALPHA.check.fn.2s,  interval=c(0, 1))$root
  y.pts.lo1 <- qbeta( qbeta(  ALPHA.pt.lo1, 1:J, J:1),  k, n+1-k)
  y.pts.up1 <- qbeta( qbeta(1-ALPHA.pt.up1, 1:J, J:1),  k, n+1-k)
  y.pts.lo2 <- qbeta( qbeta(  ALPHA.pt.2s,  1:J, J:1),  k, n+1-k)
  y.pts.up2 <- qbeta( qbeta(1-ALPHA.pt.2s,  1:J, J:1),  k, n+1-k)
  ret <- cbind(winning.bids=wins, UCB.lower.1s=y.pts.lo1, UCB.upper.1s=y.pts.up1,
               UCB.lower.2s=y.pts.lo2, UCB.upper.2s=y.pts.up2)
  if (!is.null(plot.type) && !is.na(plot.type)) {
    if (plot.type=='1s') {
      plot.bid.CDF.band(band=ret, conf.level=conf.level, two.sided=FALSE)
    } else if (plot.type=='2s') {
      plot.bid.CDF.band(band=ret, conf.level=conf.level, two.sided=TRUE)
    } else if (plot.type!='') warning("plot.type is neither '1s' nor '2s' so is ignored (no plot).")
  }
  invisible(ret)
}

# Plot fn called above; feel free to customize
# The first argument 'band' is the return value from bid.CDF.band()
# conf.level needs to be provided separately (e.g., 0.90 for 90% confidence level)
# two.sided: FALSE for one-sided, TRUE for two-sided (slightly wider but usually pretty similar)
# main: title for plot
# plot.legend: include a legend? (TRUE/FALSE)
# mar: margins argument for par()
# ...: any other arguments to pass along to plot(), like xlim, yaxs, etc.
plot.bid.CDF.band <- function(band, conf.level=NA, two.sided=FALSE, main=NULL, 
                              plot.legend=TRUE, mar=c(5.1,4.1,6.1,2.1), ...) {
  par(mar=mar)
  if (is.null(main)) {
    main <- sprintf('%s uniform confidence band%s\n',
                    ifelse(is.na(conf.level),'',sprintf('%g%%',100*conf.level)),
                    ifelse(two.sided,'','s'))
  }
  wins <- band[,'winning.bids']
  y.pts.lo1 <- band[,'UCB.lower.1s']
  y.pts.lo2 <- band[,'UCB.lower.2s']
  y.pts.up1 <- band[,'UCB.upper.1s']
  y.pts.up2 <- band[,'UCB.upper.2s']
  plot(x=range(wins), y=range(y.pts.up1,y.pts.lo1,y.pts.up2,y.pts.lo2), type='n',
       xlab="Bid value", ylab="Bid CDF", ylim=0:1, main=main, yaxt='n', ...)
  axis(side=2, at=0:4/4)
  bignum <- 100*(max(wins)-min(wins))
  if (two.sided) {
    lines(x=c(wins[1]-bignum,wins[1],wins,max(wins)+bignum), 
          y=c(0,0,y.pts.lo2,max(y.pts.lo2)),
          type='s', lwd=2, lty=2, col=1)
    lines(x=c(wins[1]-bignum,wins,max(wins),max(wins)+bignum),
          y=c(min(y.pts.up2),y.pts.up2,1,1),
          type='S', lwd=2, lty=2, col=1)
    if (plot.legend) legend('bottomright', legend='Two-sided', lty=2, lwd=2)
  } else {
    lines(x=c(wins[1]-bignum,wins[1],wins,max(wins)+bignum),
          y=c(0,0,y.pts.lo1,max(y.pts.lo1)),
          type='s', lwd=2, lty=2, col=1)
    lines(x=c(wins[1]-bignum,wins,max(wins),max(wins)+bignum),
          y=c(min(y.pts.up1),y.pts.up1,1,1),
          type='S', lwd=2, lty=2, col=1)
    if (plot.legend) {
      tmp <- legend('top', plot=FALSE, lty=2, lwd=2, ncol=2,
                    legend=c('Lower 1-sided','Upper 1-sided'))
      newxy <- c(tmp$rect$left, tmp$rect$top + 1.01*tmp$rect$h)
      legend(x=newxy[1], y=newxy[2], legend=c('Lower 1-sided','Upper 1-sided'),
             lty=2, lwd=2, ncol=2, xpd=NA)
    }
  }
}

# # Examples
# set.seed(112358)
# J <- 100; n <- 10; price <- 2
# fake.bids <- matrix(data=rnorm(J*n, mean=100, sd=10), nrow=J)
# fake.wins <- apply(X=fake.bids, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
# ret <- bid.CDF.band(J=J, n=n, price=price, wins=fake.wins, conf.level=0.90)
# # 
# set.seed(112358)
# J <- 10; n <- 4; price <- 2
# fake.bids <- matrix(data=rnorm(J*n, mean=100, sd=10), nrow=J)
# fake.wins <- apply(X=fake.bids, MARGIN=1, FUN=function(x) sort(x=x, decreasing=TRUE)[price])
# ret <- bid.CDF.band(J=J, n=n, price=price, wins=fake.wins, conf.level=0.90)

#EOF

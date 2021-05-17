# Bayesian bootstrap for stochastic dominance inference
#  with continuous distributions
# Code for the paper by David M. Kaplan and Longhao Zhuo
# Questions? Comments? kaplandm@missouri.edu

#
# Main function
# x: numeric vector
# y: CDF function, or else numeric vector
# weights.x, weights.y: sampling weights, if applicable
# max.SD.order: maximum order of stochastic dominance; e.g., 1 for first-order stochastic dominance, 2 for second-order, etc.  All orders from 1 through max.SD.order will be evaluated.
# BREP: number of Bayesian bootstrap replications
# z: only test restricted SD, over domain from -Inf to z (e.g., z=poverty line)
# eval.pts.mult: distribution functions will be evaluated at eval.pts.mult times the sample size number of points
# X.support, Y.support: for two-sample inference, if known, the supports of the X and Y distributions; e.g., if X is known to be non-negative, c(0,Inf), or if Y is known to lie between 0 and 1, c(0,1).  If unknown, leave as c(-Inf,Inf).
# interpolate: FALSE to draw step-function CDFs from the posterior, TRUE to linearly interpolate the CDF between the observed values. (Note that the step CDF SD2 the interpolated CDF!)
# na.rm: must be TRUE; NA values are removed
# Return value: for each order j in 1:max.SD.order, posterior probabilities of x SDj y and y SDj x (data.frame columns xSDy and ySDx)
#
SD.BB <- function(x,y,weights.x=NULL,weights.y=NULL,max.SD.order=1,z=NULL,BREP=99,eval.pts.mult=1,X.support=c(-Inf,Inf),Y.support=c(-Inf,Inf),interpolate=TRUE,na.rm=TRUE) {
  if (!na.rm) stop("na.rm must be TRUE")
  if (missing(x) || missing(y)) stop("arguments x and y cannot be missing")
  max.SD.order <- as.integer(max.SD.order)
  if (max.SD.order<1) stop("max.SD.order must be at least 1")
  #
  # Validate and prep x data
  tmp <- validate.and.prep.data.fn(x,weights.x)
  x <- tmp$data;  weights.x <- tmp$weights
  #
  if (is.numeric(y)) {
    # Validate and prep y data
    tmp <- validate.and.prep.data.fn(y,weights.y)
    y <- tmp$data;  weights.y <- tmp$weights
    return(SD.BB.2s(x=x,y=y,weights.x=weights.x,weights.y=weights.y,max.SD.order=max.SD.order,z=z,BREP=BREP,eval.pts.mult=eval.pts.mult,X.support=X.support,Y.support=Y.support,interpolate=interpolate))
  } else {
    if (is.character(y)) y <- get(y, mode="function", envir=parent.frame())
    if (!is.function(y)) stop("'y' must be numeric or a function or a string naming a valid function")
    return(SD.BB.1s(x=x,y=y,weights.x=weights.x,max.SD.order=max.SD.order,z=z,BREP=BREP,eval.pts.mult=eval.pts.mult,interpolate=interpolate))
  }
}

# Helper function: prep/validate data
validate.and.prep.data.fn <- function(x,wgt) {
  if (!is.numeric(x)) stop("data must be numeric")
  if (is.null(wgt)) wgt <- rep.int(1,length(x))
  if (any(wgt<0)) stop("cannot have negative entries in weights")
  wgt <- wgt[!is.na(x)];  x <- x[!is.na(x)]
  tmpx <- sort(unique(x))
  tmpw <- rep.int(NA,length(tmpx))
  for (i in 1:length(tmpw)) tmpw[i] <- sum(wgt[x==tmpx[i]])
  x <- tmpx;  wgt <- tmpw
  if (sum(wgt)<=0) stop("not enough (non-NA, positive weight) data")
  return(data.frame(data=x,weights=wgt))
}

# Helper function for 1-sample inference
SD.BB.1s <- function(x,y,weights.x,max.SD.order,z,BREP,eval.pts.mult,interpolate) {
  SD.chks <- matrix(NA,nrow=BREP,ncol=max.SD.order)
  wgts <- c(weights.x,mean(weights.x)) #extra for continuity correction, like in Banks (1988)
  n <- length(x)
  eval.pts <- rep.int(NA,(n-1)*eval.pts.mult+1)
  if (eval.pts.mult==1) {
    eval.pts <- x
  } else {
    xdiffs <- x[2:n]-x[1:(n-1)]
    baseseq <- seq(from=1,to=1+(n-2)*eval.pts.mult,by=eval.pts.mult)
    for (i in 0:(eval.pts.mult-1)) {
      eval.pts[baseseq+i] <- x[1:(n-1)]+xdiffs*(i/eval.pts.mult)
    }
    eval.pts[length(eval.pts)] <- x[n]
  }
  if (!is.null(z)) eval.pts <- eval.pts[eval.pts<=z]
  if (length(eval.pts)==0) stop("z is too small; pick a larger value, or don't specify.")
  for (b in 1:BREP) {
    # Draw from Dirichlet
    pxs <- Dirichlet.draw.fn(wgts)
    # Check SD1
    # Make Dirichlet-based CDF at eval.pts
    Fx <- rep.int(NA,(n-1)*eval.pts.mult+1)
    Fx.obs <- cumsum(pxs[1:n])
    Fxdiffs <- Fx.obs[2:n]-Fx.obs[1:(n-1)]
    if (eval.pts.mult==1) {
      Fx <- Fx.obs
    } else if (interpolate) {
      for (i in 0:(eval.pts.mult-1)) {
        Fx[baseseq+i] <- Fx.obs[1:(n-1)]+Fxdiffs*(i/eval.pts.mult)
      }
      Fx[length(Fx)] <- Fx.obs[n]
    } else {
      Fx <- c(rep(Fx.obs[-n],each=eval.pts.mult),Fx.obs[n])
    }
    Fy <- y(eval.pts)
    SD.chks[b,1] <- all(Fx<=Fy)-all(Fx>=Fy)
    # SD.chks[b,j] <- eval(call(sprintf("SD%d.chk",j),quote(Fx),quote(Fy)))
    if (max.SD.order>=2) {
      stop("Not yet implemented: SD2")
      eval.pts.diff <- eval.pts[-1]-eval.pts[-length(eval.pts)]
      if (interpolate) {
        XXX
      } else {
        increments.x <- eval.pts.diff*Fx[-length(Fx)]
        integrals.x <- cumsum(increments.x)
      }
      increments.y <- rep.int(NA,length(increments.x))
      for (i in 1:length(increments.y)) {
        increments.y[i] <- integrate(f=y,lower=eval.pts[i],upper=eval.pts[i+1])
      }
      integrals.y <- cumsum(increments.y)
      SD.chks[b,2] <- all(integrals.x<=integrals.y)-all(integrals.x>=integrals.y)
    }
    if (max.SD.order>=3) {
      stop("Not yet implemented: SDj, j>=3")
    }
  }
  ret <- data.frame(SD.order=1:max.SD.order,xSDy=NA,ySDx=NA)
  ret$xSDy <- colMeans(SD.chks>0,na.rm=TRUE)
  ret$ySDx <- colMeans(SD.chks<0,na.rm=TRUE)
  return(ret)
}

# Helper fn for 2-sample inference
SD.BB.2s <- function(x,y,weights.x,weights.y,max.SD.order,z,BREP,eval.pts.mult,X.support,Y.support,interpolate) {
  # idea: code eval.pts[i] as Fx[ind1[i]]+prop[i]*Fxdiff[ind1[i]]
  SD.chks <- matrix(NA,nrow=BREP,ncol=max.SD.order)
  wgts.x <- c(weights.x,mean(weights.x)) #extra for continuity correction, like in Banks (1988)
  wgts.y <- c(weights.y,mean(weights.y))
  n.x <- length(x); n.y <- length(y)
  # eval.pts: since CDFs are both linear splines, only need to evaluate at kinks; to avoid relying on extrapolation, start at max(x[1],y[1]) and end at min(x[n.x],y[n.y])
  eval.pts <- sort(unique(c(x[x>=y[1] & x<=y[n.y]],y[y>=x[1] & y<=x[n.x]])))
  if (!is.null(z)) eval.pts <- eval.pts[eval.pts<=z]
  if (length(eval.pts)==0) stop("z is too small; pick a larger value, or don't specify.")
  eval.indx <- eval.indy <- eval.propx <- eval.propy <- rep.int(NA,length(eval.pts))
  x <- c(x,max(x)+1);  y <- c(y,max(y)+1)  #for propx, propy below; can't use Inf b/c 0*Inf=NaN
  for (i in 1:length(eval.pts)) {
    eval.indx[i] <- max(which(x<=eval.pts[i]))
    eval.indy[i] <- max(which(y<=eval.pts[i]))
    eval.propx[i] <- (eval.pts[i]-x[eval.indx[i]]) / (x[eval.indx[i]+1]-x[eval.indx[i]])
    eval.propy[i] <- (eval.pts[i]-y[eval.indy[i]]) / (y[eval.indy[i]+1]-y[eval.indy[i]])
  }
  for (b in 1:BREP) {
    # Draw from Dirichlet
    pxs <- Dirichlet.draw.fn(wgts.x)
    pys <- Dirichlet.draw.fn(wgts.y)
    # Fx and Fy (posterior-sampled CDFs evaluated at eval.pts) are linear on each segment [i,i+1], but may have kinks at any Fx[i] or Fy[i]
    Fx <- Fy <- rep.int(NA,length(eval.pts))
    Fx.obs <- cumsum(pxs[1:n.x]); Fy.obs <- cumsum(pys[1:n.y])
    tmpxdiff <- c(Fx.obs[2:n.x]-Fx.obs[1:(n.x-1)], 0)
    tmpydiff <- c(Fy.obs[2:n.y]-Fy.obs[1:(n.y-1)], 0)
    Fx <- Fx.obs[eval.indx]
    Fy <- Fy.obs[eval.indy]
    if (interpolate) {
      Fx <- Fx + eval.propx*tmpxdiff[eval.indx]
      Fy <- Fy + eval.propy*tmpydiff[eval.indy]
    }
    # SD1 check
    SD.chks[b,1] <- all(Fx<=Fy)-all(Fx>=Fy)
    # SD.chks[b,1] <- SD1.chk(Fx,Fy) #eval(call(sprintf("SD%d.chk",j),quote(Fx),quote(Fy)))
    if (max.SD.order>=2) {
      # Integrate: if interpolate, then sideways trapezoids, "height" (width)=eval.pts[i+1]-eval.pts[i], "widths" (heights)=Fx[i]-Fy[i] and Fx[i+1]-Fy[i+1]
      tmp <- length(eval.pts)
      if (interpolate) {
        trapezoids <- (1/2)*(Fx[1:(tmp-1)]-Fy[1:(tmp-1)]+Fx[2:tmp]-Fy[2:tmp])*(eval.pts[2:tmp]-eval.pts[1:(tmp-1)])
      } else { #rectangles; rem: Fx and Fy are cadlag, so use *left* eval pt
        trapezoids <- (Fx[1:(tmp-1)]-Fy[1:(tmp-1)])*(eval.pts[2:tmp]-eval.pts[1:(tmp-1)])
      }
      tmp <- cumsum(trapezoids)
      SD.chks[b,2] <- all(tmp<=0)-all(tmp>=0)
    }
    if (max.SD.order>=3) {
      # Davidson & Duclos (2000), Section 2, eqn (2), etc.
      # for all z, same sign of \int_{-\infty}^{z}(z-t)[G(t)-F(t)]dt
      if (interpolate) { # XXX DOUBLE-CHECK
        Fxdiff <- Fx[-1]-Fx[-length(Fx)]
        Fydiff <- Fy[-1]-Fy[-length(Fy)]
        eval.diff <- eval.pts[-1]-eval.pts[-length(eval.pts)]
        m <- (Fxdiff-Fydiff)/eval.diff
        a <- Fx[-1]-Fy[-1]-m*eval.pts[-1]
        integrals <- rep.int(NA,length(eval.pts)-1)
        for (i in 1:length(integrals)) {
          z <- eval.pts[i+1]
          # [zax+(1/2)(zm-a)x^2-(m/3)x^3]_{min}^{max}
          increments <- z*a[1:i]*eval.pts[2:(i+1)]+(1/2)*(z*m[1:i]-a[1:i])*eval.pts[2:(i+1)]^2 - (m[1:i]/3)*eval.pts[2:(i+1)]^3
          increments <- 
            with(list(tmpa=a[1:i],tmpm=m[1:i],tmpx=eval.pts[2:(i+1)]),
                 z*tmpa*tmpx +(1/2)*(z*tmpm-tmpa)*tmpx^2 -(tmpm/3)*tmpx^3) -
            with(list(tmpa=a[1:i],tmpm=m[1:i],tmpx=eval.pts[1:i]),
                 z*tmpa*tmpx +(1/2)*(z*tmpm-tmpa)*tmpx^2 -(tmpm/3)*tmpx^3)
          integrals[i] <- sum(increments)
        }
        SD.chks[b,3] <- all(integrals<=0) - all(integrals>=0)
      } else { 
        #weighted version of (17) in Davidson & Duclos (2000), ignoring (s-1)!
        s <- 3
        D3diffs <- rep(NA,length(eval.pts)-1)
        for (i in 1:length(D3diffs)) {
          z <- eval.pts[i+1]
          D3x <- sum(pxs[1:eval.indx[i+1]]*(z-x[1:eval.indx[i+1]])^(s-1))
          D3y <- sum(pys[1:eval.indy[i+1]]*(z-y[1:eval.indy[i+1]])^(s-1))
          D3diffs[i] <- D3x-D3y
        }
        SD.chks[b,3] <- all(D3diffs<=0) - all(D3diffs>=0)
      }
    }
    if (max.SD.order>=4) {
      if (interpolate) {
        stop("Not yet implemented: two-sample SDj with j>=4 and interpolate; try setting interpolate=FALSE")
      } else {
        for (s in 4:max.SD.order) {
          Ds.diffs <- rep(NA,length(eval.pts)-1)
          for (i in 1:length(Ds.diffs)) {
            z <- eval.pts[i+1]
            Ds.x <- sum(pxs[1:eval.indx[i+1]]*(z-x[1:eval.indx[i+1]])^(s-1))
            Ds.y <- sum(pys[1:eval.indy[i+1]]*(z-y[1:eval.indy[i+1]])^(s-1))
            Ds.diffs[i] <- Ds.x-Ds.y
          }
          SD.chks[b,s] <- all(Ds.diffs<=0) - all(Ds.diffs>=0)
        }
      }
    }
  }
  ret <- data.frame(SD.order=1:max.SD.order,xSDy=NA,ySDx=NA)
  ret$xSDy <- colMeans(SD.chks>0,na.rm=TRUE)
  ret$ySDx <- colMeans(SD.chks<0,na.rm=TRUE)
  return(ret)
}

# Dirichlet draw
Dirichlet.draw.fn <- function(weights) {
  # Draw Gamma(a_i,1), divide by sum
  if (length(unique(weights))==1) {
    gammas <- rgamma(n=length(weights),shape=weights[1],rate=1)
  } else {
    gammas <- rep.int(NA,length(weights))
    for (i in 1:length(gammas)) gammas[i] <- rgamma(n=1,shape=weights[i],rate=1)
  }
  return(gammas/sum(gammas))
}


# Simulations: Bayesian and frequentist inference is very different for stochastic dominance
# one-sample
SD.BB.1s.sims.fn <- function(NREP=200,n=50,ALPHA=0.1,CDF=pnorm,DGP=rnorm) {
  set.seed(112358)
  rejs <- data.frame(rej.BB.xSD1y=rep.int(NA,NREP),
                     rej.BB.ySD1x=NA,
                     rej.KS.xSD1y=NA,rej.KS.ySD1x=NA)
  for (irep in 1:NREP) {
    x <- DGP(n)
    ret <- SD.BB(x,CDF)
    rejs$rej.BB.xSD1y[irep] <- (ret[1,2]<ALPHA)
    rejs$rej.BB.ySD1x[irep] <- (ret[1,3]<ALPHA)
    rejs$rej.KS.xSD1y[irep] <- (ks.test(x,CDF,alternative="greater")$p.value<ALPHA)
    rejs$rej.KS.ySD1x[irep] <- (ks.test(x,CDF,alternative="less")$p.value<ALPHA)
  }
  colMeans(rejs)
}
# two-sample
SD.BB.2s.sims.fn <- function(NREP=200,nx=50,ny=50,ALPHA=0.1,DGPx=rnorm,DGPy=rnorm) {
  set.seed(112358)
  rejs <- data.frame(rej.BB.xSD1y=rep.int(NA,NREP),
                     rej.BB.ySD1x=NA,
                     rej.KS.xSD1y=NA,rej.KS.ySD1x=NA,
                     Pr.xSD1y=NA,Pr.ySD1x=NA)
  for (irep in 1:NREP) {
    x <- DGPx(nx); y <- DGPy(ny)
    ret <- SD.BB(x,y)
    rejs$Pr.xSD1y[irep] <- ret[1,2]
    rejs$Pr.ySD1x[irep] <- ret[1,3]
    rejs$rej.BB.xSD1y[irep] <- (ret[1,2]<ALPHA)
    rejs$rej.BB.ySD1x[irep] <- (ret[1,3]<ALPHA)
    rejs$rej.KS.xSD1y[irep] <- (ks.test(x,y,alternative="greater")$p.value<ALPHA)
    rejs$rej.KS.ySD1x[irep] <- (ks.test(x,y,alternative="less")$p.value<ALPHA)
  }
  colMeans(rejs)
}
# SD.BB.2s.sims.fn()
# SD.BB.2s.sims.fn(DGPy=function(n)rnorm(n,0.5,1))
# SD.BB.2s.sims.fn(DGPy=function(n)rnorm(n,1,1)) #obvious to both

SD.BB.2s.sims.SDj.fn <- function(NREP=200,nx=50,ny=50,ALPHA=0.1,DGPx=rnorm,DGPy=rnorm,max.j=2) {
  set.seed(112358)
  j.by.type.by.rep <- array(NA,dim=c(max.j,4,NREP))
  rejs <- data.frame(rej.BB.xSD2y=rep.int(NA,NREP),
                     rej.BB.ySD2x=NA,
                     rej.BB.xSD1y=NA,rej.BB.ySD1x=NA,
                     Pr.xSD1y=NA,Pr.ySD1x=NA,
                     Pr.xSD2y=NA,Pr.ySD2x=NA)
  for (irep in 1:NREP) {
    x <- DGPx(nx); y <- DGPy(ny)
    ret <- SD.BB(x,y,max.SD.order=2)
    rejs$Pr.xSD1y[irep] <- ret[1,2]
    rejs$Pr.ySD1x[irep] <- ret[1,3]
    rejs$Pr.xSD2y[irep] <- ret[2,2]
    rejs$Pr.ySD2x[irep] <- ret[2,3]
    rejs$rej.BB.xSD2y[irep] <- (ret[2,2]<ALPHA)
    rejs$rej.BB.ySD2x[irep] <- (ret[2,3]<ALPHA)
    rejs$rej.BB.xSD1y[irep] <- (ret[1,2]<ALPHA)
    rejs$rej.BB.ySD1x[irep] <- (ret[1,3]<ALPHA)
    j.by.type.by.rep[1:max.j,1:4,irep] <- as.matrix(cbind(1*(ret[,2:3]<ALPHA),ret[,2:3]))
  }
  colMeans(rejs)
  j.by.type.avg <- rowMeans(j.by.type.by.rep,dims=2)
  data.frame(j=1:max.j,rej.xSDjy=j.by.type.avg[,1],rej.ySDjx=j.by.type.avg[,2],Pr.xSDjy=j.by.type.avg[,3],Pr.ySDjx=j.by.type.avg[,4])
}
# SD.BB.2s.sims.SDj.fn()
# SD.BB.2s.sims.SDj.fn(nx=500,ny=500) #~50% rej SD2 each way
# SD.BB.2s.sims.SDj.fn(DGPy=function(n)rnorm(n,0.5,1)) #ySD1x
# SD.BB.2s.sims.SDj.fn(DGPy=function(n)rnorm(n,0,0.5)) #ySD2x

#EOF
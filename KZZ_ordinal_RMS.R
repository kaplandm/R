# Implementation of Andrews and Barwick (2012) RMS test
# for ordinal SD1
# Original by Longhao Zhuo; revision by kaplandm@missouri.edu

library(quadprog)

# Test of SD1 (s=rep(1,length(Xc))) or SC (s=c(rep(1,j),rep(-1,length(Xc)-j)))
# Input: Xc is vector of counts for each category value in the X sample; similarly for Y
# ALPHA is the level of the test, like 0.05 for 5%
# RMSREP is the number of bootstrap replications
RMS.ordinal.test <- function(Xc, Yc, s, ALPHA, RMSREP=400) {
  if (missing(s)) stop("Must provide s; for H0:X SD1 Y, s=rep(1,length(Xcts))")
  if (missing(ALPHA)) stop("Must provide ALPHA, like ALPHA=0.1 for level 10% test.")
  if (missing(Xc) || missing(Yc)) stop("Must provide Xc and Yc (counts of each category).")
  if (length(Yc)!=length(Xc)) stop("Must have length(Xc)==length(Yc).")
  nX <- sum(Xc);  Fx <- cumsum(Xc)/nX;  Fx <- Fx[-length(Fx)]
  nY <- sum(Yc);  Fy <- cumsum(Yc)/nY;  Fy <- Fy[-length(Fy)]
  s <- s[1:length(Fx)]
  mbar <- (Fy-Fx)*s
  # \hat\Sigma defined just before (2.3)
  Sigmahat <- matrix(data=NA, nrow=length(Fx), ncol=length(Fx))
  for (i in 1:nrow(Sigmahat)) {
    for (j in 1:i) {
      # Assuming X and Y independent, and iid samples, and scale by nX
      Sigmahat[i,j] <- Sigmahat[j,i] <-
        (Fy[j]*(1-Fy[i])*nX/nY + Fx[j]*(1-Fx[i]))*s[i]*s[j]
    }
  }
  # Apply equation (2.6) from Andrews and Barwick (2012)
  Dhat <- diag(diag(Sigmahat))
  Dhatsqrtinv <- diag(1/sqrt(diag(Sigmahat)))
  Omegahat <- Dhatsqrtinv %*% Sigmahat %*% Dhatsqrtinv
  Sigmatilde <- Sigmahat  +  max(0, 0.012 - det(Omegahat)) * Dhat
  # 
  deltahat <- min(Omegahat)  # (2.9), (2.10)
  tmp <- kappa.eta.fn(delta=deltahat, P=length(Fx))  # (2.10)
  kappahat <- tmp['KAPPA']
  etahat <- tmp['ETA']
  phat <- which(sqrt(nX)*mbar/sqrt(diag(Sigmahat)) <= kappahat)  # (2.7)
  if (length(phat)==0) return(FALSE)  #phat <- 1
  Sigmatilde.sel.inv <- solve(Sigmatilde[phat,phat])
  targmin <- solve.QP(Dmat=Sigmatilde.sel.inv,
                      dvec=t(nX^(1/2)*mbar[phat]) %*% Sigmatilde.sel.inv,
                      Amat=diag(length(phat)), bvec=rep(0,length(phat)),
                      factorized=FALSE)$solution
  AQLRhat <- c(t(nX^(1/2)*mbar[phat]-targmin) %*%
               Sigmatilde.sel.inv %*% (nX^(1/2)*mbar[phat]-targmin))
  # 
  # Bootstrap world (basically repeat same as above)
  AQLRstars <- rep(NA,RMSREP)
  Xcstars <- t(rmultinom(n=RMSREP, size=nX, prob=Xc/sum(Xc)))
  Ycstars <- t(rmultinom(n=RMSREP, size=nY, prob=Yc/sum(Yc)))
  for (r in 1:RMSREP) {
    Xcstar <- Xcstars[r,];  Ycstar <- Ycstars[r,]
    Fxstar <- cumsum(Xcstar)/sum(Xcstar)
    Fxstar <- Fxstar[-length(Fxstar)]
    Fystar <- cumsum(Ycstar)/sum(Ycstar)
    Fystar <- Fystar[-length(Fystar)]
    mbarstar <- (Fystar-Fxstar)*s - mbar
    Sigmahatstar <- matrix(data=NA, nrow=length(Fxstar), ncol=length(Fxstar))
    for (i in 1:nrow(Sigmahatstar)) {
      for (j in 1:i) {
        Sigmahatstar[i,j] <- Sigmahatstar[j,i] <- 
          (Fystar[j]*(1-Fystar[i])*nX/nY + Fxstar[j]*(1-Fxstar[i]))*s[i]*s[j]
      }
    }
    Dhatstar <- diag(diag(Sigmahatstar))
    Dhatsqrtinvstar <- diag(1/sqrt(diag(Sigmahatstar)))
    Omegahatstar <- Dhatsqrtinvstar %*% Sigmahatstar %*% Dhatsqrtinvstar
    Sigmatildestar <- Sigmahatstar  +  max(0, 0.012 - det(Omegahatstar)) * Dhatstar
    Sigmatilde.sel.invstar <- solve(Sigmatildestar[phat,phat])
    if (any(is.na(Sigmatilde.sel.invstar) | is.nan(Sigmatilde.sel.invstar))) {
      AQLRstars[r] <- Inf
      next
    }
    targminstar <- 
      solve.QP(Dmat=Sigmatilde.sel.invstar,
               dvec=t(nX^(1/2)*mbarstar[phat]) %*% Sigmatilde.sel.invstar,
               Amat=diag(length(phat)), bvec=rep(0,length(phat)),
               factorized=FALSE)$solution
    AQLRstars[r] <- t(nX^(1/2)*mbarstar[phat]-targminstar) %*%
      Sigmatilde.sel.invstar %*% (nX^(1/2)*mbarstar[phat]-targminstar)
  }
  cv <- quantile(x=AQLRstars, probs=1-ALPHA, type=6) + etahat
  return(c(AQLRhat > cv))
}


# Table I in their paper (p. 2811)
kappa.eta.fn <-  function(delta, P) {
  if (delta < -.975){
    KAPPA <- 2.9; ETAAUTO <- .025;
  }else if ( delta >= -.975 && delta < -.95){
    KAPPA <- 2.9; ETAAUTO <- .026;
  }else if ( delta >= -.95 && delta < -.90){
    KAPPA <- 2.9; ETAAUTO <- .021;
  }else if ( delta >= -.90 && delta < -.85){
    KAPPA <- 2.8; ETAAUTO <- .027;
  }else if ( delta >= -.85 && delta < -.80){
    KAPPA <- 2.7; ETAAUTO <- .062;
  }else if ( delta >= -.80 && delta < -.75){
    KAPPA <- 2.6; ETAAUTO <- .104;    
  }else if ( delta >= -.75 && delta < -.70){
    KAPPA <- 2.6; ETAAUTO <- .103;
  }else if ( delta >= -.70 && delta < -.65){
    KAPPA <- 2.5; ETAAUTO <- .131;
  }else if ( delta >= -.65 && delta < -.60){
    KAPPA <- 2.5; ETAAUTO <- .122;
  }else if ( delta >= -.60 && delta < -.55){
    KAPPA <- 2.5; ETAAUTO <- .113;
  }else if ( delta >= -.55 && delta < -.50){
    KAPPA <- 2.5; ETAAUTO <- .104;
  }else if ( delta >= -.50 && delta < -.45){
    KAPPA <- 2.4; ETAAUTO <- .124;
  }else if ( delta >= -.45 && delta < -.40){
    KAPPA <- 2.2; ETAAUTO <- .158;
  }else if ( delta >= -.40 && delta < -.35){
    KAPPA <- 2.2; ETAAUTO <- .133;
  }else if ( delta >= -.35 && delta < -.30){
    KAPPA <- 2.1; ETAAUTO <- .138;
  }else if ( delta >= -.30 && delta < -.25){
    KAPPA <- 2.1; ETAAUTO <- .111;
  }else if ( delta >= -.25 && delta < -.20){
    KAPPA <- 2.1; ETAAUTO <- .082;
  }else if ( delta >= -.20 && delta < -.15){
    KAPPA <- 2.0; ETAAUTO <- .083;
  }else if ( delta >= -.15 && delta < -.10){
    KAPPA <- 2.0; ETAAUTO <- .074;
  }else if ( delta >= -.10 && delta < -.05){
    KAPPA <- 1.9; ETAAUTO <- .082;
  }else if ( delta >= -.05 && delta < 0.00){
    KAPPA <- 1.8; ETAAUTO <- .075;
  }else if ( delta >= 0.00 && delta < 0.05){
    KAPPA <- 1.5; ETAAUTO <- .114;
  }else if ( delta >= 0.05 && delta < 0.10){
    KAPPA <- 1.4; ETAAUTO <- .112;
  }else if ( delta >= 0.10 && delta < 0.15){
    KAPPA <- 1.4; ETAAUTO <- .083;
  }else if ( delta >= 0.15 && delta < 0.20){
    KAPPA <- 1.3; ETAAUTO <- .089;
  }else if ( delta >= 0.20 && delta < 0.25){
    KAPPA <- 1.3; ETAAUTO <- .058;
  }else if ( delta >= 0.25 && delta < 0.30){
    KAPPA <- 1.2; ETAAUTO <- .055;
  }else if ( delta >= 0.30 && delta < 0.35){
    KAPPA <- 1.1; ETAAUTO <- .044;
  }else if ( delta >= 0.35 && delta < 0.40){
    KAPPA <- 1.0; ETAAUTO <- .040;
  }else if ( delta >= 0.40 && delta < 0.45){
    KAPPA <- 0.8; ETAAUTO <- .051;
  }else if ( delta >= 0.45 && delta < 0.50){
    KAPPA <- 0.8; ETAAUTO <- .023;
  }else if ( delta >= 0.50 && delta < 0.55){
    KAPPA <- 0.6; ETAAUTO <- .033;
  }else if ( delta >= 0.55 && delta < 0.60){
    KAPPA <- 0.6; ETAAUTO <- .013;
  }else if ( delta >= 0.60 && delta < 0.65){
    KAPPA <- 0.4; ETAAUTO <- .016;
  }else if ( delta >= 0.65 && delta < 0.70){
    KAPPA <- 0.4; ETAAUTO <- .000;
  }else if ( delta >= 0.70 && delta < 0.75){
    KAPPA <- 0.2; ETAAUTO <- .003;
  }else if ( delta >= 0.75 && delta < 0.80){
    KAPPA <- 0.0; ETAAUTO <- .002;
  }else if ( delta >= 0.80 && delta < 0.85){
    KAPPA <- 0.0; ETAAUTO <- .000;
  }else if ( delta >= 0.85 && delta < 0.90){
    KAPPA <- 0.0; ETAAUTO <- .000;
  }else if ( delta >= 0.90 && delta < 0.95){
    KAPPA <- 0.0; ETAAUTO <- .000;
  }else if ( delta >= 0.95 && delta < 0.975){
    KAPPA <- 0.0; ETAAUTO <- .000;
  }else if ( delta >= 0.975 && delta < 0.99){
    KAPPA <- 0.0; ETAAUTO <- .000;
  }else if ( delta >= 0.99){       
    KAPPA <- 0.0; ETAAUTO <- .000;
  }
  
  # per (2.10) in their paper: THE AUTOETA VALUES ABOVE HAVE BEEN 
  #   DETERMINED FOR P==2 (for which eta2=0). 
  #   FOR P >= 3, WE HAVE TO INCREASE THEM ALL BY THE AMOUNT INDICATED,
  #   i.e., add eta2 (whose values are also in Table I)
  if (P==2) {
    # already fine
  }else if ( P == 3){
    ETAAUTO <- ETAAUTO + .15;
  }else if ( P == 4){
    ETAAUTO <- ETAAUTO + .17;
  }else if ( P == 5){
    ETAAUTO <- ETAAUTO + .24;
  }else if ( P == 6){
    ETAAUTO <- ETAAUTO + .31;
  }else if ( P == 7){
    ETAAUTO <- ETAAUTO + .33;
  }else if ( P == 8){
    ETAAUTO <- ETAAUTO + .37;
  }else if ( P == 9){
    ETAAUTO <- ETAAUTO + .45;
  }else if ( P == 10){
    ETAAUTO <- ETAAUTO + .50
  } else stop("P>10 not supported")
  return(c(KAPPA=KAPPA,ETA=ETAAUTO))
}

# EOF

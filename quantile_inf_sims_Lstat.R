# Simulations for Goldman and Kaplan (2017), "Fractional order statistic approximation for nonparametric conditional quantile inference"
# Questions? Comments? kaplandm@missouri.edu

prob.loaded <- exists("quantile.inf")
success <-
  tryCatch({source('https://raw.githubusercontent.com/kaplandm/R/main/quantile_inf.R');TRUE},
           error=function(w) FALSE)
if (!success) {
  success <- tryCatch({source("quantile_inf.R");TRUE}, error=function(w) FALSE)
  if (!success) {
    if (prob.loaded) {
      warning("Couldn't load quantile_inf.R, but it seems like you already did.")
    } else {
      stop("Failed to source() quantile_inf.R from web or local file.  You may download and source() it yourself, or at least make sure it's in your getwd().  Currently available at https://github.com/kaplandm/R")
    }
  }
}

SAVEFLAG <- TRUE  # TRUE=save results to file; FALSE=display in console only.

if (SAVEFLAG) {
  OUTFILE <- paste0(gsub("-","_",Sys.Date()),"_quantile_inf_sims_Lstat.txt")
  NREPLIC <- 10000
  BREP1 <- 99
  BREP2 <- 100 #0=not studentized.  1=kernel PDF studentized.  >1=BS-var studentized.
  L.ONLY.FLAG <- FALSE
  # For calib: NREPLIC <- 1e5; BREP1 <- BREP2 <- 0; L.ONLY.FLAG <- TRUE
} else {
  OUTFILE <- ""
  NREPLIC <- 1e5
  BREP1 <- 0
  BREP2 <- 0
  L.ONLY.FLAG <- TRUE
}

PARALLEL <- 4 #number of CPUs; if >1, need to "run" this file (quantile_inf_sims.R), not "source"
if (PARALLEL>1) {
  if (!require("parallel") || !require("foreach")) {warning("Install package foreach in order to run PARALLEL."); PARALLEL <- 0}
  if (!require("doParallel")) {warning("Install package doParallel in order to run PARALLEL."); PARALLEL <- 0}
  PARALLEL <- tryCatch({workers <- makeCluster(PARALLEL); registerDoParallel(workers); on.exit(stopCluster(workers),add=TRUE); PARALLEL}, error=function(Err){warning(sprintf("Error creating %d clusters",PARALLEL)); 0})
}

STARTTIME <- Sys.time()
cat(sprintf("NREPLIC=%d,BREP1=%d,BREP2=%d; start time is %s",
            NREPLIC,BREP1,BREP2, format(STARTTIME, "%X, %A, %d %b %Y")),
    file=OUTFILE,append=TRUE,sep='\n')

n <- 25;  p <- 0.5;  ALPHA <- 0.05
CASE1a <- list(p=p,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
CASE1b <- list(p=p,ALPHA=ALPHA,n=n,rfn=function(n)rt(n,1),qfn=function(p)qt(p,1),Fstr="Cauchy")
CASE1c <- list(p=p,ALPHA=ALPHA,n=n,rfn=runif,qfn=qunif,Fstr="Uniform")
CASE1d <- list(p=p,ALPHA=ALPHA,n=n,rfn=rexp,qfn=qexp,Fstr="Exponential")
CASE1e <- list(p=p,ALPHA=ALPHA,n=n,rfn=rlnorm,qfn=qlnorm,Fstr="Log-normal")
n <- 9
CASE1f <- list(p=p,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
CASE1g <- list(p=p,ALPHA=ALPHA,n=n,rfn=function(n)rt(n,1),qfn=function(p)qt(p,1),Fstr="Cauchy")
CASE1h <- list(p=p,ALPHA=ALPHA,n=n,rfn=runif,qfn=qunif,Fstr="Uniform")
CASE1i <- list(p=p,ALPHA=ALPHA,n=n,rfn=rexp,qfn=qexp,Fstr="Exponential")
CASE1j <- list(p=p,ALPHA=ALPHA,n=n,rfn=rlnorm,qfn=qlnorm,Fstr="Log-normal")
n <- 25;  p <- 0.2
CASE1k <- list(p=p,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
CASE1l <- list(p=p,ALPHA=ALPHA,n=n,rfn=function(n)rt(n,1),qfn=function(p)qt(p,1),Fstr="Cauchy")
CASE1m <- list(p=p,ALPHA=ALPHA,n=n,rfn=runif,qfn=qunif,Fstr="Uniform")
CASE1n <- list(p=p,ALPHA=ALPHA,n=n,rfn=rexp,qfn=qexp,Fstr="Exponential")
CASE1o <- list(p=p,ALPHA=ALPHA,n=n,rfn=rlnorm,qfn=qlnorm,Fstr="Log-normal")
n <- 99;  p <- 0.05
CASE1p <- list(p=p,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
n <- 99;  p <- 0.95
CASE1q <- list(p=p,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
n <- 99;  p <- 0.037
CASE1r <- list(p=p,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
CASE1s <- list(p=p,ALPHA=ALPHA,n=n,rfn=function(n)rt(n,1),qfn=function(p)qt(p,1),Fstr="Cauchy")
CASE1t <- list(p=p,ALPHA=ALPHA,n=n,rfn=runif,qfn=qunif,Fstr="Uniform")
CASE1u <- list(p=p,ALPHA=ALPHA,n=n,rfn=rexp,qfn=qexp,Fstr="Exponential")
CASE1v <- list(p=p,ALPHA=ALPHA,n=n,rfn=rlnorm,qfn=qlnorm,Fstr="Log-normal")
n <- 10
CASE1w <- list(p=0.35,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
CASE1x <- list(p=0.40,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
CASE1y <- list(p=0.45,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")
CASE1z <- list(p=0.50,ALPHA=ALPHA,n=n,rfn=rnorm,qfn=qnorm,Fstr="Normal")

if (SAVEFLAG) {
  CASES <- list(CASE1w,CASE1x,CASE1y,CASE1z) #Calib; make sure to change NREPLIC/BREP1/BREP2 above
  CASES <- list(CASE1k,CASE1l,CASE1m,CASE1n,CASE1o,
                CASE1p,CASE1q,
                CASE1r,CASE1s,CASE1t,CASE1u,CASE1v,
                CASE1f,CASE1g,CASE1h,CASE1i,CASE1j,
                CASE1a,CASE1b,CASE1c,CASE1d,CASE1e)
} else {
  CASES <- list(CASE1y,CASE1z) #CASE1w,CASE1x,
}

for (icase in 1:length(CASES)) {
  set.seed(112358)
  CASE <- CASES[[icase]]
  p <- CASE$p;  ALPHA <- CASE$ALPHA;  n <- CASE$n;  Fstr <- CASE$Fstr
  rfn <- CASE$rfn;  Q0 <- CASE$qfn(p)
  Xs <- matrix(rfn(NREPLIC*n),nrow=NREPLIC)
  
  tmp <- matrix(NA,NREPLIC,1)
  methnames <- c("L-stat","Calib","BH","Norm","K15","BStsym","BSt","BS","BSsym")
  # methnames <- c("L-stat","Norm","K15","BStsym","BSt","BS","BSsym")
  CIs <- data.frame(dummy=tmp)
  for (nm in methnames) {
    CIs[[paste0(nm,".lo")]] <- tmp
    CIs[[paste0(nm,".hi")]] <- tmp
  }
  CIs <- CIs[,-1]
  
  # Beran and Hall (1993), page 647; note their \alpha is our 1-ALPHA
  # Also note typo: their r in Q_U and Q_L should be (r-1)
  BHt <- floor((n+1)*p);  BHalpha <- 1-ALPHA
  if (abs(BHt-(n+1)*p)<=0.00001) EXTRA <- 0 else EXTRA <- 1
  found <- FALSE
  for (BHr in 1:min(BHt-1,n-1-BHt-EXTRA)) {
    BHrho <- pbinom(BHt+BHr-1+EXTRA,n,p) - pbinom(BHt-BHr-1,n,p)
    if (BHrho>=BHalpha) {
      BHrho1 <- pbinom(BHt+(BHr-1)-1+EXTRA,n,p) - pbinom(BHt-(BHr-1)-1,n,p)
      BHpi <- (BHalpha-BHrho1) / (BHrho-BHrho1)
      #       CIs[["BH.lo"]][irep] <- BHpi*X[BHt-BHr] + (1-BHpi)*X[BHt-BHr+1]
      #       CIs[["BH.hi"]][irep] <- (1-BHpi)*X[BHt+BHr-1+EXTRA] + BHpi*X[BHt+BHr+EXTRA]
      found <- TRUE
      break
    }
  }
  if (!found) BHt <- BHr <- BHpi <- NA
  
  # L-stat/Hutson, uncalibrated
  Lul <- quantile.inf.CIul(p,n,ALPHA/2)
  Luh <- quantile.inf.CIuh(p,n,ALPHA/2)
  #
  # Calibrated L-stat/Huston
  z <- qnorm(1-ALPHA/2)
  epsl <- (n+1)*Lul - floor((n+1)*Lul)
  a.calib <- ALPHA/2 + (epsl*(1-epsl)*z*exp(-z^2/2))/(sqrt(2*pi)*p*(1-p)*n)
  Lul.c <- quantile.inf.CIul(p,n,a.calib)
  epsh <- (n+1)*Luh - floor((n+1)*Luh)
  a.calib <- ALPHA/2 + (epsh*(1-epsh)*z*exp(-z^2/2))/(sqrt(2*pi)*p*(1-p)*n)
  Luh.c <- quantile.inf.CIuh(p,n,a.calib)
  
  LOOP.STARTTIME <- Sys.time()
  for (irep in 1:NREPLIC) {
    cat(sprintf("rep %d of %d (%s since loop start at %s)\r",irep,NREPLIC,capture.output(print(Sys.time()-LOOP.STARTTIME)),LOOP.STARTTIME), file="")
    flush.console() #otherwise, sometimes won't print until loop is done...
    
    X <- sort(Xs[irep,]) #need to sort for quantile.inf.interp() calls below
    
    # L-stat with numerical integration
    # tmp <- quantile.inf(X=X,p=p,ALPHA=ALPHA,ONESIDED=0,METHOD.TYPE='single')
    CIs[["L-stat.lo"]][irep] <- quantile.inf.interp(X,Lul) #tmp$CI$lo
    CIs[["L-stat.hi"]][irep] <- quantile.inf.interp(X,Luh) #tmp$CI$hi
    
    # L-stat with calibration
    CIs[["Calib.lo"]][irep] <- quantile.inf.interp(X,Lul.c)
    CIs[["Calib.hi"]][irep] <- quantile.inf.interp(X,Luh.c)
    
    if (!L.ONLY.FLAG) {
      # Beran and Hall (1993), page 647; note their \alpha is our 1-ALPHA
      # Also note typo: their r in Q_U and Q_L should be (r-1)
      if (!is.na(BHr)) {
        CIs[["BH.lo"]][irep] <- BHpi*X[BHt-BHr] + (1-BHpi)*X[BHt-BHr+1]
        CIs[["BH.hi"]][irep] <- (1-BHpi)*X[BHt+BHr-1+EXTRA] + BHpi*X[BHt+BHr+EXTRA]
      }
      
      # Kaplan (2015) Edgeworth-based approach
      tmp <- quantile.inf.single.K(X=X,p=p,ALPHA=ALPHA,ONESIDED=0)
      CIs[["K15.lo"]][irep] <- tmp$lo
      CIs[["K15.hi"]][irep] <- tmp$hi
      
      # Bootstrap
      if (BREP1>1) { #else, leave as NA
        Qhat <- quantile.inf.interp(X,p)
        Tstarfn <- function(X,p,n,BREP2) {
          Xstar <- sort(X[sample(1:n,n,replace=TRUE)])
          Qhatstar <- quantile.inf.interp(Xstar,p)
          if (BREP2==0) {
            sdstar2 <- 1
          } else if (BREP2==1) {
            fXhatstar <- density(Xstar,n=1,from=Qhatstar,to=Qhatstar)$y
            sdstar2 <- sqrt(p*(1-p)/(n*fXhatstar^2))
          } else {
            Qstars2 <- rep.int(NA,BREP2)
            for (iBS2 in 1:BREP2) {
              Xstar2 <- sort(Xstar[sample(1:n,n,replace=TRUE)])
              Qstars2[iBS2] <- quantile.inf.interp(Xstar2,p)
            }
            sdstar2 <- sd(Qstars2)
          }
          return(c((Qhatstar-Qhat)/sdstar2, Qhatstar-Qhat))
        }
        Tstars <- matrix(NA,BREP1,2)
        if (PARALLEL>1) {
          clusterSetRNGStream(workers)
          Tstars <- foreach(ipar=1:BREP1,.combine=rbind,.inorder=FALSE) %dopar% {
            Tstarfn(X,p,n,BREP2)
          }
        } else {
          for (iBS in 1:BREP1) {
            Tstars[iBS,] <- Tstarfn(X,p,n,BREP2)
          }
        }
        BStcvsym <- quantile(abs(Tstars[,1]),1-ALPHA,na.rm=TRUE,type=6)
        BStcvs <- quantile(Tstars[,1],c(ALPHA/2,1-ALPHA/2),na.rm=TRUE,type=6)
        BScvsym <- quantile(abs(Tstars[,2]),1-ALPHA,na.rm=TRUE,type=6)
        BScvs <- quantile(Tstars[,2],c(ALPHA/2,1-ALPHA/2),na.rm=TRUE,type=6)
        if (BREP2==0) {
          sdstar <- 1
        } else if (BREP2==1) {
          fXhat <- density(X,n=1,from=Qhat,to=Qhat)$y
          sdstar <- sqrt(p*(1-p)/(n*fXhat^2))
        } else {
          Qstars <- rep.int(NA,BREP2)
          for (iBS2 in 1:BREP2) {
            Xstar <- sort(X[sample(1:n,n,replace=TRUE)])
            Qstars[iBS2] <- quantile.inf.interp(Xstar,p)
          }
          sdstar <- sd(Qstars)
        }
        CIs[["BSt.lo"]][irep] <- Qhat - BStcvs[2]*sdstar
        CIs[["BSt.hi"]][irep] <- Qhat - BStcvs[1]*sdstar
        CIs[["BStsym.lo"]][irep] <- Qhat - BStcvsym*sdstar
        CIs[["BStsym.hi"]][irep] <- Qhat + BStcvsym*sdstar
        CIs[["BS.lo"]][irep] <- Qhat - BScvs[2]
        CIs[["BS.hi"]][irep] <- Qhat - BScvs[1]
        CIs[["BSsym.lo"]][irep] <- Qhat - BScvsym
        CIs[["BSsym.hi"]][irep] <- Qhat + BScvsym
      }
      
      # Normal approx
      Qhat <- quantile.inf.interp(X,p)
      fXhat <- density(X,n=1,from=Qhat,to=Qhat)$y
      SE <- sqrt(p*(1-p)/(n*fXhat^2))
      CIs[["Norm.lo"]][irep] <- Qhat - qnorm(1-ALPHA/2)*SE
      CIs[["Norm.hi"]][irep] <- Qhat - qnorm(ALPHA/2)*SE
    }
  }
  
  # Compute CP, median/mean length
  CPs <- TLs <- THs <- medlens <- meanlens <- data.frame(dummy=NA)
  for (nm in methnames) {
    CPs[[nm]] <- mean(CIs[[paste0(nm,".lo")]]<=Q0 & Q0<=CIs[[paste0(nm,".hi")]],na.rm=TRUE)
    TLs[[nm]] <- mean(CIs[[paste0(nm,".hi")]]<Q0,na.rm=TRUE)
    THs[[nm]] <- mean(CIs[[paste0(nm,".lo")]]>Q0,na.rm=TRUE)
    medlens[[nm]] <- median(CIs[[paste0(nm,".hi")]] - CIs[[paste0(nm,".lo")]],na.rm=TRUE)
    meanlens[[nm]] <- mean(CIs[[paste0(nm,".hi")]] - CIs[[paste0(nm,".lo")]],na.rm=TRUE)
  }
  CPs <- CPs[,-1];  medlens <- medlens[,-1];  meanlens <- meanlens[,-1]
  cat(sprintf("\np=%g,ALPHA=%g,n=%d,Q0=%g,rfn=%s\n",
              p,ALPHA,n,Q0,as.character(as.expression(body(rfn)))),
      file="",append=TRUE,sep='')
  maxnchar <- max(nchar(methnames))

  tmpstr <- sprintf("$%3d$ & $%-5g$ & %15s & %%-%ds & %%5.3f & %%5.3f & %%5.3f & %%4.2f \\\\\n",
                    n,p,Fstr,maxnchar)
  cat(" $n$  &   $p$   &       $F$       & method &   CP  &   TL  &   TH  & medlen\\\\\n",
      file="",append=T,sep='')
  cat("\\hline\n",file=OUTFILE,append=T,sep='')
  for (nm in methnames[1:(length(methnames)-3)]) {
    cat(sprintf(tmpstr, nm, CPs[[nm]], TLs[[nm]], THs[[nm]], medlens[[nm]]),
        file=OUTFILE,append=TRUE,sep='')
  }
  for (nm in methnames) {
    cat(sprintf(tmpstr, nm, CPs[[nm]], TLs[[nm]], THs[[nm]], medlens[[nm]]),
        file="",append=TRUE,sep='')
  }
}
tmpt <- as.numeric(Sys.time()-STARTTIME,units="secs")
tmps <- sprintf("Total time elapsed=%d seconds (i.e., %dmin; i.e., %4.2fhrs)\n\n\n",as.integer(tmpt),as.integer(tmpt/60),tmpt/3600)
cat(tmps,file=OUTFILE,sep="",append=TRUE)

if (PARALLEL>1) stopCluster(workers)

#EOF
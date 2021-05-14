# Feedback: KaplanDM@Missouri.edu
# Simulations for Kaplan (2014) "Nonparametric inference on quantile marginal effects"

#################################################
# SIMULATIONS FOR PAPER

source("qme_inf.R") # from https://github.com/kaplandm/R
# require(boot)
STARTTIME <- Sys.time()
OUTFILE <- "" #"2014_08_19.txt" ######## SET THIS TO SOME VALUE LIKE "out.txt" IN ORDER TO SAVE #####
SAVE.FLAG <- ifelse(OUTFILE=="",FALSE,TRUE);  TRUEFUNOUT.FLAG <- FALSE
tmp <- sprintf("Start time is %s",format(STARTTIME, "%X, %A, %d %b %Y"))
if (OUTFILE=="") {
  cat(tmp,sep='\n')
} else {
  cat(tmp,file=OUTFILE,sep="\n",append=TRUE)
}

set.seed(112358) #for replicability
ONESIDED <- 0 #0:two-sided; 1/-1:upper/lower one-sided
#
TRIM0 <- 0.05;  NUM0 <- 11
CASE1a <- data.frame(p=0.50,n=400,Findex=010,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE1b <- data.frame(p=0.50,n=400,Findex=012,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE1c <- data.frame(p=0.50,n=400,Findex=014,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE1d <- data.frame(p=0.50,n=400,Findex=016,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE1e <- data.frame(p=0.50,n=400,Findex=011,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE1f <- data.frame(p=0.50,n=400,Findex=013,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE1g <- data.frame(p=0.50,n=400,Findex=015,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE1h <- data.frame(p=0.50,n=400,Findex=017,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE2  <- data.frame(p=0.95,n=1e3,Findex=010,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE3  <- data.frame(p=0.05,n=5e3,Findex=022,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE4  <- data.frame(p=0.50,n=6e3,Findex=034,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE5  <- data.frame(p=0.50,n=1e3,Findex=030,ALPHA=0.05,NUM.X0=NUM0,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASE6  <- data.frame(p=0.95,n=1e3,Findex=020,ALPHA=0.05,NUM.X0=5   ,TRIM.X0=TRIM0,XUNIF.FLAG=TRUE,HET.MULT=1,sig0=0.2)
CASES <- rbind(CASE1a) #exploration
NREPLIC <- 200; BREPLIC <- 2
if (SAVE.FLAG) { #paper
  CASES <- rbind(CASE1a,CASE1b,CASE1c,CASE1d,
                 CASE2,CASE3,CASE4,CASE5,CASE6)
  BREPLIC <- 299
  NREPLIC <- 5000
}
#
QTE.NO.BS <- TRUE #TRUE: always use Kaplan rather than BS, no matter n$t/n$c
QTE.GK.ONLY <- FALSE #TRUE: don't use Kaplan, just put NA
QTE.NO.GK <- FALSE #TRUE: never use GK (for comparison...)
LINEADJ <- 0.95;  CUBEADJ <- 1.75;  QEST.FLAG <- FALSE
HMULT <- 1;  DENOM.MIN.fn <- function(n) { (1e6/n)^(1/3) }
HEST.FLAG <- TRUE;  LOCAL <- TRUE
nms <- c("IDEAL1h","line1h","cube1h","line1hC91","cube1hC91") #quad1h
#          ,"IDEAL2h","quad2h","cube2h"
#       ,"BScube1h","BScube2h"
#       ,"BStcube1h","BStcube2h"
polynms <- nms[c(which(substr(nms,1,4)=="quad"),which(substr(nms,1,4)=="cube"),which(substr(nms,1,4)=="line"))]
#
#For the non-IDEAL method (currently local poly).  Could change: h, f0hat, normal vs. BS
other.h.fn <- function(n,d,sQ,h.IDEAL) h.IDEAL * n^(-1/(sQ+d)) / n^(-7/(7*d+18))
#
d <- 1 #only continuous var is scalar X
#

# bstatfn <- function(data,indices,deg,h,p,Qinv22,fbw,stud) {
#   Yloc <- data[indices,1];  Xloc <- data[indices,2]
#   return(tryCatch({
#     if (deg==2) {
#       X1 <- Xloc/h;  X2 <- Xloc^2/h^2
#       tmp <- rq(Yloc~X1+X2,tau=p) #auto-includes constant term
#     } else if (deg==3) {
#       X1 <- Xloc/h;  X2 <- Xloc^2/h^2;  X3 <- Xloc^3/h^3
#       tmp <- rq(Yloc~X1+X2+X3,tau=p)
#     }
#     est <- tmp$coefficients[2] / h
#     if (!stud) { return(est) }
#     f0hat <- fitted(npcdens(bws=fbw,tau=p,cxkertype='gaussian',cxkorder=2,cykertype='gaussian',cykorder=2,exdat=0,eydat=tmp$coefficients[1]))
#     other.var <- p*(1-p) * f0hat^(-2) * Qinv22  #    other.sd <- sqrt(other.var[2,2])
#     other.SE <- sqrt(other.var) / (h*sqrt(length(Yloc)))
#     return(c(est,other.SE^2))
# #     tstat <- (est-H0)/other.SE }
#     },
#     error=function(rqERR) { return(NA) }
#   ))
# }

CUBE1hopt <- QUAD1hopt <- LINE1hopt <- lineC91hopt <- cubeC91hopt <- FALSE #CPE-opt (vs. multiple of GK)
# for (n in ns) {
#   for (Findex in Findices) {
#     for (p in ps) {
for (icase in 1:dim(CASES)[1]) {
  set.seed(112358) #for replicability
  list[p,n,Findex,ALPHA,NUM.X0,TRIM.X0,XUNIF.FLAG,HET.MULT,sig0] <- CASES[icase,]
      NORM.APPROX <- FALSE
      NORM.APPROX <- (p==0.5)
      GAMMA.EST <- FALSE
      POLYINF <- 'boot'  #nid,rank,iid,ker,boot
      CV <- qnorm(p=1-ALPHA/2,mean=0,sd=1)
  
 #For any graphs
FILENAME.SUFFIX <- sprintf("_F%03d_a%02g_n%d_p%02d_het%g_xunif%d_numx%d_trimx%02d_nrep%d",Findex,100*ALPHA,n,as.integer(round(100*p)),HET.MULT,as.integer(XUNIF.FLAG),NUM.X0,as.integer(100*TRIM.X0),NREPLIC)

cat(sprintf("It is %s",format(Sys.time(), "%X, %A, %d %b %Y")),sep='\n')
if (OUTFILE!="") {
  cat(sprintf("\nIt is %s",format(Sys.time(), "%X, %A, %d %b %Y")),
      sep='\n',append=TRUE,file=OUTFILE)
}
cat(sprintf("Findex=%03d, HET.MULT=%g, sig0=%g, ALPHA=%g, n=%d, p=%g, numx=%d, trimx=%g, NREPLIC=%d, XUNIF.FLAG=%s, HEST.FLAG=%s, LOCAL=%s, NORM.APPROX=%s, GAMMA.EST=%s, QTE.NO.GK=%s, QTE.NO.BS=%s, QTE.GK.ONLY=%s, HMULT=%g, POLYINF=%s, BREPLIC=%d\n",
            Findex,HET.MULT,sig0,ALPHA,n,p,NUM.X0,TRIM.X0,NREPLIC,XUNIF.FLAG,HEST.FLAG,LOCAL,NORM.APPROX,GAMMA.EST,QTE.NO.GK,QTE.NO.BS,QTE.GK.ONLY,HMULT,POLYINF,BREPLIC))
if (OUTFILE!="") {
  cat(sprintf("Findex=%03d, HET.MULT=%g, sig0=%g, ALPHA=%g, n=%d, p=%g, numx=%d, trimx=%g, NREPLIC=%d, XUNIF.FLAG=%s, HEST.FLAG=%s, LOCAL=%s, NORM.APPROX=%s, GAMMA.EST=%s, QTE.NO.GK=%s, QTE.NO.BS=%s, QTE.GK.ONLY=%s, HMULT=%g, POLYINF=%s, BREPLIC=%d\n",
              Findex,HET.MULT,sig0,ALPHA,n,p,NUM.X0,TRIM.X0,NREPLIC,XUNIF.FLAG,HEST.FLAG,LOCAL,NORM.APPROX,GAMMA.EST,QTE.NO.GK,QTE.NO.BS,QTE.GK.ONLY,HMULT,POLYINF,BREPLIC),
    sep='\n',append=TRUE,file=OUTFILE)
}

x0.ps <- TRIM.X0 + (1-2*TRIM.X0)*c(0:(NUM.X0-1))/(NUM.X0-1) #X points of interest, as quantiles of X dist'n
X.MAX <- ifelse(XUNIF.FLAG,10,1)
if (XUNIF.FLAG) {
  rfn <- function(n) { runif(n,0,X.MAX) }
  f.X.fn <- function(Xs) rep.int(1/X.MAX,length(Xs))
  fp.X.fn <- function(Xs) rep.int(0,length(Xs))
  x0s <- qunif(x0.ps,0,X.MAX)
} else {
  rfn <- function(n) { rnorm(n,0,1) }
  X.MAX <- 1 #used for scaling later; =1 means don't scale
  f.X.fn <- function(Xs) dnorm(Xs)
  fp.X.fn <- function(Xs) -Xs*dnorm(Xs)
  x0s <- qnorm(x0.ps,0,1)
}
ns <- n
Xs <- replicate(max(ns),rfn(NREPLIC)) #NREPLIC rows, n cols; X~U(0,X.MAX)
unirands <- replicate(max(ns),runif(NREPLIC,0,1)) #NREPLIC rows, n cols; U(0,1)
if (NREPLIC==1) { 
  Xs <- matrix(Xs,1)
  unirands <- matrix(unirands,1) 
}
for (i in 1:NREPLIC) { #checking data not silly
  tmp <- range(Xs[i,1:min(ns)])
  for (k in 1:100) {
    if (tmp[1]>=x0s[1] || tmp[2]<=x0s[NUM.X0]) {
      Xs[i,] <- rfn(max(ns))
      tmp <- range(Xs[i,1:min(ns)])
    } else {
      break
    }
  }
  if (tmp[1]>=x0s[1] || tmp[2]<=x0s[NUM.X0]) {
    stop("Try another seed--not always drawing with x0s within range(X).")
  }
}

if (floor((Findex-100*floor(Findex/100))/10)==1) {
  Y.fn <- function(X) (X*(1-X))^(1/2)*sin(2*pi*(1+2^(-7/5))/(X+2^(-7/5)))
  Y.fn.expr.txt <- expression(sqrt(x(1-x))*sin*bgroup("(",2*pi*(1+2^{-7/5}) / (x+2^{-7/5}),")"))
  Y.fn.expr <- quote((X*(1-X))^(1/2)*sin(2*pi*(1+2^(-7/5))/(X+2^(-7/5))))
  Yd.fn <- function(X) ((1-2*X)*sin((2*(4+2^(3/5))*pi)/(4*X+2^(3/5))) + (16*(4+2^(3/5))*pi*(X-1)*X*cos((2*(4+2^(3/5))*pi)/(4*X+2^(3/5))))/(4*X+2^(3/5))^2) /
    (2*sqrt(X*(1-X))) #analytic derivative of Y.fn; matches Qp.fn
} else if (floor((Findex-100*floor(Findex/100))/10)==2) {
  if (XUNIF.FLAG) {
    Y.fn <- function(X) log(X)
    Yd.fn <- function(X) 1/X
    Y.fn.expr <- quote(log(X))
    Y.fn.expr.txt <- sprintf("log(X/%g)",X.MAX) #expression(log(X))
  } else {
    Y.fn <- function(X) exp(-X)
    Yd.fn <- function(X) -exp(-X)
    Y.fn.expr <- quote(exp(-X))
    Y.fn.expr.txt <- sprintf("exp(-X/%g)",X.MAX) 
  }
} else if (floor((Findex-100*floor(Findex/100))/10)==3) {
  slope <- floor(Findex/100)
  Y.fn <- function(X) slope*X #rep.int(0,length(X))
  Y.fn.expr <- bquote(.(slope)*X)
  Y.fn.expr.txt <- sprintf("%gX",slope) #expression(bquote(.(slope)*X))
  Yd.fn <- function(X) rep.int(slope,length(X))
} else {
  stop(sprintf("Not supported: Findex=%g",Findex))
}
y0s <- Y.fn(x0s/X.MAX) #no add'l effect of heterosk. since p-quantile of U is zero
yd0s <- Yd.fn(x0s/X.MAX)/X.MAX #no add'l effect of heterosk. since p-quantile of U is zero
#
Qp.fn <- function(X) eval(D(Y.fn.expr,'X'))
Qpp.fn <- function(X) eval(D(D(Y.fn.expr,'X'),'X'))
Qppp.fn <- function(X) eval(D(D(D(Y.fn.expr,'X'),'X'),'X'))
stopifnot(max(abs(Yd.fn(x0s/X.MAX)/X.MAX-Qp.fn(x0s/X.MAX)/X.MAX))<0.00001)


Ferr <- Findex-10*floor(Findex/10)
#Homoskedastic if Findex is even, heteroskedastic if odd
if ((Ferr %% 2) == 0) {
  sig.fn <- function(X) rep.int(sig0,length(X)) #homoskedastic
  sigp.fn <-function(X) rep.int(0,length(X)) 
  sig.fn.expr <- bquote(.(sig0)*1)
} else {
  sig.fn <- function(X) sig0*(1+HET.MULT*abs(X)) #heteroskedastic
  sigp.fn <-function(X) rep.int(sig0,length(X)) 
  sig.fn.expr <- bquote(.(sig0)*(1+.(HET.MULT)*X)) #rqss: 1+X
  sig.fn.expr <- substitute(S0*group("(",1+HM*group("|",X,"|"),")"),list(S0=sig0,HM=HET.MULT)) #rqss: 1+X
}
#Different error distribution shapes
if (any(Ferr==c(0,1))) { #normal
  Us <- qnorm(unirands[,1:n],0,1) - qnorm(p,0,1)
  if (!HEST.FLAG) {
    f.Y.fn <- function(p,Xs) dnorm(qnorm(p,0,sig.fn(Xs)),0,sig.fn(Xs))
    Fp.fn.expr <- D(bquote(pnorm((Y-.(Y.fn.expr))/.(sig.fn.expr) +.(qnorm(p,0,1)))),'X')
    Fp.fn.chk <- function(y,x) dnorm((y-Y.fn(x))/sig.fn(x)) *(-sig.fn(x)*Qp.fn(x) -(y-Y.fn(x))*sigp.fn(x))/sig.fn(x)^2
  }
} else if (any(Ferr==c(2,3))) { #t3
  Us <- qt(unirands[,1:n],3) - qt(p,3)
  if (!HEST.FLAG) {
    f.Y.fn <- function(p,Xs) dt(qt(p,3),3)/sig.fn(Xs)
    Fp.fn.expr <- D(bquote(pt((Y-.(Y.fn.expr))/.(sig.fn.expr) +.(qt(p,3)),3)),'X')
  }
} else if (any(Ferr==c(4,5))) { #t1 (Cauchy)
  Us <- qt(unirands[,1:n],1) - qt(p,1)
  if (!HEST.FLAG) {
    f.Y.fn <- function(p,Xs) dt(qt(p,1),1)/sig.fn(Xs)
    Fp.fn.expr <- D(bquote(pt((Y-.(Y.fn.expr))/.(sig.fn.expr) +.(qt(p,1)),1)),'X')
  }
} else if (any(Ferr==c(6,7))) { #chi2(3)
  Us <- qchisq(unirands[,1:n],3) - qchisq(p,3)
} else stop(sprintf("Not implemented: Findex=%g",Findex))
if (!HEST.FLAG) {
  Fp.fn <- function(Y,X) eval(Fp.fn.expr)
  Fpp.fn.expr <- D(Fp.fn.expr,'X')
  Fpp.fn <- function(Y,X) eval(Fpp.fn.expr)
  Fppp.fn.expr <- D(Fpp.fn.expr,'X')
  Fppp.fn <- function(Y,X) eval(Fppp.fn.expr)
}

# Calculations for infeasible bandwidth
if (!HEST.FLAG) {
  if (Findex==030) {
    hs.opt <- rep.int(Inf,NUM.X0)
  } else {
    f.Xs <- f.X.fn(x0s) #already accounts for X.MAX
    fp.Xs <- fp.X.fn(x0s) #already accounts for X.MAX
    Qps <- Qp.fn(x0s/X.MAX)/X.MAX
    stopifnot(max(abs(Qps-yd0s))<0.1) #check
    Qpps <- Qpp.fn(x0s/X.MAX)/X.MAX^2
#     f.Ys <- f.Y.fn(p,x0s) #this is only approximate, but hopefully good enough.
#     Fps.chk <- -f.Ys*Qps #Fp.fn(x0s)
    Fps <- Fp.fn(y0s,x0s/X.MAX)/X.MAX #in case X.MAX>1
#     stopifnot(max(abs(Fps.chk-Fps))<0.1) #check Fp.fn formulated properly
    Fpps <- Fpp.fn(y0s,x0s/X.MAX)/X.MAX^2
  zz1.pt <- qnorm(1-ALPHA); zz2.pt <- qnorm(1-ALPHA/2)
    hs.pt <- as.vector(quantile.inf.np.bandwidths(p=p,n=n,X=c(-Inf,Inf),x0s=x0s,f.X.hats=f.Xs,fp.X.hats=fp.Xs,Fp.hats=Fps,Fpp.hats=Fpps,ONESIDED=ONESIDED,zz1=ifelse(ONESIDED==0,zz2.pt,zz1.pt),DENOM.MIN=DENOM.MIN.fn(n)/X.MAX^3)) #X only used for max/min--done later in NREPLIC loop anyway (and need to divide by 2 first, too).
    hs.pt <- hs.pt * n^(4/((18+7*d)*(2+d))) #different optimal rate, per Kaplan (2014)
    hs.pt <- quantile.inf.np.bandwidth.correct(hs=hs.pt,x0s=x0s,HCOR=HCOR)
    hs.pt <- hs.pt / 2  #since together windows span x0-2h to x0+2h, not just +/-h
    hs.opt <- hs.pt * HMULT
    rm(hs.pt)
  }
}

#Now can generate Y data, depending on U dist'n
Ys <- Y.fn(Xs[,1:n]/X.MAX) + pmax(sig.fn(Xs[,1:n]/X.MAX),0)*Us
if (NREPLIC==1) { Ys <- matrix(Ys,1) }

#Storage variables for big loop
stores <- NULL
tmp <- array(data=NA,dim=c(NREPLIC,NUM.X0))
for (nm in nms) {
  stores[[nm]] <- list(hs=tmp,CI.lo=tmp,CI.hi=tmp,name=nm)
}
Q3 <- cbind(c(1,0,1/3,0),c(0,1/3,0,1/5),c(1/3,0,1/5,0),c(0,1/5,0,1/7)) #Chaudhuri(1991, p. 765)
Qinv3 <- solve(Q3) #cubic
Q2 <- Q3[1:3,1:3] #cbind(c(1,0,2/3),c(0,2/3,0),c(2/3,0,2/5))
Qinv2 <- solve(Q2) #quadratic
Q1 <- Q3[1:2,1:2];  Qinv1 <- solve(Q1)
#
tmp <- matrix(0,NREPLIC,NUM.X0)
methref <- c("GoldmanKaplan","Kaplan","bootstrap-perc-t")
meth.cnt <- list(GK=tmp,K=tmp,BS=tmp,NONE=tmp)
GKNLs <- GKNRs <- matrix(NA,NREPLIC,NUM.X0)
bs.total.time <- IDEAL.total.time <- IDEALnb.total.time <- Sys.time()-Sys.time()
########### NREPLIC LOOP #############
LOOP.STARTTIME <- Sys.time()
for (irep in 1:NREPLIC) {
  cat(sprintf("rep %d of %d (%s since loop start at %s)\r",irep,NREPLIC,capture.output(print(Sys.time()-LOOP.STARTTIME)),LOOP.STARTTIME))
  flush.console() #otherwise, sometimes won't print until loop is done...
  X <- Xs[irep,1:n];  Y <- Ys[irep,1:n]
  if (!HEST.FLAG) {
    MAXHS <- apply(cbind(x0s-min(X),max(X)-x0s),1,min) #assuming scalar regressor
    hs <- apply(cbind(MAXHS/2,hs.opt),1,min)
  } else {
    hs <- NULL
  }
  #
  IDEAL.start.time <- Sys.time()
  useNA <- tryCatch({
    tmp <- qme.inf(Y=Y,X=X,x0s=x0s,p=p,ALPHA=ALPHA,NORM.APPROX=NORM.APPROX,HMULT=HMULT,hs=hs,GAMMA.EST=GAMMA.EST, QTE.GK.ONLY=QTE.GK.ONLY, QTE.NO.BS=QTE.NO.BS, QTE.NO.GK=QTE.NO.GK, LOCAL=LOCAL, DENOM.MIN=DENOM.MIN.fn(n))
    FALSE}, error=function(ERR) {TRUE})
  if (useNA || !any(!is.na(tmp$CI.lo))) {
    meth.cnt[['NONE']][irep,] <- 1
    next
  }
  IDEAL.elapsed.time <- (Sys.time()-IDEAL.start.time)
  IDEAL.total.time <- IDEAL.total.time + IDEAL.elapsed.time
  # Time *without* bandwidth calculation (since required for local poly, too)
  IDEALnb.start.time <- Sys.time()
#   qme.inf(Y=Y,X=X,x0s=x0s,p=p,ALPHA=ALPHA,NORM.APPROX=NORM.APPROX,HMULT=HMULT,hs=tmp$hs,GAMMA.EST=GAMMA.EST, QTE.GK.ONLY=QTE.GK.ONLY, QTE.NO.BS=QTE.NO.BS, QTE.NO.GK=QTE.NO.GK) #LOCAL ignored since tmp$hs provided
  IDEALnb.elapsed.time <- (Sys.time()-IDEALnb.start.time)
  IDEALnb.total.time <- IDEALnb.total.time + IDEALnb.elapsed.time
  #
#   if (is.na(tmp$CI.lo[1]*tmp$CI.lo[length(tmp$CI.lo)])) { next } #8s/skip-rep...not actually much shorter
  for (im in 1:length(methref)) {
    meth.cnt[[im]][irep,methref[im]==tmp$methname] <- 1
  }
  meth.cnt[['NONE']][irep,] <- as.integer(is.na(tmp$methname))
  stores[["IDEAL1h"]]$hs[irep,] <- tmp$hs  #IDEAL CPE-optimal plug-in
  stores[["IDEAL1h"]]$CI.lo[irep,] <- tmp$CI.lo
  stores[["IDEAL1h"]]$CI.hi[irep,] <- tmp$CI.hi
  GKNLs[irep,] <- ifelse(is.na(tmp$CI.lo),NA,tmp$NL)
  GKNRs[irep,] <- ifelse(is.na(tmp$CI.lo),NA,tmp$NR)
  #
  #Other bandwidth setting
  lineh <- lineC91h <- LINEADJ*2*stores[["IDEAL1h"]]$hs[irep,] 
  cubeh <- cubeC91h <- CUBEADJ*2*stores[["IDEAL1h"]]$hs[irep,] 
  quadh <- 3*stores[["IDEAL1h"]]$hs[irep,] #### XXX ######
  lineopth <- 2*other.h.fn(n=n,d=d,sQ=2,h.IDEAL=stores[["IDEAL1h"]]$hs[irep,])
  quadopth <- 2*other.h.fn(n=n,d=d,sQ=3,h.IDEAL=stores[["IDEAL1h"]]$hs[irep,])
  cubeopth <- 2*other.h.fn(n=n,d=d,sQ=4,h.IDEAL=stores[["IDEAL1h"]]$hs[irep,])
  if (LINE1hopt) { lineh <- lineopth }
  if (QUAD1hopt) { lineh <- quadopth }
  if (CUBE1hopt) { lineh <- cubeopth }
  if (lineC91hopt)   { lineC91h  <- lineopth }
  if (cubeC91hopt)   { cubeC91h  <- cubeopth }
  #make sure h doesn't extend past min or max, to prevent bias
  maxhs <- apply(cbind(x0s-min(X),max(X)-x0s),1,min)
  lineh <- apply(cbind(maxhs,lineh),1,min)
  quadh <- apply(cbind(maxhs,quadh),1,min)
  cubeh <- apply(cbind(maxhs,cubeh),1,min)
  lineC91h <- apply(cbind(maxhs,lineC91h),1,min)
  cubeC91h <- apply(cbind(maxhs,cubeC91h),1,min)
  if ("quad2h" %in% nms) { stores[["quad2h"]]$hs[irep,] <- quadh }
  if ("cube2h" %in% nms) { stores[["cube2h"]]$hs[irep,] <- cubeh }
  if ("cube1h" %in% nms) { stores[["cube1h"]]$hs[irep,] <- cubeh }
  if ("quad1h" %in% nms) { stores[["quad1h"]]$hs[irep,] <- quadh }
  if ("line1h" %in% nms) { stores[["line1h"]]$hs[irep,] <- lineh }
  if ("line1hC91" %in% nms) { stores[["line1hC91"]]$hs[irep,] <- lineC91h }
  if ("cube1hC91" %in% nms) { stores[["cube1hC91"]]$hs[irep,] <- cubeC91h }
  #
  if ("IDEAL2h" %in% nms) { 
    stores[["IDEAL2h"]]$hs[irep,] <- cubeopth / 2 #to use same data
    tmp <- qme.inf(Y=Y,X=X,x0s=x0s,p=p,ALPHA=ALPHA,NORM.APPROX=NORM.APPROX,
                   hs=stores[["IDEAL2h"]]$hs[irep,],GAMMA.EST=GAMMA.EST,
                   QTE.GK.ONLY=QTE.GK.ONLY, QTE.NO.BS=QTE.NO.BS)
    stores[["IDEAL2h"]]$CI.lo[irep,] <- tmp$CI.lo
    stores[["IDEAL2h"]]$CI.hi[irep,] <- tmp$CI.hi
  }
  #
  #Local polynomial calculations (see Chaudhuri 1991)
  bs.start.time <- Sys.time()
  for (ix in 1:NUM.X0) {
   f0hat.bw <- NA
   for (nm in polynms) {
    nmsuff <- substr(nm,nchar(nm)-2,nchar(nm))
    if (POLYINF=='boot' && BREPLIC==1 && nmsuff!="C91") {
      stores[[nm]]$CI.lo[irep,ix] <- stores[[nm]]$CI.hi[irep,ix] <- NA
      next
    } else {
      BSnm <- paste0("BS",nm,collapse="")
      BStnm <- paste0("BSt",nm,collapse="")
      h <- stores[[nm]]$hs[irep,ix]
      iloc <- (x0s[ix]-h<=X) & (X<=x0s[ix]+h)
      Yloc <- Y[iloc];  Xloc <- X[iloc]-x0s[ix]  #center X s.t. x0=0, as in Chaudhuri (1991)
      if (length(Yloc)<2) {
        stores[[nm]]$CI.lo[irep,ix] <- stores[[nm]]$CI.hi[irep,ix] <- NA
        next
      }
      Xdes <- cbind(1,Xloc/h,Xloc^2/h^2,Xloc^3/h^3,Xloc^4/h^4,Xloc^5/h^5,Xloc^6/h^6)
      EXs <- colMeans(Xdes)
      Q3hat <- cbind(c(1,EXs[2],EXs[3],EXs[4]),c(EXs[2],EXs[3],EXs[4],EXs[5]),c(EXs[3],EXs[4],EXs[5],EXs[6]),c(EXs[4],EXs[5],EXs[6],EXs[7])) #estimate Qn of Q in Chaudhuri(1991, p. 765)
      Q3hatinv <- tryCatch(solve(Q3hat),error=function(ER) Qinv3)
      useNA <- tryCatch({
        if (substr(nm,1,4)=="quad") {
          rq.ret <- rq(Yloc~Xdes[,2]+Xdes[,3],tau=p) #auto-includes constant term
          Qinv <- Qinv2;  Qinvhat <- Q3hatinv[1:3,1:3]
          deg <- 2
        } else if (substr(nm,1,4)=="line") {
          rq.ret <- rq(Yloc~Xdes[,2],tau=p) #auto-includes constant term
          Qinv <- Qinv1;  Qinvhat <- Q3hatinv[1:2,1:2]
          deg <- 1
        } else if (substr(nm,1,4)=="cube") {
          rq.ret <- rq(Yloc~Xdes[,2]+Xdes[,3]+Xdes[,4],tau=p)
          Qinv <- Qinv3;  Qinvhat <- Q3hatinv
          deg <- 3
        } else {
          stop("Only linear, quadratic, and cubic local polynomials implemented.")
        }
        FALSE},
        error=function(rqERR) { return(TRUE) }
      )
      if (useNA) {
        stores[[nm]]$CI.lo[irep,ix] <- stores[[nm]]$CI.hi[irep,ix] <- NA
        if (any(nms==BSnm)) {
          stores[[BSnm]]$CI.lo[irep,ix] <- stores[[BSnm]]$CI.hi[irep,ix] <- NA
        }
        if (any(nms==BStnm)) {
          stores[[BStnm]]$CI.lo[irep,ix] <- stores[[BStnm]]$CI.hi[irep,ix] <- NA
        }
      } else {
        other.est <- rq.ret$coefficients[2] / h
        if (nmsuff=="C91") {
#           f0hat.bw <- npcdensbw(formula=Yloc~Xloc,bwmethod="normal-reference",ckertype='gaussian',ckerorder=2)
#           f0hat <- fitted(npcdens(bws=f0hat.bw,tau=p,cxkertype='gaussian',cxkorder=2,cykertype='gaussian',cykorder=2,exdat=0,eydat=rq.ret$coefficients[1]))
          f0hat <- kde(Yloc,eval.points=rq.ret$coefficients[1],h=hns(Yloc,deriv.order=0))$estimate
          # #         Uhat <- Yloc - predict.rq(tmp)
          # #         f0hat <- kde(Uhat,eval.points=0,h=hns(Uhat,deriv.order=0))$estimate
          if (QEST.FLAG) { Q <- Qinvhat } else { Q <- Qinv }
          other.var <- p*(1-p) * f0hat^(-2) * Q  #C91 Prop 4.2
          other.SE <- sqrt(other.var[2,2]) / (h*sqrt(length(Yloc)))
          stores[[nm]]$CI.lo[irep,ix] <- other.est - CV * other.SE
          stores[[nm]]$CI.hi[irep,ix] <- other.est + CV * other.SE 
        } else {
          tryCatch({if (POLYINF=='rank') {
            stores[[nm]]$CI.lo[irep,ix] <- summary(rq.ret,se="rank")$coefficients[2,2]
            stores[[nm]]$CI.hi[irep,ix] <- summary(rq.ret,se="rank")$coefficients[2,3]
          } else {
            other.SE <- summary(rq.ret,se=POLYINF,R=BREPLIC)$coefficients[2,2] / h
            stores[[nm]]$CI.lo[irep,ix] <- other.est - CV * other.SE
            stores[[nm]]$CI.hi[irep,ix] <- other.est + CV * other.SE 
          }
          }, error=function(ERR){
            #cat("ERR",file="");
            warning(ERR$message)
            stores[[nm]]$CI.lo[irep,ix] <- stores[[nm]]$CI.hi[irep,ix] <- NA
            #                    tryCatch({other.SE <- summary(rq.ret,se="boot",R=BREPLIC)$coefficients[2,2] / h
            #                              stores[[nm]]$CI.lo[irep,ix] <- other.est - CV * other.SE
            #                              stores[[nm]]$CI.hi[irep,ix] <- other.est + CV * other.SE}, 
            #                             error=function(ERR){cat("err",file="");warning(ERR$message)})
          })
          #
          if (any(nms==BSnm) && !any(nms==BStnm)) {
            bootobj <- boot(data=cbind(Yloc,Xloc),statistic=bstatfn,R=BREPLIC,
                            deg=deg,h=h,p=p,Qinv22=Qinv[2,2],fbw=f0hat.bw,stud=FALSE)
            bci <- boot.ci(bootobj,conf=1-ALPHA,type="basic")
            stores[[BSnm]]$CI.lo[irep,ix] <- bci$basic[1,4]
            stores[[BSnm]]$CI.hi[irep,ix] <- bci$basic[1,5]
          }
          if (any(nms==BStnm)) {
            bootobj <- boot(data=cbind(Yloc,Xloc),statistic=bstatfn,R=BREPLIC,
                            deg=deg,h=h,p=p,Qinv22=Qinv[2,2],fbw=f0hat.bw,stud=TRUE)
            bci <- boot.ci(bootobj,conf=1-ALPHA,type=c("stud","basic"))
            stores[[BStnm]]$CI.lo[irep,ix] <- bci$student[1,4]
            stores[[BStnm]]$CI.hi[irep,ix] <- bci$student[1,5]
            if (any(nms==BSnm)) {
              stores[[BSnm]]$CI.lo[irep,ix] <- bci$basic[1,4]
              stores[[BSnm]]$CI.hi[irep,ix] <- bci$basic[1,5]
            }
          }
        }
      }
    }
   }
  }
  bs.time.elapsed <- Sys.time()-bs.start.time
  bs.total.time <- bs.total.time+bs.time.elapsed
} #end NREPLIC loop
#
cat(sprintf("Total time for BS, IDEAL (non-bandwidth), IDEAL (total): %s, %s, %s\n",format(bs.total.time),format(IDEALnb.total.time),format(IDEAL.total.time)))
if (OUTFILE!="") {
  cat(sprintf("Total time for BS, IDEAL (non-bandwidth), IDEAL (total): %s, %s, %s\n",format(bs.total.time),format(IDEALnb.total.time),format(IDEAL.total.time)),
      append=TRUE,sep='',file=OUTFILE)
}
#
nonNA1h <- base::colSums(x=!is.na(stores[["IDEAL1h"]]$CI.lo))
cat(paste0("1h non-NA: ",paste0(sprintf("%d(%d,%d,%d) ",nonNA1h,colSums(meth.cnt$GK),colSums(meth.cnt$K),colSums(meth.cnt$BS)),collapse=""),collapse=""),
    file=OUTFILE,sep='\n',append=TRUE)
if (prod(1==(meth.cnt$GK+meth.cnt$K+meth.cnt$BS+meth.cnt$NONE))!=1) {
  cat("CHECK meth.cnt: NOT EXHAUSTIVE")
}
if ("IDEAL2h" %in% nms) {
  nonNA2h <- base::colSums(x=!is.na(stores[["IDEAL2h"]]$CI.lo))
  cat(paste0("2h non-NA: ",paste0(sprintf("%d ",nonNA2h),collapse=""),collapse=""),
      file=OUTFILE,sep='\n',append=TRUE)
}
#
cat(paste0("gk.med.hs: ",paste0(sprintf("%6.3f ",apply(stores[["IDEAL1h"]]$hs,2,median,na.rm=TRUE)),collapse=""),collapse=""),
    file=OUTFILE,sep='\n',append=TRUE)
cat(paste0("gk.p90.hs: ",paste0(sprintf("%6.3f ",apply(stores[["IDEAL1h"]]$hs,2,quantile,na.rm=TRUE,p=0.90,type=7)),collapse=""),collapse=""),
    file=OUTFILE,sep='\n',append=TRUE)
#
cat(paste0("gk.med.NL: ",paste0(sprintf("%6g ",apply(GKNLs,2,quantile,na.rm=TRUE,p=0.5,type=1)),collapse=""),collapse=""),
    file=OUTFILE,sep='\n',append=TRUE)
#
cat(paste0("gk.med.NR: ",paste0(sprintf("%6g ",apply(GKNRs,2,quantile,na.rm=TRUE,p=0.5,type=1)),collapse=""),collapse=""),
    file=OUTFILE,sep='\n',append=TRUE)
#
#Calculate and output coverage probabilities (CP), median CI length
hstrs <- "1h"
lenmax <- 0
for (nm in nms[substr(nms,1,3) %in% c("IDE","lin")]) {
  lens <- stores[[nm]]$CI.hi-stores[[nm]]$CI.lo
  lens <- ifelse(is.na(stores[["IDEAL1h"]]$CI.lo),NA,lens)
  if ("IDEAL2h" %in% nms) { 
    lens <- ifelse(is.na(stores[["IDEAL2h"]]$CI.lo),NA,lens)
  }
  medlen <- apply(X=lens,MARGIN=2,FUN=median,na.rm=TRUE)
#   ml <- apply(X=stores[[nm]]$CI.hi-stores[[nm]]$CI.lo,MARGIN=2,FUN=median,na.rm=TRUE)
  medlen <- ifelse(medlen==Inf,0,medlen)
  if (any(!is.na(medlen))) {
    lenmax <- max(lenmax,max(medlen,na.rm=TRUE),na.rm=TRUE)
  }
}
lenmax <- 1.1 * lenmax
for (hstr in hstrs) {
  if (SAVE.FLAG) pdf(file=sub(".txt",paste0("_medlen",FILENAME.SUFFIX,".pdf"),OUTFILE)) #else draw to current device
  par(family='serif',mar=c(5.0,6.0,6.0,2.1))
  tmprx <- range(x0s)
  tmprx <- round(c(sum(tmprx)/2-(tmprx[2]-tmprx[1])/2/0.92,sum(tmprx)/2+(tmprx[2]-tmprx[1])/2/0.92), digits=2)
  plot(tmprx,c(0,lenmax),type='n', 
       cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
       main=paste0(hstr,FILENAME.SUFFIX,collapse=""), 
       xlab="X", ylab="Median CI Length")
  leg <- NULL;  ltys <- NULL
#   for (nm in nms[substr(nms,nchar(nms)-1,nchar(nms))==hstr]) {
  for (nm in nms) {
    if (!length(grep(hstr,nms))>0) next
    pch <- NA
    lens <- stores[[nm]]$CI.hi-stores[[nm]]$CI.lo
    lens <- ifelse(is.na(stores[[paste0("IDEAL",hstr,collapse="")]]$CI.lo),NA,lens)
    medlen <- apply(X=lens,MARGIN=2,FUN=median,na.rm=TRUE)
    medlen <- ifelse(medlen==Inf,lenmax*1000,medlen)
    f4 <- substr(nm,1,3)
    NOPRINT <- FALSE
    if (f4=="IDE") {
      lines(x0s,medlen,lwd=2,lty=1,pch=pch,col=1,type='o',cex=0.7)
      leg <- c(leg,"New");  ltys <- c(ltys,1)
    } else if (f4=="qua") {
      lines(x0s,medlen,lwd=2,lty=2,pch=pch,col=1,type='o',cex=0.7)
      leg <- c(leg,sprintf("Quad/%s",POLYINF));  ltys <- c(ltys,2)
    } else if (f4=="cub") {
      if (substr(nm,nchar(nm)-2,nchar(nm))=="C91") {
        lines(x0s,medlen,lwd=2,lty=6,pch=pch,col=1,type='o',cex=0.7)
        leg <- c(leg,sprintf("Cubic/C91"));  ltys <- c(ltys,6)
      } else if (POLYINF!='boot' || BREPLIC>1) {
        lines(x0s,medlen,lwd=2,lty=4,pch=pch,col=1,type='o',cex=0.7)
        leg <- c(leg,sprintf("Cubic/%s",POLYINF));  ltys <- c(ltys,4)
      } else {
        NOPRINT <- TRUE
      }
    } else if (f4=="lin") {
      if (substr(nm,nchar(nm)-2,nchar(nm))=="C91") {
        lines(x0s,medlen,lwd=2,lty=2,pch=pch,col=1,type='o',cex=0.7)
        leg <- c(leg,sprintf("Linear/C91"));  ltys <- c(ltys,2)
      } else if (POLYINF!='boot' || BREPLIC>1) {
        lines(x0s,medlen,lwd=2,lty=3,pch=pch,col=1,type='o',cex=0.7)
        leg <- c(leg,sprintf("Linear/%s",POLYINF));  ltys <- c(ltys,3)
      } else {
        NOPRINT <- TRUE
      }
    } else if (f4=="BSc") {
      lines(x0s,medlen,lwd=2,lty=5,pch=pch,col=1,type='o',cex=0.7)
      leg <- c(leg,"Cubic/boot");  ltys <- c(ltys,5)
    }
    if (!NOPRINT) {
      #     cat(nm,sep='\n',append=TRUE,file=OUTFILE)
      cat(paste0(nm,"len <- c(",
                 paste0(sprintf("%6.2f",medlen),collapse=","),")",collapse=""),
          sep='\n',append=TRUE,file=OUTFILE)
    }
  }
  legend("top", leg, inset=c(0.0,0.0), bty='n',
         col=1, pch=NA, lty=ltys, lwd=2, cex=1.3)
  if (SAVE.FLAG) dev.off()
}
#
if (any(substr(nms,nchar(nms)-1,nchar(nms))=="2h")) { hstrs <- c(hstrs,"2h") }
for (hstr in hstrs) {
  if (SAVE.FLAG) pdf(file=sub(".txt",paste0("_CP",FILENAME.SUFFIX,".pdf"),OUTFILE)) #else draw to current device
  par(family='serif',mar=c(5.0,6.0,6.0,2.1))
  tmprx <- range(x0s)
  tmprx <- round(c(sum(tmprx)/2-(tmprx[2]-tmprx[1])/2/0.92,sum(tmprx)/2+(tmprx[2]-tmprx[1])/2/0.92), digits=2)
  plot(tmprx,0:1,type='n', 
       ylim=c(0.8,1.0), #JBES R1 suggestion
       cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
       main=paste0(hstr,FILENAME.SUFFIX,collapse=""), 
       xlab="X", ylab="Coverage Probability")
  lines(tmprx,c(1-ALPHA,1-ALPHA),lty=1,lwd=1,col=1) #nominal
  leg <- NULL;  ltys <- NULL
#   for (nm in nms[substr(nms,nchar(nms)-1,nchar(nms))==hstr]) {
  for (nm in nms) {
    if (!length(grep(hstr,nms))>0) next
    cov <- (stores[[nm]]$CI.lo<=(matrix(1,NREPLIC,1)%x%t(yd0s))) &
      (stores[[nm]]$CI.hi>=(matrix(1,NREPLIC,1)%x%t(yd0s)))
    cov <- ifelse(is.na(stores[[paste0("IDEAL",hstr,collapse="")]]$CI.lo),NA,cov)
    cp <- colMeans(cov,na.rm=TRUE)
    pch <- NA
    f4 <- substr(nm,1,3)
    NOPRINT <- FALSE
    if (f4=="IDE") {
      lines(x0s,cp,lwd=2,lty=1,pch=pch,col=1,type='o',cex=0.7)
      leg <- c(leg,"New");  ltys <- c(ltys,1)
    } else if (f4=="qua") {
      lines(x0s,cp,lwd=2,lty=2,pch=pch,col=1,type='o',cex=0.7)
      leg <- c(leg,sprintf("Quad/%s",POLYINF));  ltys <- c(ltys,2)
    } else if (f4=="cub") {
      if (substr(nm,nchar(nm)-2,nchar(nm))=="C91") {
        lines(x0s,cp,lwd=2,lty=6,pch=pch,col=1,type='o',cex=0.7)
        leg <- c(leg,sprintf("Cubic/C91"));  ltys <- c(ltys,6)
      } else if (POLYINF!='boot' || BREPLIC>1) {
        lines(x0s,cp,lwd=2,lty=4,pch=pch,col=1,type='o',cex=0.7)
        leg <- c(leg,sprintf("Cubic/%s",POLYINF));  ltys <- c(ltys,4)
      } else {
        NOPRINT <- TRUE
      }
    } else if (f4=="lin") {
      if (substr(nm,nchar(nm)-2,nchar(nm))=="C91") {
        lines(x0s,cp,lwd=2,lty=2,pch=pch,col=1,type='o',cex=0.7)
        leg <- c(leg,sprintf("Linear/C91"));  ltys <- c(ltys,2)
      } else if (POLYINF!='boot' || BREPLIC>1) {
        lines(x0s,cp,lwd=2,lty=3,pch=pch,col=1,type='o',cex=0.7)
        leg <- c(leg,sprintf("Linear/%s",POLYINF));  ltys <- c(ltys,3)
      } else {
        NOPRINT <- TRUE
      }
    } else if (f4=="BSc") {
      lines(x0s,cp,lwd=2,lty=5,pch=pch,col=1,type='o',cex=0.7)
      leg <- c(leg,"Cubic/boot");  ltys <- c(ltys,5)
    }
    if (!NOPRINT) {
      cat(paste0(nm,"cp <- c(",
                 paste0(sprintf("%5.3f",cp),collapse=","),")",collapse=""),
          sep='\n',append=TRUE,file=OUTFILE)
    }
  }
  legend("bottomright", leg, inset=c(0.0,0.0), bty='n',
         col=1, pch=NA, lty=ltys, lwd=2, cex=1.3)
  if (SAVE.FLAG) dev.off()
}
      cat("\n"); if (OUTFILE!="") {cat("\n",file=OUTFILE,append=TRUE)}

#       } #end of loop over p

    # GRAPH TRUE FUNCTION AND x0s
    if (TRUEFUNOUT.FLAG) {
      if (SAVE.FLAG)  pdf(file=sprintf("F%03d_NX%02d_TRIM%02d_true_fn.pdf",Findex,NUM.X0,as.integer(100*TRIM.X0)),pointsize=12, width=11, height=7)
      par(family="serif",mar=c(5.0,6.0,6.0,2.1))
      tmprx <- range(x0s)
      tmprx <- round(c(sum(tmprx)/2-(tmprx[2]-tmprx[1])/2/0.92,sum(tmprx)/2+(tmprx[2]-tmprx[1])/2/0.92), digits=2)
      xx <- seq(from=tmprx[1],to=tmprx[2],by=(tmprx[2]-tmprx[1])/1000)
      yy <- Y.fn(xx/X.MAX);  tmpry <- range(Y.fn(xx[xx>=x0s[1] & xx<=x0s[NUM.X0]]/X.MAX))
      plot(xx,yy,type="l",lwd=4,xlab="X",ylab=expression(Q[Y*symbol("|")*X](p*symbol(";")*x)),xlim=tmprx,ylim=tmpry, mgp=c(3,1,0),col=1,lty=1,main=Y.fn.expr.txt, cex.main=2,cex.lab=2,cex.axis=2)
      points(x0s,Y.fn(x0s/X.MAX),lwd=2,pch=1,col=1,cex=2)
      if (SAVE.FLAG) dev.off()
    }
#   } #end of loop over Findex
# } #end of loop over n
} #end of loop over CASE

tmpt <- as.numeric(Sys.time()-STARTTIME,units="secs")
tmps <- sprintf("Total time elapsed=%d seconds (i.e., %dmin; i.e., %4.2fhrs); Current time: %s\n",as.integer(tmpt),as.integer(round(tmpt/60)),tmpt/3600,format(Sys.time(), "%X, %A, %d %b %Y"))
cat(tmps)
if (OUTFILE!="") {
  cat(tmps,file=OUTFILE,sep="\n",append=TRUE)
}


if (FALSE) {
# TO REDRAW LENGTH GRAPH FROM TEXT RESULTS:
#first, run code from output to store IDEAL1hlen, etc.
#second, store x0s based on trim.x0 and num.x0
pdf(file="CUSTOM_LEN_GRAPH.pdf")
par(family='serif',mar=c(5.0,6.0,6.0,2.1))
plot(c(0,10),c(0,6),type='n', 
     cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
     main="CUSTOM LEN GRAPH", 
     xlab="X", ylab="Median CI Length")
lines(x0s,IDEAL1hlen,  lwd=2,lty=1,pch=NA,col=1,type='o',cex=0.7)
lines(x0s,line1hlen,   lwd=2,lty=3,pch=NA,col=1,type='o',cex=0.7)
lines(x0s,cube1hlen,   lwd=2,lty=4,pch=NA,col=1,type='o',cex=0.7)
lines(x0s,line1hC91len,lwd=2,lty=2,pch=NA,col=1,type='o',cex=0.7)
lines(x0s,cube1hC91len,lwd=2,lty=6,pch=NA,col=1,type='o',cex=0.7)
legend("top", c("New","Linear/boot","Cubic/boot","Linear/C91","Cubic/C91"), inset=c(0.0,0.0), bty='n',
       col=1, pch=NA, lty=c(1,3,4,2,6), lwd=2, cex=1.3)
dev.off()
#
# TO REDRAW CP GRAPH FROM TEXT RESULTS:
#first, run code from output to store IDEAL1hcp, etc.
#second, store x0s based on trim.x0 and num.x0
pdf(file="CUSTOM_CP_GRAPH.pdf")
par(family='serif',mar=c(5.0,6.0,6.0,2.1))
plot(c(0,10),c(0,1),type='n', 
     ylim=c(0.8,1.0), #JBES R1 suggestion
     cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
     main="CUSTOM CP GRAPH", 
     xlab="X", ylab="Coverage Probability")
lines(c(0,10),c(0.95,0.95),lty=1,col=1,lwd=1)
lines(x0s,IDEAL1hcp,  lwd=2,lty=1,pch=NA,col=1,type='o',cex=0.7)
lines(x0s,line1hcp,   lwd=2,lty=3,pch=NA,col=1,type='o',cex=0.7)
lines(x0s,cube1hcp,   lwd=2,lty=4,pch=NA,col=1,type='o',cex=0.7)
lines(x0s,line1hC91cp,lwd=2,lty=2,pch=NA,col=1,type='o',cex=0.7)
lines(x0s,cube1hC91cp,lwd=2,lty=6,pch=NA,col=1,type='o',cex=0.7)
legend("bottomright", c("New","Linear/boot","Cubic/boot","Linear/C91","Cubic/C91"), inset=c(0.0,0.0), bty='n',
       col=1, pch=NA, lty=c(1,3,4,2,6), lwd=2, cex=1.3)
dev.off()

} #END IF

# COMPUTATION TIME
comp.time.fn <- function(reps,NO.BS=FALSE,OUTFILE="") {
  require("quantreg");  source("qme_inf.R") # from https://github.com/kaplandm/R
  p <- 0.5;  BREPLIC <- 299;  NREPLIC <- reps
  ALPHA <- 0.05;  CV <- qnorm(1-ALPHA/2)
  cat(sprintf("NREPLIC=%d,BREPLIC=%d,p=%g,ALPHA=%5.3f",NREPLIC,BREPLIC,p,ALPHA),
      sep='\n',file=OUTFILE,append=TRUE)
  rxfn <- function(n) runif(n)
  ryfn <- function(n) runif(n)
  ns <- c(1e4,1e5,1e6) #1e4,1e5 and 1e6 (much slower)
  for (NUM.X0 in c(1,10,100)) {
    cat(sprintf("NUM.X0=%d",NUM.X0),sep='\n',file=OUTFILE,append=TRUE)
    tGKtots <- tBStots <- rep.int(0,length(ns))
    for (nind in 1:length(ns)) {
      n <- ns[nind]
      for (irep in 1:NREPLIC) {
        X <- rxfn(n);  Y <- ryfn(n)
        x0s <- quantile(x=X,probs=seq(from=0.25,to=0.75,length.out=NUM.X0))
        ttot <- system.time({tmp <- qme.inf(Y=Y,X=X,x0s=x0s,p=p,ALPHA=0.05,ONESIDED=0)})[1]
        tGK <- system.time({tmp <- qme.inf(Y=Y,X=X,x0s=x0s,p=p,ALPHA=0.05,ONESIDED=0,hs=tmp$hs)})[1]
        if (!NO.BS) {
          tBS <- system.time({maxhs <- apply(cbind(x0s-min(X),max(X)-x0s),1,min)
                              cubehs <- 2*1.75*tmp$hs
                              cubeh <- apply(cbind(maxhs,cubehs),1,min)
                              CI.lo <- CI.hi <- rep.int(NA,NUM.X0)
                              for (ix in 1:NUM.X0) {
                                h <- cubehs[ix]
                                iloc <- (x0s[ix]-h<=X) & (X<=x0s[ix]+h)
                                Yloc <- Y[iloc];  Xloc <- X[iloc]-x0s[ix]  #center X s.t. x0=0, as in Chaudhuri (1991)
                                if (length(Yloc)<2) {
                                  CI.lo[ix] <- CI.hi[ix] <- NA
                                  next
                                }
                                Xdes <- cbind(1,Xloc/h,Xloc^2/h^2,Xloc^3/h^3)
                                useNA <- tryCatch({
                                  rq.ret <- rq(Yloc~Xdes[,2]+Xdes[,3]+Xdes[,4],tau=p)
                                  FALSE},
                                  error=function(rqERR) { return(TRUE) }
                                )
                                if (useNA) {
                                  CI.lo[ix] <- CI.hi[ix] <- NA
                                } else {
                                  other.est <- rq.ret$coefficients[2] / h
                                  tryCatch({
                                    other.SE <- summary(rq.ret,se='boot',R=BREPLIC)$coefficients[2,2] / h
                                    CI.lo[ix] <- other.est - CV * other.SE
                                    CI.hi[ix] <- other.est + CV * other.SE 
                                  }, error=function(ERR){
                                    warning(ERR$message)
                                    CI.lo[ix] <- CI.hi[ix] <- NA
                                  })
                                }
                              }
          })[1]
        } else {
          tBS <- NA
        }
        tGKtots[nind] <- ttot+tGKtots[nind]
        tBStots[nind] <- tBS+ttot-tGK+tBStots[nind]
      } #irep loop
    } #loop over n
    cat(paste0("Bootstrap ",paste0(sprintf("& %7.1f ",tBStots/NREPLIC),collapse=''),' \\\\',collapse=''),sep='\n',file=OUTFILE,append=TRUE)
    cat(paste0("New       ",paste0(sprintf("& %7.1f ",tGKtots/NREPLIC),collapse=''),' \\\\',collapse=''),sep='\n',file=OUTFILE,append=TRUE)
  } #loop over NUM.X0
  cat(sprintf("n=%d",ns),sep=' ',file=OUTFILE,append=TRUE)
  cat("\n",file=OUTFILE,append=TRUE)
  cat(sprintf("rxfn=%s, ryfn=%s",as.character(as.expression(body(rxfn))),as.character(as.expression(body(ryfn)))),file=OUTFILE,append=TRUE)
}

#EOF
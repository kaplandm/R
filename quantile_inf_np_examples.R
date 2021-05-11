source("quantile_inf_np.R")

# Examples of using quantile.inf.np(), with both simulated and observed data

##################
# From paper     #
##################
# Quantile Engel curves from UK LCFS 2012, http://esds.ac.uk/doi/?sn=7472#2

SAVE.FLAG <- TRUE
SUBMAX <- Inf #max subsample size for CI comparison

require(foreign) #to load .dta files
require(KernSmooth) #for dpill()
require(fields) #for qsreg()
source("quantile_inf_np.R")
raw <- vector('list',12)
for (i in 1:12) {
  raw[[i]] <- read.dta(sprintf("20%02d_dvhh_ukanon.dta",i),convert.factors=FALSE)
}
moneycols <- c('P630tp','P539t','P601t','P607t','P611t',
               'P604t','P609t','P603t','P605t') 
moneycols <- c('P630tp','P537t','P538t','P539t','P541t')
keepcols <- c('A062','Gorx','A121','Year',
              'A023','A024','A025','A026','A027',
              'A033','A034','A035','A036','A037',
              moneycols)
for (i in 10:12) {
  names(raw[[i]])[names(raw[[i]]) %in% keepcols]  <- 
    tolower(names(raw[[i]])[names(raw[[i]]) %in% keepcols])
}
keepcols <- tolower(keepcols);  moneycols <- tolower(moneycols)
#
dats <- vector('list',12)
for (i in 1:12) {
  dats[[i]] <- raw[[i]][,keepcols]
}
#
CPI <- c(94.2,95.4,96.7,98.0,100,102.3,104.7,108.5,110.8,114.5,119.6,123.0) #2001-2012; 2005=100
for (i in 1:11) {
  dats[[i]][,moneycols] <- dats[[i]][,moneycols]*CPI[length(CPI)]/CPI[i]
}
#
dat <- dats[[1]]
for (i in 2:12) { dat <- rbind(dat,dats[[i]]) }
#
f.under60 <- (dat$a025+dat$a026+dat$a027+dat$a035+dat$a036+dat$a037==0)
f.under65 <- (dat$a026+dat$a027+dat$a036+dat$a037==0)
f.1m <- (dat$a062 %in% c(1,3,5))
f.1w <- (dat$a062 %in% c(2,4,6))
f.1a <- f.1m | f.1w
f.1m1w <- dat$a062 %in% c(7,9,11,13)
f.2a <- (dat$a062 %in% c(7:17))
f.0c <- (dat$a062 %in% c(1:2,7:8,18,23,26))
f.1c <- (dat$a062 %in% c(3,4,9:10,19,24))
f.2c <- (dat$a062 %in% c(5,6,11:12,20))
f.2plusc <- (dat$a062 %in% c(5,6,11:17,20:22,25,27))
f.London <- (dat$gorx==7)
f.SE <- (dat$gorx==8)
f.E <- (dat$gorx==6)
f.SW <- (dat$gorx==9)
f.Mid <- (dat$gorx %in% 4:5)
f.England <- (dat$gorx %in% 1:9)
f.Wales <- (dat$gorx %in% 10)
f.Scotland <- (dat$gorx %in% 11)
f.NI <- (dat$gorx %in% 12)
f.renter <- (dat$a121 %in% c(1:4,8))
#
filter <- (dat$p630tp!=0)
filter <- filter & f.1m1w & (f.London | f.SE)
#filter <- filter & f.0c
sum(filter)

cats <- rbind(c('p538t','p537t','p539t','p541t','other'),
              c('Food','Fuel, light, and power','Clothing/footwear','Alcohol','Other'),
              c('food','fuel','clothing','alc','other'))
dat[['other']] <- dat$p630tp
for (catind in 1:4) {
  dat[['other']] <- dat[['other']] - dat[[cats[1,catind]]]
}
expf <- dat$p630tp[filter]
lnexpf <- log(expf);  lnexpf2 <- lnexpf^2
lnexpf.orig <- lnexpf
x0s <- quantile(lnexpf,p=1:9/10)
ps <- c(0.5,0.75,0.9) #1:3/4
ALPHA <- 0.10

for (catind in 1:5) {
  set.seed(112358)
  w <- dat[[cats[1,catind]]][filter]/expf
  w <- ifelse(w<0,0,w)
  w.orig <- w
  #
  # source("quantile_inf.R")
  tmp50 <- quantile.inf(X=w,p=0.50,ALPHA=0.01)
  tmp75 <- quantile.inf(X=w,p=0.75,ALPHA=0.01)
  tmp90 <- quantile.inf(X=w,p=0.90,ALPHA=0.01)
  cat(sprintf("L-stat & %-8s & (%6.4f,%6.4f) & (%6.4f,%6.4f) & (%6.4f,%6.4f) \\\\\n", 
              cats[3,catind], tmp50$CI[1], tmp50$CI[2], tmp75$CI[1], tmp75$CI[2], tmp90$CI[1], tmp90$CI[2]),
      file=ifelse(SAVE.FLAG,"quantile_inf_np_ex_uncond.txt",""),append=TRUE)
  tmp50 <- quantile.inf(X=w,p=0.50,ALPHA=0.01,SINGLE.CALIB=TRUE)
  tmp75 <- quantile.inf(X=w,p=0.75,ALPHA=0.01,SINGLE.CALIB=TRUE)
  tmp90 <- quantile.inf(X=w,p=0.90,ALPHA=0.01,SINGLE.CALIB=TRUE)
  cat(sprintf("Calib  & %-8s & (%6.4f,%6.4f) & (%6.4f,%6.4f) & (%6.4f,%6.4f) \\\\\n", 
              cats[3,catind], tmp50$CI[1], tmp50$CI[2], tmp75$CI[1], tmp75$CI[2], tmp90$CI[1], tmp90$CI[2]),
      file=ifelse(SAVE.FLAG,"quantile_inf_np_ex_uncond.txt",""),append=TRUE)
  #
  if (FALSE) {
    par(family='serif',mar=c(5.0,6.0,6.0,2.1))
    plot(lnexpf,ifelse(w==0,runif(length(w),-0.01,0),w),main=cats[2,catind],
         pch=16,cex=0.1,ylim=c(0,1),xlim=range(x0s)*c(0.9,1.1),
         cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
         xlab="Log total expenditure",ylab="Budget share")
    lines(c(-99,99),c(0,0),col=1,lwd=1)
  }
  if (FALSE) {
    subind <- sample(1:length(lnexpf),2000,T)
    if (SAVE.FLAG) pdf(file=sprintf("quantile_inf_np_ex_%s_subscatter.pdf",cats[3,catind])) #else draw to current device
    par(family='serif',mar=c(5.0,6.0,6.0,2.1))
    plot(lnexpf[subind],ifelse(w[subind]==0,runif(length(w[subind]),-0.01,0),w[subind]),
         main=sprintf("%s:\nsubsample of %d",cats[2,catind],length(subind)),
         pch=16,cex=0.1,ylim=c(0,1),xlim=range(x0s)*c(0.9,1.1),
         cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
         xlab="Log total expenditure",ylab="Log budget share")
    lines(c(-99,99),c(0,0),col=1,lwd=1)
    if (SAVE.FLAG) dev.off()
  }
  CI.pt.los <- CI.pt.his <- CI.jt.los <- CI.jt.his <- matrix(NA,length(x0s),length(ps))
  rqs <- vector('list',length(ps))
  rqs.lin <- vector('list',length(ps))
  for (pind in 1:length(ps)) {
    ret <- quantile.inf.np(Y=w,X=lnexpf,x0s=x0s,p=ps[pind],ALPHA=ALPHA,JOINT.FLAG=TRUE,ONESIDED=0,METHOD.TYPE='single')
    CI.pt.los[,pind] <- ret[[1]]$CI.pointwise.lows
    CI.jt.los[,pind] <- ret[[1]]$CI.joint.lows
    CI.pt.his[,pind] <- ret[[1]]$CI.pointwise.highs
    CI.jt.his[,pind] <- ret[[1]]$CI.joint.highs
    rqs[[pind]] <- rq(w~lnexpf+lnexpf2,tau=ps[pind])
    rqs.lin[[pind]] <- rq(w~lnexpf,tau=ps[pind])
  }
  xx <- seq(from=min(x0s)*0.9,to=max(x0s)*1.1,by=0.001)
  if (FALSE) {
    if (SAVE.FLAG) pdf(file=sprintf("quantile_inf_np_ex_pt_%s.pdf",cats[3,catind])) #else draw to current device
    par(family='serif',mar=c(5.0,6.0,6.0,2.1))
    tmp <- sort(c(CI.pt.los,CI.pt.his))
    plot(x0s[c(1,length(x0s))],quantile(tmp,p=c(0.0,1.0)),type='n',
         cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
         main=sprintf("Pointwise: %s",cats[2,catind]), 
         xlab="log expenditure", ylab="budget share")
    for (pind in 1:length(ps)) {
      lines(x0s,CI.pt.los[,pind],lwd=2,lty=pind,pch=pind,col=1,type='o',cex=1.5)
      lines(x0s,CI.pt.his[,pind],lwd=2,lty=pind,pch=pind,col=1,type='o',cex=1.5)
      yy <- predict.rq(rqs[[pind]],newdata=data.frame(lnexpf=xx,lnexpf2=xx^2))
      lines(xx,yy,lwd=2,lty=pind,pch=NA,col=1,type='o',cex=1.5)
    }
    legend("top",sprintf("%4.2f-quantile",ps),inset=c(0.0,0.0), bty='n',
           col=1, pch=1:length(ps), lty=1:length(ps), lwd=2, cex=1.3)
    if (SAVE.FLAG) dev.off()
  }
  #
  for (i in 1:2) {
    if (SAVE.FLAG) {#else draw to current device
      if (i==1) {
        pdf(file=sprintf("quantile_inf_np_ex_jt_%s.pdf",cats[3,catind])) 
      } else {
        pdf(file=sprintf("quantile_inf_np_ex_jt_lin_%s.pdf",cats[3,catind]))
      }
    }
    par(family='serif',mar=c(5.0,6.0,6.0,2.1))
    tmp <- sort(c(CI.jt.los,CI.jt.his))
    plot(x0s[c(1,length(x0s))],quantile(tmp,p=c(0.0,1.0)),type='n',
         cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
         main=sprintf("Joint: %s",cats[2,catind]), 
         xlab="log expenditure", ylab="budget share")
    for (pind in 1:length(ps)) {
      lines(x0s,CI.jt.los[,pind],lwd=2,lty=pind,pch=pind,col=1,type='o',cex=1.5)
      lines(x0s,CI.jt.his[,pind],lwd=2,lty=pind,pch=pind,col=1,type='o',cex=1.5)
      if (i==1) {
        yy <- predict.rq(rqs[[pind]],newdata=data.frame(lnexpf=xx,lnexpf2=xx^2))
      } else {
        yy <- predict.rq(rqs.lin[[pind]],newdata=data.frame(lnexpf=xx,lnexpf2=xx^2))
      }
      lines(xx,yy,lwd=2,lty=pind,pch=NA,col=1,type='o',cex=1.5)
    }
    legend("top",sprintf("%4.2f-quantile",ps),inset=c(0.0,0.0), bty='n',
           col=1, pch=1:length(ps), lty=1:length(ps), lwd=2, cex=1.3)
    if (SAVE.FLAG) dev.off()
  }
  
  ### Compare w/ other method, and include nonparametric estimates
  SUBSIZE <- min(SUBMAX,length(w))
  rsamp <- sample(1:length(w), SUBSIZE, replace=FALSE)
  w <- w.orig[rsamp];  lnexpf <- lnexpf.orig[rsamp]  #restored below
  if (SAVE.FLAG) pdf(file=sprintf("quantile_inf_np_ex_jt_%s_comparison_sub%d.pdf",cats[3,catind],SUBSIZE))
  RUNSIM <- FALSE;  source("Fan_Liu_2015_CI_code.R")
  FL.CIs.jt <- matrix(NA,nrow=2*length(ps),ncol=length(x0s))
  qsreg.est <- matrix(NA,nrow=length(ps),ncol=length(xx))
  rqss.est <- matrix(NA,nrow=length(ps),ncol=length(xx))
  #
  n <- length(w)
  FLdelta <- 1/5+1/20;  intKsq <- 5/7;  mu2K <- 1/7;  intKpsq <- 15/7
  UNIFcv <- (log(2)-log(abs(log(1-ALPHA))))/sqrt(2*FLdelta*log(n)) + 
    sqrt(2*FLdelta*log(n)) + log(intKpsq/(4*pi*intKsq)) / sqrt(2*FLdelta*log(n)) #FL eqn (20)
  if (UNIFcv<qnorm(1-(ALPHA/length(x0s))/2)) stop("Should use FL uniform band (tighter than Bonferroni)")
  #
  for (pind in 1:length(ps)) {
    p <- ps[pind]
    for (ix in 1:length(x0s)) {
      FL.CIs.jt[2*pind-1:0,ix] <- FL.CI.fn(lnexpf,w,x0s[ix],p,ALPHA/length(x0s),UNIF=FALSE,ORIGINAL=TRUE)
    }
    tmp <- qsreg(lnexpf,w,alpha=p)
    qsreg.est[pind,] <- predict(tmp,derivative=0,x=xx)
    g <- function(lam,y,x) AIC(rqss(y ~ qss(x, lambda=lam), tau=p),k = -1)
    tmp <- sample(1:length(w),min(1000,length(w)),replace=FALSE)
    lamstar <- optimize(g, interval = c(0.001, 0.5), x=lnexpf[tmp], y=w[tmp])
    tmp <- rqss(w~qss(lnexpf, lambda=lamstar$min), tau=p)
    rqss.est[pind,] <- predict.rqss(tmp,newdata=data.frame(lnexpf=xx))
    #
    ret <- quantile.inf.np(Y=w,X=lnexpf,x0s=x0s,p=ps[pind],ALPHA=ALPHA,JOINT.FLAG=TRUE,ONESIDED=0,METHOD.TYPE='single')
    CI.pt.los[,pind] <- ret[[1]]$CI.pointwise.lows
    CI.jt.los[,pind] <- ret[[1]]$CI.joint.lows
    CI.pt.his[,pind] <- ret[[1]]$CI.pointwise.highs
    CI.jt.his[,pind] <- ret[[1]]$CI.joint.highs
  }
  #
  par(family='serif',mar=c(5.0,6.0,6.0,2.1))
  tmp <- sort(c(FL.CIs.jt,CI.jt.los,CI.jt.his))
  plot(x0s[c(1,length(x0s))],quantile(tmp,p=c(0.0,1.0)),type='n',
       cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
       main=sprintf("Joint: %s (subsample of %d)",cats[2,catind],SUBSIZE), 
       xlab="log expenditure", ylab="budget share")
  COLS <- c("#A6BDDB","#014636")
  for (pind in c(1,3)) { #length(ps)
    lines(xx,qsreg.est[pind,],     lwd=1,lty=1,pch=NA,col=1,type='o',cex=1) #qsreg.est or rqss.est
    lines(x0s,FL.CIs.jt[2*pind-1,],lwd=2,lty=2,pch=1,col=COLS[1],type='o',cex=1.5)
    lines(x0s,FL.CIs.jt[2*pind,],  lwd=2,lty=2,pch=1,col=COLS[1],type='o',cex=1.5)
    lines(x0s,CI.jt.los[,pind],lwd=2,lty=1,pch=4,col=COLS[2],type='o',cex=1.5)
    lines(x0s,CI.jt.his[,pind],lwd=2,lty=1,pch=4,col=COLS[2],type='o',cex=1.5)
  }
  legpos <- "topright"
  if (catind>=4) legpos <- "topleft"
  legend(legpos,c("estimate","Fan-Liu","L-stat"),inset=c(0.0,0.0), bty='n',
         col=c("#000000",COLS[1],COLS[2]), pch=c(NA,1,4), lty=c(1,2,1), lwd=c(1,2,2), cex=1.3)
  if (SAVE.FLAG) dev.off()
  #
  lnexpf <- lnexpf.orig;  w <- w.orig
}
cat(sum(filter))


stop("End of UK LFCS examples.")


##################
# Simulated data #
##################
set.seed(112358)
NUM.X0 <- 47
TRIM.X0 <- 0.04 #trim off top and bottom TRIM portion of X0 points
#
ONESIDED <- 0 #0:two-sided; 1/-1:upper/lower one-sided
ALPHA <- 0.05
n <- 1000 #rqss vignette: n=400
p <- 0.5 #quantile, b/w 0 and 1
#
x0.ps <- TRIM.X0 + (1-2*TRIM.X0)*c(0:(NUM.X0-1))/(NUM.X0-1) #X points of interest, as quantiles of X dist'n
Findex <- 010
#From rqss vignette (from Ruppert, Wand, and Carroll (2003, \S17.5.1))
X <- runif(n,0,1)
x0s <- qunif(x0.ps,0,1)
Y.fn <- function(X) (X*(1-X))^(1/2)*sin(2*pi*(1+2^(-7/5))/(X+2^(-7/5)))
sig0 <- 0.2
#Homoskedastic if Findex is even, heteroskedastic if odd
if (((Findex-10*floor(Findex/10)) %% 2) == 0) {
    sig.fn <- function(X) rep.int(sig0,length(X)) #homoskedastic
} else {
    sig.fn <- function(X) sig0*(1+X) #heteroskedastic
}
#Different distribution shapes for 010/011 vs. 012/013 vs. ...
if (floor((Findex-10*floor(Findex/10))/2)==0) {
    U <- rnorm(n,0,1) - qnorm(p,0,1) #Gaussian (010/011)
    # f.Y.fn <- function(p,Xs) dnorm(qnorm(p,0,sig.fn(Xs)),0,sig.fn(Xs))
} else {stop("Oops.")}
Y <- Y.fn(X)+sig.fn(X)*U
y0s <- Y.fn(x0s) #no add'l effect of heterosk. since p-quantile of U is zero
iqr0s <- sig.fn(x0s)*(qnorm(0.75,0,1)-qnorm(0.25,0,1))

# EXAMPLE S.1a: NO DISCRETE X, single quantile
source("quantile_inf_np.R")
system.time(
rets <- quantile.inf.np(Y=Y,X=X,x0s=x0s, p=p,ALPHA=ALPHA,JOINT.FLAG=TRUE,ONESIDED=0,KERNEL.TYPE='gaussian')
)
plot.quantile.inf.np(rets,plot.data=F,CI.interpolate=T)

# EXAMPLE S.1b: NO DISCRETE X, joint upper & lower quartile
source("quantile_inf_np.R")
x0sj <- x0s[seq(from=2,to=length(x0s)-1,by=3)]
system.time(
rets <- quantile.inf.np(Y=Y,X=X,x0s=x0sj, p=c(0.75,0.25),ALPHA=ALPHA,JOINT.FLAG=TRUE,METHOD.TYPE='joint', BETA.BLOCK.SIZE=2500, BETA.BLOCKS=1)
)
#quartz()
par(family="serif")
plot(X,Y,pch=20,cex=0.25,main='Pointwise (by x0) L-stat 95% CIs for\nconditional joint upper and lower quartiles',mgp=c(2.1,0.5,0), cex.main=1.8, cex.axis=1.5, cex.lab=2)
points(c(x0sj,x0sj),c(rets[[1]]$CI.pointwise.lows[,1],rets[[1]]$CI.pointwise.highs[,1]),pch=2,col=4,lwd=2,cex=0.6)
points(c(x0sj,x0sj),c(rets[[1]]$CI.pointwise.lows[,2],rets[[1]]$CI.pointwise.highs[,2]),pch=6,col=3,lwd=2,cex=0.6)
legend('bottomright', inset=0.01, legend=c('Data','Upper quartile CI','Lower quartile CI'), pch=c(20,2,6), col=c(1,4,3), lwd=2, lty=NA, cex=1.5)

# EXAMPLE S.1c: NO DISCRETE X, lincom (IQR) inference
source("quantile_inf_np.R")
x0sl <- x0s[seq(from=1,to=length(x0s),by=5)]
iqr0sl <- iqr0s[seq(from=1,to=length(iqr0s),by=5)]
system.time(
retsl <- quantile.inf.np(Y=Y,X=X,x0s=x0sl, p=c(0.75,0.25),ALPHA=ALPHA,JOINT.FLAG=TRUE,METHOD.TYPE='lincom',PSI=c(1,-1), BETA.BLOCK.SIZE=6000, BETA.BLOCKS=1)
)
#quartz()
par(family="serif")
plot(x0sl,iqr0sl,pch=20,cex=0.25,main='Pointwise L-stat 95% CIs for\nconditional interquartile range (IQR)',mgp=c(2.1,0.5,0), cex.main=1.8, cex.axis=1.5, cex.lab=2, col=1, type='p', ylim=c(0,1), xlab='X',ylab='IQR')
points(c(x0sl,x0sl),c(retsl[[1]]$CI.pointwise.lows,retsl[[1]]$CI.pointwise.highs),pch=3,col=3,lwd=2,cex=0.6)
legend('top', inset=0.01, legend=c('True IQR','L-stat CI'), pch=c(20,3), col=c(1,3), lwd=2, lty=NA, cex=1.5)

# EXAMPLE S.2: SINGLE BINARY DUMMY
#Data are sorted so that dummy=1 gets generally higher values.
X2 <- cbind(X,as.integer(Y>Y.fn(X))) #X2 <- cbind(X,c(rep.int(0,n/2),rep.int(1,n/2)))
source("quantile_inf_np.R")
rets2 <- quantile.inf.np(Y=Y,X=X2,x0s=x0s, p=p,ALPHA=ALPHA,JOINT.FLAG=TRUE,ONESIDED=0,KERNEL.TYPE='gaussian')
plot.quantile.inf.np(rets2,plot.data=F,CI.interpolate=T,x.axis.title="Simulated cts X")

# EXAMPLE S.3: SECOND BINARY DUMMY.
#Data are simply replicated so that the graphs for (0,0) and (0,1) should be identical, and (1,0) and (1,1) should also be identical.
Y3 <- c(Y,Y)
X3 <- cbind(rbind(X2,X2),c(rep.int(0,n),rep.int(1,n)))
source("quantile_inf_np.R")
rets3 <- quantile.inf.np(Y=Y3,X=X3,x0s=x0s, p=p,ALPHA=ALPHA,JOINT.FLAG=TRUE,ONESIDED=0,KERNEL.TYPE='gaussian')
plot.quantile.inf.np(rets3,plot.data=F,CI.interpolate=T,x.axis.title="Simulated cts X")

# EXAMPLE S.4: TREATMENT, NO DISCRETE
source("quantile_inf_np.R")
X2 <- cbind(X,as.integer(Y>Y.fn(X))) #X2 <- cbind(X,c(rep.int(0,n/2),rep.int(1,n/2)))
rets4 <- quantile.inf.np(Y=Y,X=X2,x0s=x0s, p=p,ALPHA=ALPHA,JOINT.FLAG=T,ONESIDED=0,METHOD.TYPE='qte')
plot.quantile.inf.np(rets4,plot.data=T,CI.interpolate=T,CI.type='both', title='L-stat 95% CI for\nconditional median treatment effect')
# Plot w/ data on secondary y-axis
#quartz()
par(family='serif', mar=c(5,4,4,5)+.1)
CIjls <- rets4[[1]]$CI.joint.lows
CIjhs <- rets4[[1]]$CI.joint.highs
plot(c(x0s,x0s),c(rets4[[1]]$CI.joint.lows,rets4[[1]]$CI.joint.highs), type="n",main='L-stat 95% CI (interpolated for visual ease) for\nconditional median treatment effect', mgp=c(2.1,0.5,0), cex.main=1.8, cex.axis=1.5, cex.lab=2, xlab='X', ylab='QTE', ylim=c(min(CIjls[abs(CIjls)!=Inf],na.rm=T)-0.7,max(CIjhs[abs(CIjhs)!=Inf],na.rm=T)+0.0))
points(type='o',x=x0s,y=rets4[[1]]$CI.pointwise.lows, pch=3,col=3, lwd=2,lty=3)
points(type='o',x=x0s,y=rets4[[1]]$CI.pointwise.highs, pch=3,col=3, lwd=2,lty=3)
points(type='o',x=x0s,y=rets4[[1]]$CI.joint.lows, pch=4,col=4, lwd=2,lty=3)
points(type='o',x=x0s,y=rets4[[1]]$CI.joint.highs, pch=4,col=4, lwd=2,lty=3)
par(family='serif',new=TRUE)
plot(type='n',x=X2[,1],y=Y, col=1, xaxt="n",yaxt="n",xlab="",ylab="", ylim=c(-1.0,1.5))
points(X2[X2[,2]==1,1],Y[X2[,2]==1], col=1, pch=84, cex=0.4)
points(X2[X2[,2]==0,1],Y[X2[,2]==0], col=1, pch=67, cex=0.4)
axis(4, cex.axis=1.5)
mtext("Y",side=4,line=3, cex=2)
legend('bottomright',legend=c('Data (treated)','Data (control)','Pointwise','Joint'), pch=c(84,67,3,4), col=c(1,1,3,4), lwd=2, lty=NA, inset=0.01, cex=1.5, bg='white', y.intersp=0.9)

# EXAMPLE S.5: TREATMENT, ONE DISCRETE
source("quantile_inf_np.R")
X2 <- cbind(X,as.integer(Y>Y.fn(X))) #X2 <- cbind(X,c(rep.int(0,n/2),rep.int(1,n/2)))
X3 <- cbind(rbind(X2,X2),c(rep.int(0,n),rep.int(1,n)))
rets5 <- quantile.inf.np(Y=c(Y,-Y),X=X3,x0s=x0s, p=p,ALPHA=ALPHA,JOINT.FLAG=F,ONESIDED=0,METHOD.TYPE='qte')
plot.quantile.inf.np(rets5,plot.data=T,CI.interpolate=T,title='L-stat 95% CI for\nconditional median treatment effect',y.axis.title='Y (for data); QTE (for CI)')







##################
# Empirical data #
##################

#
# Example E.0: value-at-risk (VaR)
#
source("quantile_inf_np.R")
#Graphing fn for both 2a and 2b
graph.VaR <- function(X.cts.all,Y.cts.all,rets,rets.VaR.unconditional,X.cts.varname="X",plot.data=FALSE,ylim=c(-.08,.08),xlim=c(-.05,.05),legpos="topright",datacex=0.1,lwd=1,datapch=1) {
  pt.lo <- vector("list",length(rets));  pt.lo.ss <- vector("list",length(rets))
  pt.hi <- vector("list",length(rets));  pt.hi.ss <- vector("list",length(rets))
  jt.lo <- vector("list",length(rets));  jt.lo.ss <- vector("list",length(rets))
  jt.hi <- vector("list",length(rets));  jt.hi.ss <- vector("list",length(rets))
  #
  for (tmpi in 1:length(rets)) {
    pt.lo[[tmpi]] <- rets[[tmpi]]$CI.pointwise.lows
    pt.hi[[tmpi]] <- rets[[tmpi]]$CI.pointwise.highs
    jt.lo[[tmpi]] <- rets[[tmpi]]$CI.joint.lows
    jt.hi[[tmpi]] <- rets[[tmpi]]$CI.joint.highs
  }
  X.cts <- X.cts.all;  Y.cts <- Y
  ALPHA <- rets[[1]]$ALPHA;  p <- rets[[1]]$p
  x0s <- vector("list",length(rets))
  for (tmpi in 1:length(rets)) {x0s[[tmpi]]<-rets[[tmpi]]$x0s}
  #Interpolated bands
  for (tmpi in 1:length(rets)) {
    pt.lo.ss[[tmpi]] <- smooth.spline(x0s[[tmpi]],pt.lo[[tmpi]], all.knots=TRUE,penalty=0)
    pt.hi.ss[[tmpi]] <- smooth.spline(x0s[[tmpi]],pt.hi[[tmpi]], all.knots=TRUE,penalty=0)
    jt.lo.ss[[tmpi]] <- smooth.spline(x0s[[tmpi]],jt.lo[[tmpi]], all.knots=TRUE,penalty=0)
    jt.hi.ss[[tmpi]] <- smooth.spline(x0s[[tmpi]],jt.hi[[tmpi]], all.knots=TRUE,penalty=0)
  }
  #
  #X prep
  xr <- range(x0s)
  xx <- seq(from=xr[1],to=xr[2],by=(xr[2]-xr[1])/length(Y.cts.all))

 #Draw graph
 #quartz()
 par(family="serif")
 plot(x=X.cts.all,y=Y.cts.all,type="n",pch=1,col=1, main=sprintf("L-stat %g%% CI for %g-quantile",(1-ALPHA)*100,p), xlab=X.cts.varname, ylab="Daily return, ln(P[t]/P[t-1])", cex.main=2.0, mgp=c(2.1,0.5,0), cex.axis=1.5, cex.lab=2, ylim=ylim, xlim=xlim) 
 #
 #Data: small dots
 if (plot.data) {
   points(X.cts.all,Y.cts.all,type='p',pch=datapch,col=1,cex=datacex)
   legtxt<-c('Data','Unconditional'); legcol<-c(1,2); legpch<-c(datapch,NA); leglty<-c(0,2)
 } else {
   legtxt<-c('Unconditional'); legcol<-c(2); legpch<-c(NA); leglty<-c(2)
 }
 #Unconditional: horizontal dotted lines
 lines(c(-1,1),c(rets.VaR.unconditional$CI$lo,rets.VaR.unconditional$CI$lo),type='l',col=2,lty=2,lwd=5)
 lines(c(-1,1),c(rets.VaR.unconditional$CI$hi,rets.VaR.unconditional$CI$hi),type='l',col=2,lty=2,lwd=5)
 #Conditional
 for (tmpi in 1:length(rets)) {
   points(x0s[[tmpi]],pt.lo[[tmpi]],type='p',col=tmpi+2,pch=tmpi+2,cex=1.5)
   points(x0s[[tmpi]],pt.hi[[tmpi]],type='p',col=tmpi+2,pch=tmpi+2,cex=1.5)
    lines(predict(pt.lo.ss[[tmpi]],xx),type='l',col=tmpi+2,lwd=lwd,lty=3)
    lines(predict(pt.hi.ss[[tmpi]],xx),type='l',col=tmpi+2,lwd=lwd,lty=3)
   if (length(rets)==1) {
     points(x0s[[tmpi]],jt.lo[[tmpi]],type='p',col=tmpi+3,pch=tmpi+3,cex=1.5)
     points(x0s[[tmpi]],jt.hi[[tmpi]],type='p',col=tmpi+3,pch=tmpi+3,cex=1.5)
     lines(predict(jt.lo.ss[[tmpi]],xx),type='l',col=tmpi+3,lwd=lwd,lty=3)
     lines(predict(jt.hi.ss[[tmpi]],xx),type='l',col=tmpi+3,lwd=lwd,lty=3)
   }
   legtxt<-c(legtxt,paste0("Pointwise ",as.character(rets[[tmpi]]$X.discrete)))
   legcol<-c(legcol,tmpi+2)
   legpch<-c(legpch,tmpi+2); leglty<-c(leglty,3)
 }
 #
 if (length(rets)==1) {
  if (plot.data) {
    legend(legpos, legend=c('Data','Unconditional','Pointwise','Joint'), inset=0.01, col=c(1,2,3,4), pch=c(datapch,NA,3,4), lty=c(0,2,3,3), lwd=lwd/1.5, cex=1.5, y.intersp=0.9, bg='white')
  } else {
    legend(legpos, legend=c('Unconditional','Pointwise','Joint'), inset=0.01, col=c(2,3,4), pch=c(NA,3,4), lty=c(2,3,3), lwd=lwd/1.5, cex=1.5, y.intersp=0.9, bg='white')
  }
 } else {
   legend(legpos, legend=legtxt, inset=0.01, col=legcol, pch=legpch, lty=leglty, lwd=lwd/1.5, cex=1.5, y.intersp=0.9, bg='white')
 }
}





#
# Example E.0a: VaR
#
#Chernozhukov and Fernandez-Val (2011)
#Their data: http://restud.oxfordjournals.org/content/78/2/559/suppl/DC1
#But: they drop observations when return is zero (why?); see Example E.2b for including zeroes
#September 1986 to November 1998 for rows 1:2527
#Y (and X3) is daily returns computed as ln(P_t/P_{t-1}) in terms of daily prices P_t
#Data input from Chernozhukov & Fernandez-Val (2011),
# available @ https://drive.google.com/file/d/0B-_LUSJVBv20aWRweldoWTN0eU0/view?usp=sharing
VaR.data.raw <- read.table("CF2011_VaR_data.txt",colClasses='numeric')
Y       <- VaR.data.raw[1528:2527,5] #VaR.data.raw[1528:2527, 5]
Xt      <- VaR.data.raw[1528:2527,1:4] #VaR.data.raw[1528:2527, 1:4]
X1.pos      <- (Xt[ ,2] >= 0) * Xt[ ,2]
X1.neg      <- -(Xt[ ,2] <= 0) * Xt[ ,2]
X2.pos      <- (Xt[ ,3] >= 0) * Xt[ ,3]
X2.neg      <- -(Xt[ ,3] <= 0) * Xt[ ,3]
X3.pos      <- (Xt[ ,4] >= 0) * Xt[ ,4]
X3.neg      <- -(Xt[ ,4] <= 0) * Xt[ ,4]
X1      <- Xt[ ,2]
X2      <- Xt[ ,3]
#<--end data input from C & F-V (2011)
p <- 0.05;  ALPHA <- 0.1
#Unconditional VaR
rets.VaR.unconditional <- quantile.inf(X=Y,p=p,ALPHA=ALPHA)
#Conditional VaR
#Pick a cts X to use (of 3 options)
X.cts.ind <- 2;  X.cts.varname <- "Yesterday's oil spot price return"; x0s=seq(from=-0.05,to=0.05,by=0.005)
X.cts.ind <- 3;  X.cts.varname <- "Yesterday's DJIA return"; x0s=seq(from=-0.015,to=0.02,by=0.005)
X.cts.ind <- 4;  X.cts.varname <- "Yesterday's OXY return"; x0s=seq(from=-0.025,to=0.025,by=0.005)
#
rets.VaR <- quantile.inf.np(Y=Y,X=Xt[,X.cts.ind],x0s=x0s,p=p,ALPHA=ALPHA)
tmph <- rets.VaR[[1]]$bandwidth.joint
tmpN <- rets.VaR[[1]]$effective.N.joint
rets.VaR <- quantile.inf.np(Y=Y,X=Xt[,X.cts.ind],x0s=x0s,p=p,ALPHA=ALPHA,hs.joint=ifelse(tmpN<20,0.0025,tmph))
#
#Graphing
Y.cts.all <- Y
X.cts.all <- Xt[,X.cts.ind]
#
graph.VaR(X.cts.all,Y.cts.all,rets.VaR,rets.VaR.unconditional,X.cts.varname)



#
# Example E.0b: VaR
#
#Google finance: http://www.google.com/finance/historical?cid=26001&startdate=Aug+29%2C+1986&enddate=Aug+31%2C+2012&num=200
#.csv available @ https://drive.google.com/file/d/0B-_LUSJVBv20anFLdTZpQlJfbms/view?usp=sharing
VaR.data.raw <- read.csv("OXY_daily_data_no_holidays.csv")
Y <- VaR.data.raw$DailyLnRet
X <- VaR.data.raw$LagLnRet
X.year <- VaR.data.raw$Year
X.per <- (X.year>=1990 & X.year<=1999) + 2*(X.year>=2000 & X.year<=2009)
X.per.str <- sprintf("%d",X.per)
X.per.str[X.per==1] <- "1990s"
X.per.str[X.per==2] <- "2000s"
p <- 0.05;  ALPHA <- 0.1
#Unconditional VaR
rets.VaR.unconditional <- quantile.inf(X=Y,p=p,ALPHA=ALPHA)
#
#Conditional VaR w/ cts X only
rets.VaR <- quantile.inf.np(Y=Y,X=X,x0s=quantile(X,seq(from=0.16-min(p,1-p),to=0.89+min(p,1-p),by=0.02)),p=p,ALPHA=ALPHA)
graph.VaR(X.cts.all=X,Y.cts.all=Y,rets=rets.VaR,rets.VaR.unconditional,X.cts.varname="Yesterday's return",ylim=c(-0.05,0.05),xlim=c(-0.03,.03),legpos="topleft",lwd=4,datacex=0.4,datapch=16)
#
#Conditional VaR w/ cts X and by year
filter <- as.logical(X.per!=0)
source("quantile_inf_np.R")
rets.VaR.yearly <- quantile.inf.np(Y=Y[filter],X=cbind(X[filter],X.per[filter]),x0s=quantile(X[filter],seq(from=0.16-min(p,1-p),to=0.87+min(p,1-p),by=0.05)),p=p,ALPHA=ALPHA)
rets.VaR.yearly[[1]]$X.discrete <- ifelse(rets.VaR.yearly[[1]]$X.discrete==1,"1990s","2000s")
rets.VaR.yearly[[2]]$X.discrete <- ifelse(rets.VaR.yearly[[2]]$X.discrete==1,"1990s","2000s")
graph.VaR(X,Y,rets.VaR.yearly,rets.VaR.unconditional,"Yesterday's return",ylim=c(-0.05,0.05),xlim=c(-0.03,0.03),legpos="topleft",lwd=4,datacex=0.5,datapch=16,plot.data=TRUE)
plot2.quantile.inf.np(rets.VaR.yearly,plot.data=FALSE,CI.type='joint',x.axis.title="Yesterday's return",y.axis.title="Daily return, ln(P[t]/P[t-1])")








#
# Example E.1: motorcycle crash data
#
library(MASS)
source("quantile_inf_np.R")
rets.mcycle <- quantile.inf.np(Y=mcycle[,2],X=mcycle[,1],x0s=5:50,p=0.5,ALPHA=0.1)
plot.quantile.inf.np(rets.mcycle,plot.data=FALSE,CI.interpolate=T,x.axis.title="Time (ms after impact)",y.axis.title="Acceleration (g)")
#points(mcycle[,1],mcycle[,2],pch=16,cex=0.1)



#
# Example E.2: birthweight: smokers vs. non-smokers (filtered on <=HS edu, married, white)
#
#Load data per Chernozhukov & Fernandez-Val (2011) code
#Data: http://restud.oxfordjournals.org/content/78/2/559/suppl/DC1
library(foreign) #for .dta input
nat <- read.dta("Birthweights_data.dta",convert.factors=FALSE)
nat$educ <- nat$ed_hs + 2*nat$ed_smcol + 3*nat$ed_col
nat <- nat[order(-nat$smoke),]
attach(nat) #make variables accessible by just using names
#weight(=depvar,in kg to start),id,birmon(=6,June only),married(0/1),{tri1,tri2,tri3,novisit},m_wtgain(lbs,integer),{ed_hs,ed_smcol,ed_col},mom_age,black(0/1),boy(0/1),smoke(0/1),cigsper(integer but grouped at 5's)
weight.g          <- 1000*weight; #was in kg(?), now g
m.wt.gain.lbs.net <- m_wtgain - 2.20462*weight
#
source("quantile_inf_np.R")
ALPHA <- 0.1;  p <- 0.10
filter <- as.logical((ed_hs+(educ==0)) * married * (1-black))
system.time(rets.birthweight <- quantile.inf.np(Y=weight.g[filter],X=cbind(mom_age[filter],smoke[filter]),x0s=seq(from=22,to=37,by=1),p=p,ALPHA=ALPHA,JOINT.FLAG=TRUE)) #~15 seconds
rets.birthweight[[1]]$X.discrete <- ifelse(rets.birthweight[[1]]$X.discrete==1,"(smoker)","(non-smoker)")
rets.birthweight[[2]]$X.discrete <- ifelse(rets.birthweight[[2]]$X.discrete==1,"(smoker)","(non-smoker)")
plot2.quantile.inf.np(rets.birthweight,plot.data=FALSE,CI.interpolate=TRUE,CI.type='joint',x.axis.title="Mother's age (years)",y.axis.title="Birthweight (g)",title=sprintf("L-stat %g%% CI for %g-quantile\n Mother: married, white, HS grad or less",100*(1-ALPHA),p))

rets.birthweight[[1]]$bandwidth.joint
rets.birthweight[[2]]$bandwidth.joint
#When the bandwidths are greater than 1 -> not smoothing over age at all, just using observations with that exact age.  How does this compare to Li and Racine's recommendations for when to smooth over ordinal discrete variables?

detach(nat)


#
# Example E.3: US wage inequality
#
library(foreign)
#Data: http://economics.mit.edu/faculty/angrist/data1/data/angchefer06
cps80 <- read.dta("census80.dta",convert.factors=FALSE)
cps90 <- read.dta("census90.dta",convert.factors=FALSE)
cps00 <- read.dta("census00.dta",convert.factors=FALSE)
source("quantile_inf_np.R")
ALPHA <- 0.1
for (i in 1:3) {
if (i==1) {cps<-cps80} else if (i==2) {cps<-cps90} else {cps<-cps00}
filter <- as.logical((cps$educ==12)) #cps00$black * 
#
p <- 0.75
system.time(rets.wages.hi <- quantile.inf.np(Y=cps$logwk[filter],X=cbind(cps$exper[filter],cps$black[filter]),x0s=seq(from=14,to=35,by=1),p=p,ALPHA=ALPHA))
rets.wages.hi[[1]]$X.discrete <- ifelse(rets.wages.hi[[1]]$X.discrete==1,"black","white")
rets.wages.hi[[2]]$X.discrete <- ifelse(rets.wages.hi[[2]]$X.discrete==1,"black","white")
# plot2.quantile.inf.np(rets.wages.hi,plot.data=T,CI.interpolate=T,x.axis.title="Experience (years)",y.axis.title="Log wage")
#
rets.wages.all <- vector("list",4L)
rets.wages.all[[1]] <- rets.wages.hi[[1]]
rets.wages.all[[2]] <- rets.wages.hi[[2]]
rets.wages.all[[1]]$X.discrete <- sprintf("%s/%g",rets.wages.all[[1]]$X.discrete,p)
rets.wages.all[[2]]$X.discrete <- sprintf("%s/%g",rets.wages.all[[2]]$X.discrete,p)
#
p <- 0.25
system.time(rets.wages.lo <- quantile.inf.np(Y=cps$logwk[filter],X=cbind(cps$exper[filter],cps$black[filter]),x0s=seq(from=14,to=35,by=1),p=p,ALPHA=ALPHA)) #~10 seconds
rets.wages.lo[[1]]$X.discrete <- ifelse(rets.wages.lo[[1]]$X.discrete==1,"black","white")
rets.wages.lo[[2]]$X.discrete <- ifelse(rets.wages.lo[[2]]$X.discrete==1,"black","white")
# plot2.quantile.inf.np(rets.wages.lo,plot.data=T,CI.type='joint',CI.interpolate=T,x.axis.title="Experience (years)",y.axis.title="Log wage")
#
rets.wages.all[[3]] <- rets.wages.lo[[1]]
rets.wages.all[[4]] <- rets.wages.lo[[2]]
rets.wages.all[[3]]$X.discrete <- sprintf("%s/%g",rets.wages.all[[3]]$X.discrete,p)
rets.wages.all[[4]]$X.discrete <- sprintf("%s/%g",rets.wages.all[[4]]$X.discrete,p)
#
yr <- ifelse(i==1,"1980",ifelse(i==2,"1990","2000"))
plot2.quantile.inf.np(rets.wages.all,plot.data=FALSE,CI.type='joint',CI.interpolate=TRUE,x.axis.title="Experience (years)",y.axis.title="Log wage",title=sprintf("L-stat %g%% CI\n Male HS Grads Ages 40-49 in %s Census",100*(1-ALPHA),yr),ylim=c(4,7),legpos='bottom')

#Any smoothing over experience here (i.e. bandwidth>1)?  (See comment at end of birthweight analysis re: smoothing over age.)  Looks like yes this time.
for (k in 1:4) {cat(rets.wages.all[[k]]$bandwidth.joint);cat('\n')}
}





#
# Example E.4a: health outcome (hemoglobin) vs. HH expenditure from IFLS, w/ or w/o HH head edu
#
#For exp_hb_edu_merge_files.do file and resulting exp_hb_edu.dta file: https://drive.google.com/drive/folders/0B-_LUSJVBv20bTNmX3NYeFJKZ2c?usp=sharing
#Data from: Indonesia Family Life Survey (IFLS) 2007, prepped with exp_hb_edu_merge_files.do
#http://www.rand.org/labor/FLS/IFLS.html
#Similar to Li, Lin, and Racine (2013), but not restricted to Sundanese, and possibly other differences with filtering/accounting; results look quite different than w/ their dataset.
#Also need to restrict to one person per hh (say, youngest) to make iid assumption more plausible (else, probably strong intra-household correlation, e.g. from similar genetics).  Even then, could argue for spatial correlation.
library(foreign) #for .dta input
source("quantile_inf_np.R")
ifls.raw <- read.dta("exp_hb_edu.dta")
#limit to children
ifls.kids <- ifls.raw[ifls.raw$age<=15 & ifls.raw$age>=2,]
# #limit to youngest member of each hh
# ifls <- ifls.raw[order(ifls.raw$hhid07,ifls.raw$age),]
#limit to oldest child of each hh
ifls <- ifls.kids[order(ifls.kids$hhid07,-ifls.kids$age),]
ifls <- ifls[!duplicated(ifls$hhid07),]
#log expenditure; filter
logexp <- log(ifls$exp_hh_pc_yr)
filter <- (!is.na(ifls$hb_adj) & !is.na(ifls$exp_hh_pc_yr)) &
         # (1-ifls$male) #&
          ifls$age<=15
#By education (dl06: highest level attended)
#02=elementary,03=jr high general,4=jr hi vocational,5=sr hi/gen,6=sr hi/voc;
#60=college,61=univ(bachelor),62=univ(master),63=univ(phd); 11=adult edu A,12=adult edu B,15=adult edu C; 13=open univ; 14=islamic (pesantren),17=disabled,72=islamic elementary,73=islamic jr hi,74=islamic sr hi,90=kindergarten,98=don't know, other=95
edu.hi <- as.numeric(ifls$dl06_hh>=60 & ifls$dl06_hh<=63) #+
          # as.numeric(ifls$dl06_hh>=5 & ifls$dl06_hh<=6) # +...==74
filter0 <- !is.na(ifls$hb_adj) & !is.na(ifls$exp_hh_pc_yr)
filter1 <- !(ifls$dl06_hh>=11 & ifls$dl06_hh<=13) & (ifls$dl06_hh!=15) & (ifls$dl06_hh!=17) & (ifls$dl06_hh!=95)
filter1[is.na(filter1)] <- FALSE
filter2 <- (ifls$dl06_hh<=4 | edu.hi==1)
filter2[is.na(filter2)] <- FALSE
# filter2 <- (ifls$dl06_hh==2 | edu.hi==1) # | ifls$dl06_hh==72
# filter2[is.na(filter2)] <- FALSE
filter <- filter0 & filter2 & (ifls$age<=15)

#
ALPHA<-0.1
# for (p in c(0.15)) { #0.23, 0.03
p <- 0.15

x0s.p01 <- seq(from=13.5,to=17.25,by=0.25) #near severe cutoff
x0s.p05 <- seq(from=13.25,to=17.5,by=0.25)   #near moderate cutoff
x0s.p10 <- seq(from=12.75,to=20.25,by=0.25)      #near XXXXXXX
x0s.p10 <- seq(from=13.25,to=17.25,by=0.25)      #near XXXXXXX
x0s.p15 <- seq(from=12.75,to=17.25,by=0.25)      #near XXXXXXX
x0s.p25 <- seq(from=12.5,to=21.25,by=0.25)      #near XXXXXXX
x0s.p25 <- seq(from=13.25,to=17.25,by=0.25)      #near XXXXXXX
if (p<=0.01) {x0s<-x0s.p01} else if (p<=0.05) {x0s<-x0s.p05} else if (p<=0.1) {x0s<-x0s.p10} else if (p<=0.15) {x0s<-x0s.p15} else {x0s<-x0s.p25}

# p<-0.15;  ALPHA<-0.1

tmp <- quantile.inf(X=list(t=ifls$hb_adj[filter & edu.hi],c=ifls$hb_adj[filter & !edu.hi]), p=p, ONESIDED=0, ALPHA=ALPHA, METHOD.TYPE='qte', NORM.APPROX=T) #identical CI if NORM.APPROX=F
tmp 
cbind(quantile(ifls$hb_adj[filter & edu.hi],p),quantile(ifls$hb_adj[filter & !edu.hi],p))

x0s.p01 <- seq(from=13.5,to=17.25,by=0.25) #near severe cutoff
x0s.p05 <- seq(from=13.25,to=19.75,by=0.25)   #near moderate cutoff
x0s.p10 <- seq(from=13.25,to=16.75,by=0.5)      #near XXXXXXX
x0s.p25 <- seq(from=13.5,to=17.5,by=0.25)      #near XXXXXXX
if (p==0.01) {x0s<-x0s.p01} else if (p==0.05) {x0s<-x0s.p05} else if (p==0.1) {x0s<-x0s.p10} else {x0s<-x0s.p25}
#
# Threshold plotting fn; thresholds from http://www.who.int/vmnis/indicators/haemoglobin.pdf
plot.thresh.fn <- function(x) {
  lines(x=c(-100,100),y=c(13,13),lwd=4,col="navy") #non-anemia
  lines(x=c(-100,100),y=c(11,11),lwd=3,col="navy") #cutoff between "mild" and "moderate" anemia
  lines(x=c(-100,100),y=c(8,8),lwd=3,col="navy") #cutoff for "severe"
  text(x,20,"Non-anemia",pos=3,cex=1.4)
  text(x,11.1,"Mild",pos=3,cex=1.4)
  text(x,9,"Moderate",pos=3,cex=1.4)
  text(x,6,"Severe",pos=3,cex=1.4)
}
#
# Separate single-quantile CIs
system.time(rets.ifls.hb.edu <- quantile.inf.np(Y=ifls$hb_adj[filter],X=cbind(logexp[filter],edu.hi[filter]),x0s=x0s,p=p,ALPHA=ALPHA)) #~6 seconds
rets.ifls.hb.edu[[1]]$X.discrete <- ifelse(rets.ifls.hb.edu[[1]]$X.discrete==1,"(high edu)","(low edu)")
rets.ifls.hb.edu[[2]]$X.discrete <- ifelse(rets.ifls.hb.edu[[2]]$X.discrete==1,"(high edu)","(low edu)")
plot2.quantile.inf.np(rets.ifls.hb.edu,plot.data=F,CI.type='pointwise',CI.interpolate=TRUE,x.axis.title="Ln per capita expenditure (Rp/yr)",y.axis.title="Adjusted hemoglobin (g/dl)",title=sprintf("L-stat %g%% CI for %g-quantile of children",100*(1-ALPHA),p),xlim=c(12,22),ylim=c(0,22), datacex=0.15)
plot.thresh.fn(21)
#
# CQTE inference
source('quantile_inf_np.R')
system.time(rets2.ifls.hb.edu <- quantile.inf.np(Y=ifls$hb_adj[filter],X=cbind(logexp[filter],edu.hi[filter]),x0s=x0s,p=p,ALPHA=ALPHA,METHOD.TYPE='qte')) #~17 seconds
plot.quantile.inf.np(rets2.ifls.hb.edu,plot.data=T,CI.type='both',CI.interpolate=FALSE,x.axis.title="Ln per capita expenditure (Rp/yr)",y.axis.title="Adjusted hemoglobin difference (g/dl)",title=sprintf("L-stat %g%% CI for %g-quantile treatment\neffect of HH head's education, in children",100*(1-ALPHA),p),xlim=c(12,22),ylim=c(-1,22), datacex=0.25)
plot.thresh.fn(21)
lines(x=c(-100,18),y=c(0,0),lwd=1,col=1)
Y <- ifls$hb_adj[filter & logexp<17];  X <- logexp[filter & logexp<17];  Tr <- edu.hi[filter & logexp<17]
Y.t <- Y[!!Tr];  X.t <- X[!!Tr];  Y.c <- Y[!Tr];  X.c <- X[!Tr]
lines(sort(X.t),predict.rqss(rqss(Y.t~qss(X.t),tau=p),data.frame(X.t=X.t))[order(X.t)],col=6,lwd=2)
lines(sort(X.c),predict.rqss(rqss(Y.c~qss(X.c),tau=p),data.frame(X.c=X.c))[order(X.c)],col=5,lwd=2)


#pooled
Y <- ifls$hb_adj[filter]
X <- logexp[filter]
system.time(rets.ifls.hb <- quantile.inf.np(Y=Y,X=X,x0s=x0s,p=p,ALPHA=ALPHA))
# system.time(rets.ifls.hb <- quantile.inf.np(Y=ifls$hb_adj[filter],X=logexp[filter],x0s=x0s,p=p,ALPHA=ALPHA))
plot.quantile.inf.np(rets.ifls.hb,plot.data=TRUE,CI.type='pointwise',CI.interpolate=TRUE,x.axis.title="Ln per capita expenditure (Rp/yr)",y.axis.title="Adjusted hemoglobin (g/dl)",title=sprintf("L-stat %g%% CI for %g-quantile",100*(1-ALPHA),p),xlim=c(12,22),ylim=c(0,22), datacex=0.1)
#
# ### BEGIN RQSS (5min) ###
  # n <- length(Y)
  # g <- function(lam,y,x) AIC(rqss(y ~ qss(x, lambda = lam)),k = -1) #k=-1 does BIC (look @ code for AIC.rqss)
  # cat(sprintf("n=%d: ",n))
  # cat(paste0("BIC=",sprintf("%8.1f",system.time(
    # lamstar <- optimize(g, interval = c(0.001, 0.5), x = X, y = Y) )[1]),";"))
  # cat(paste0("fit=",sprintf("%8.1f",system.time(
    # rqss.fit <- rqss(Y~qss(X, lambda=lamstar$min), tau=p) )[1]) ,";"))
  # flush.console()
  # cat("\n")
  # tmp <- predict.rqss(rqss.fit,newdata=data.frame(X=x0s),interval="confidence",level=1-ALPHA)
  # rqss.ci.pt <- tmp[,2:3]
  # # plot(rqss.fit,bands="both")
# lines(type='o',x=x0s,y=rqss.ci.pt[,1],lwd=1,col=4,lty=4,pch=4)
# lines(type='o',x=x0s,y=rqss.ci.pt[,2],lwd=1,col=4,lty=4,pch=4)
# ### END RQSS ###

#thresholds from http://www.who.int/vmnis/indicators/haemoglobin.pdf
lines(x=c(10,30),y=c(13,13),lwd=4,col="navy") #non-anemia
lines(x=c(10,30),y=c(11,11),lwd=3,col="navy") #cutoff between "mild" and "moderate" anemia
lines(x=c(10,30),y=c(8,8),lwd=3,col="navy") #cutoff for "severe"
text(21,20,"Non-anemia",pos=3,cex=1.4)
text(21,11.1,"Mild",pos=3,cex=1.4)
text(21,9,"Moderate",pos=3,cex=1.4)
text(21,6,"Severe",pos=3,cex=1.4)
# }
# rets.ifls.hb[[1]]$pointwise.method.Lstat
# rets.ifls.hb[[1]]$joint.method.Lstat




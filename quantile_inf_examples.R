# Load code
prob.loaded <- exists("quantile.inf")
success <-
  tryCatch({source('https://raw.githubusercontent.com/kaplandm/R/main/quantile_inf.R');TRUE},
           error=function(w) FALSE)
if (!success) {
  success <-
    tryCatch({source("quantile_inf.R");TRUE}, error=function(w) FALSE)
  if (!success) {
    if (prob.loaded) {
      warning("Couldn't load quantile_inf.R, but it seems like you already did.")
    } else {
      stop("Failed to source() quantile_inf.R from web or local file.  You may download and source() it yourself, or at least make sure it's in your getwd().  Currently available at https://github.com/kaplandm/R")
    }
  }
}

# Examples of using quantile.inf(), with both simulated and observed data

##################
# Simulated data #
##################

# EXAMPLE S.1: single and joint
set.seed(112358)
X <- rnorm(123)
CI.p50 <- quantile.inf(X=X,p=0.50,NORM.APPROX=FALSE,SPACING.FLAG=TRUE)
CI.p95 <- quantile.inf(X=X,p=0.95,SPACING.FLAG=TRUE)
system.time(CI.p50p95 <- quantile.inf(X=X,p=c(0.50,0.95),METHOD.TYPE='joint',SPACING.FLAG=TRUE,BETA.BLOCK.SIZE=10^4))
cbind(rbind(CI.p50$CI,CI.p95$CI),CI.p50p95$CI)
#
nrep <- 400
CI.p50 <- CI.p95 <- mat.or.vec(nrep,2)
CI.joint <- mat.or.vec(nrep,4)
set.seed(112358)
for (irep in 1:nrep) {
  X <- rnorm(123)
  CI.p50[irep,] <- unlist(quantile.inf(X=X,p=0.50,NORM.APPROX=FALSE,SPACING.FLAG=TRUE)$CI)
  CI.p95[irep,] <- unlist(quantile.inf(X=X,p=0.95,SPACING.FLAG=TRUE)$CI)
  tmp <- quantile.inf(X=X,p=c(0.50,0.95),METHOD.TYPE='joint',SPACING.FLAG=TRUE,BETA.BLOCK.SIZE=10^3)$CI
  CI.joint[irep,] <- unlist(c(tmp[1,],tmp[2,]))
}
mean(CI.p50[,1]<qnorm(0.5) & qnorm(0.5)<CI.p50[,2]) #nominal: 95%
mean(CI.p95[,1]<qnorm(0.95) & qnorm(0.95)<CI.p95[,2]) #nominal: 95%
mean(CI.p50[,1]<qnorm(0.5) & qnorm(0.5)<CI.p50[,2] & CI.p95[,1]<qnorm(0.95) & qnorm(0.95)<CI.p95[,2]) #should be below 95%, which is why we need the 'joint' METHOD.TYPE -->
mean(CI.joint[,1]<qnorm(0.5) & qnorm(0.5)<CI.joint[,2] & CI.joint[,3]<qnorm(0.95) & qnorm(0.95)<CI.joint[,4]) #nominal: 95%


# EXAMPLE S.2: QTE with different options
set.seed(112358)
X <- list(t=rt(31,3)+1, c=rt(51,3))
system.time(CI.p50a <- quantile.inf(X=X,p=0.50,METHOD.TYPE='qte', NORM.APPROX=TRUE,SPACING.FLAG=TRUE))
system.time(CI.p50b <- quantile.inf(X=X,p=0.50,METHOD.TYPE='qte', NORM.APPROX=TRUE,SPACING.FLAG=FALSE))
system.time(CI.p50c <- quantile.inf(X=X,p=0.50,METHOD.TYPE='qte', NORM.APPROX=FALSE,SPACING.FLAG=TRUE, BETA.BLOCK.SIZE=10^4))
system.time(CI.p50d <- quantile.inf(X=X,p=0.50,METHOD.TYPE='qte', NORM.APPROX=FALSE,SPACING.FLAG=FALSE, BETA.BLOCK.SIZE=10^4))
rbind(CI.p50a$CI,CI.p50b$CI,CI.p50c$CI,CI.p50d$CI)


# EXAMPLE S.3: QTE not at median
set.seed(112358)
X <- list(t=rt(51,1)+1, c=rt(71,3))
ret <- quantile.inf(X=X,p=0.75,METHOD.TYPE='qte', NORM.APPROX=TRUE,SPACING.FLAG=FALSE)
ret
#
nrep <- 1600
QTE <- 1
CI.p50 <- CI.p75 <- mat.or.vec(nrep,2)
set.seed(112358)
for (irep in 1:nrep) {
  X <- list(t=rt(31,3)+QTE, c=rt(51,3))
  CI.p50[irep,] <- unlist(quantile.inf(X=X,p=0.50,METHOD.TYPE='qte', NORM.APPROX=TRUE,SPACING.FLAG=TRUE)$CI)
  CI.p75[irep,] <- unlist(quantile.inf(X=X,p=0.75,METHOD.TYPE='qte', NORM.APPROX=TRUE,SPACING.FLAG=TRUE)$CI)
}
mean(CI.p50[,1]<QTE & QTE<CI.p50[,2]) #nominal: 95%
mean(CI.p75[,1]<QTE & QTE<CI.p75[,2]) #nominal: 95%


# EXAMPLE S.4: lincom (ex: inference on interquartile range)
set.seed(112358)
X <- runif(101,0,1)
system.time(ret <- quantile.inf(X=X,p=c(0.75,0.25),PSI=c(1,-1),METHOD.TYPE='lincom', SPACING.FLAG=TRUE))
ret
#
# REPEATED SIMULATION
nrep <- 200
CI.IQR <- mat.or.vec(nrep,2)
pIQR <- c(0.75,0.25)
set.seed(112358)
for (irep in 1:nrep) {
  X <- runif(101)
  CI.IQR[irep,] <- unlist(quantile.inf(X=X,p=pIQR,METHOD.TYPE='lincom', PSI=c(1,-1), BETA.BLOCK.SIZE=10^3)$CI)
}
mean(CI.IQR[,1]<(qunif(pIQR[1])-qunif(pIQR[2])) & (qunif(pIQR[1])-qunif(pIQR[2]))<CI.IQR[,2]) #nominal: 95%


# EXAMPLE S.5: treatment effect on IQR
set.seed(112358)
X <- list(t=rbeta(101,5,8), c=rbeta(101,5,8)+1)
system.time(ret <- quantile.inf(X=X,p=c(0.75,0.25),PSI=c(1,-1),METHOD.TYPE='qte'))
ret
#
# REPEATED SIMULATION
nrep <- 200
CI.tIQR <- mat.or.vec(nrep,2)
set.seed(112358)
system.time(
for (irep in 1:nrep) {
  X <- list(t=rbeta(101,5,8), c=rbeta(101,5,8)+1)
  CI.tIQR[irep,] <- unlist(quantile.inf(X=X,p=c(0.75,0.25),PSI=c(1,-1),METHOD.TYPE='qte',BETA.BLOCK.SIZE=10^3)$CI)
})
mean(CI.tIQR[,1]<0 & 0<CI.tIQR[,2]) #nominal: 95%






##################
# Empirical data #
##################

#
# Example E.1: Gneezy and List (2006)
#
## Import data: typed from Gneezy and List (2006, page 1371)
# Two matrices: noGiftLib (10 participants) and GiftLib(9 participants)
# For each participant, have productivity (#books logged) per 90-minute
# increment, measured at 90, 180, 270, and 360 minutes.  For fundraising task, units are dollars raised per period, first two periods shown.
noGiftLib <- c(56,61,58,63,  52,52,51,45,  46,44,52,42,  45,41,43,38, 
          41,29,33,25,  38,42,44,46,  37,39,38,38,  34,35,32,37, 
          32,32,28,27,  26,30,33,35)
dim(noGiftLib) <- c(4,10)
noGiftLib <- t(noGiftLib)
#
GiftLib  <- c(75,71,60,58,  64,65,63,61,  63,65,59,63,  58,40,35,31, 
          54,42,33,34,  47,35,28,25,  42,37,47,39,  37,29,30,30, 
          25,20,20,22)
dim(GiftLib) <- c(4,9)
GiftLib <- t(GiftLib)
#
noGiftFund <- c(6,7, 6,21, 20,24, 35,15, 6,25, 8,13, 0,4, 41,25, 49,51, 21,14)
dim(noGiftFund) <- c(2,10)
noGiftFund <- t(noGiftFund)
#
GiftFund <- c(35,26, 32,34, 31,20, 14,19, 27,17, 42,25, 31,11, 26,3,  15,25, 42,16, 77,19, 29,33, 28,26)
dim(GiftFund) <- c(2,13)
GiftFund <- t(GiftFund)

# Period 1 quartile treatment effects for Lib and Fund, 90% nominal 2-sided CI
APPROX <- TRUE
ps <- c(0.25,0.5,0.75)
disp <- cbind(ps,lo=NA,hi=NA)
for (i in 1:length(ps)) {
  disp[i,2:3] <- unlist(quantile.inf(X=list(t=GiftLib[,1],c=noGiftLib[,1]),p=ps[i],METHOD.TYPE='qte', NORM.APPROX=APPROX, SPACING.FLAG=TRUE, ONESIDED=0, ALPHA=0.10, BETA.BLOCK.SIZE=10^5)$CI)
}
disp
# Library, Period 2
APPROX <- TRUE
ps <- c(0.25,0.5,0.75)
disp <- cbind(ps,lo=NA,hi=NA)
for (i in 1:length(ps)) {
  disp[i,2:3] <- unlist(quantile.inf(X=list(t=GiftLib[,2],c=noGiftLib[,2]),p=ps[i],METHOD.TYPE='qte', NORM.APPROX=APPROX, SPACING.FLAG=TRUE, ONESIDED=0, ALPHA=0.10, BETA.BLOCK.SIZE=10^5)$CI)
}
disp
#
# Fundraising, Period 1
APPROX <- TRUE
ps <- c(0.25,0.5,0.75)
disp <- cbind(ps,lo=NA,hi=NA)
for (i in 1:length(ps)) {
  disp[i,2:3] <- unlist(quantile.inf(X=list(t=GiftFund[,1],c=noGiftFund[,1]),p=ps[i],METHOD.TYPE='qte', NORM.APPROX=APPROX, SPACING.FLAG=TRUE, ONESIDED=0, ALPHA=0.10, BETA.BLOCK.SIZE=10^5)$CI)
}
disp
# Fundraising, Period 2
APPROX <- TRUE
ps <- c(0.25,0.5,0.75)
disp <- cbind(ps,lo=NA,hi=NA)
for (i in 1:length(ps)) {
  disp[i,2:3] <- unlist(quantile.inf(X=list(t=GiftFund[,2],c=noGiftFund[,2]),p=ps[i],METHOD.TYPE='qte', NORM.APPROX=APPROX, SPACING.FLAG=TRUE, ONESIDED=0, ALPHA=0.10, BETA.BLOCK.SIZE=10^5)$CI)
}
disp


#
# Example E.2: birthweight: smokers vs. non-smokers (filtered on <=HS edu, married, white)
#
#Load data per Chernozhukov & Fernandez-Val (2011) code
#Data: http://restud.oxfordjournals.org/content/78/2/559/suppl/DC1
library(foreign) #for .dta input
nat <- read.dta("Birthweights_data.dta")
nat$educ <- nat$ed_hs + 2*nat$ed_smcol + 3*nat$ed_col
nat <- nat[order(-nat$smoke),]
attach(nat) #make variables accessible by just using names
#weight(=depvar,in kg to start),id,birmon(=6,June only),married(0/1),{tri1,tri2,tri3,novisit},m_wtgain(lbs,integer),{ed_hs,ed_smcol,ed_col},mom_age,black(0/1),boy(0/1),smoke(0/1),cigsper(integer but grouped at 5's)
weight.g          <- 1000*weight; #was in kg(?), now g
m.wt.gain.lbs.net <- m_wtgain - 2.20462*weight
#
ALPHA <- 0.05;  p<-0.10
filter <- as.logical((ed_hs+(educ==0)) * married * (1-black))
disc.vals <- data.frame(smoke=c(0,1,0,1),ed_hs=c(0,0,1,1))
disp <- cbind(disc.vals,lo=NA,hi=NA)
for (i in 1:4) {
  disp[i,3:4] <- unlist(quantile.inf(X=weight.g[filter & smoke==disc.vals$smoke[i] & ed_hs==disc.vals$ed_hs[i]],p=p,ALPHA=ALPHA)$CI)
}
disp
#
#difference b/w birthweights (in grams) for smokers & nonsmokers given HS grad (only), or given not completed HS
quantile.inf(X=list(t=weight.g[filter & ed_hs==1 & smoke==1],c=weight.g[filter & ed_hs==1 & smoke==0]),p=p,ALPHA=ALPHA,METHOD.TYPE='qte')$CI
quantile.inf(X=list(t=weight.g[filter & ed_hs==0 & smoke==1],c=weight.g[filter & ed_hs==0 & smoke==0]),p=p,ALPHA=ALPHA,METHOD.TYPE='qte')$CI
#
#difference b/w HS grads & non-grads, given non-smoker or given smoker
quantile.inf(X=list(t=weight.g[filter & ed_hs==1 & smoke==0],c=weight.g[filter & ed_hs==0 & smoke==0]),p=p,ALPHA=ALPHA,METHOD.TYPE='qte')$CI
quantile.inf(X=list(t=weight.g[filter & ed_hs==1 & smoke==1],c=weight.g[filter & ed_hs==0 & smoke==1]),p=p,ALPHA=ALPHA,METHOD.TYPE='qte')$CI



# EOF
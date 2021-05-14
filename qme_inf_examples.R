# Empirical example for qme_inf.R, "Nonparametric inference on quantile marginal effects" (David M. Kaplan)
# Quantile Engel elasticities from UK LCFS 2012, http://esds.ac.uk/doi/?sn=7472#2

SAVE.FLAG <- TRUE

library(foreign)
source("qme_inf.R") # from https://github.com/kaplandm/R
raw <- vector('list',12)
for (i in 1:12) {
  raw[[i]] <- read.dta(sprintf("20%02d_dvhh_ukanon.dta",i),convert.factors=FALSE)
}
moneycols <- c('P630tp','P539t','P601t','P607t','P611t',
               'P604t','P609t','P603t','P605t') 
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
filter <- filter & f.under60 & (f.1a | f.2a) & (f.0c)
#filter <- filter & f.under65 & f.renter & f.1a & f.0c
#filter <- filter & f.under60 & (f.1a) & (f.0c)
#filter <- filter & f.under60 & (f.2a) & (f.0c | f.0c)
sum(filter)

expf <- dat$p630tp[filter]
lnexpf <- log(expf)
x0s <- quantile(lnexpf,p=1:9/10)
ps <- c(0.5,0.75,0.9) #1:3/4

cats <- rbind(c('p601t','p539t','p607t','p611t',
                'p604t','p609t','p603t','p605t'),
              c('Food and non-alcoholic beverage','Alcohol','Transportation','Restaurants+Hotels',
                'Housing+utilities','Recreation','Clothing/footwear','Household goods and services'),
              c('food','alc','transport','restaurant','housing','rec','clothing','HH'))
for (catind in 1:4) {
  set.seed(112358)
  w <- dat[[cats[1,catind]]][filter]/expf
  w <- ifelse(w<0,0,w)
  lnw <- log(w) # log(0)=-Inf
  if (FALSE) {
    par(family='serif',mar=c(5.0,6.0,6.0,2.1))
    plot(lnexpf,ifelse(lnw==-Inf,max(-7,min(lnw[lnw!=-Inf])),lnw),main=cats[2,catind],
         pch=16,cex=0.1,ylim=c(-7,0),
         cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
         xlab="Log total expenditure",ylab="Log budget share")
    subind <- sample(1:length(lnexpf),2000,T)
    if (SAVE.FLAG) pdf(file=sprintf("qme_inf_ex_%s_subscatter.pdf",cats[3,catind])) #else draw to current device
    par(family='serif',mar=c(5.0,6.0,6.0,2.1))
    plot(lnexpf[subind],ifelse(lnw[subind]==-Inf,-7+runif(length(lnw),-0.2,0.2),lnw[subind]),main=sprintf("%s:\nsubsample of %d",cats[2,catind],length(subind)),
         pch=16,cex=0.1,ylim=c(-7,0),xlim=range(lnexpf),
         cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
         xlab="Log total expenditure",ylab="Log budget share")
    if (SAVE.FLAG) dev.off()
  }
  CI.los <- CI.his <- matrix(NA,length(x0s),length(ps))
  for (pind in 1:length(ps)) {
    ret <- qme.inf(Y=lnw,X=lnexpf,x0s=x0s,p=ps[pind],ALPHA=0.1,ONESIDED=0) 
                   #,hs=rep.int(0.1,length(x0s)))
    CI.los[,pind] <- 1+ret$CI.lo #transform to get dln(q)/dln(y)
    CI.his[,pind] <- 1+ret$CI.hi
  }
  if (SAVE.FLAG) pdf(file=sprintf("qme_inf_ex_%s.pdf",cats[3,catind])) #else draw to current device
  par(family='serif',mar=c(5.0,6.0,6.0,2.1))
  tmp <- sort(c(CI.los,CI.his))
  #plot(rep.int(x0s,2*length(ps)),c(CI.los,CI.his),type='n',
  plot(x0s[c(1,length(x0s))],quantile(tmp,p=c(0.3,0.7)),type='n',
       cex.main=2.0, mgp=c(3.5,1,0), cex.axis=1.5, cex.lab=2,
       main=cats[2,catind], 
       xlab="log expenditure", ylab="elasticity",ylim=c(-2,4))
  for (pind in 1:length(ps)) {
    lines(x0s,CI.los[,pind],lwd=2,lty=pind,pch=pind,col=1,type='o',cex=1.5)
    lines(x0s,CI.his[,pind],lwd=2,lty=pind,pch=pind,col=1,type='o',cex=1.5)
  }
  legend("bottom",sprintf("%4.2f-quantile",ps),inset=c(0.0,0.0), bty='n',
         col=1, pch=1:length(ps), lty=1:length(ps), lwd=2, cex=1.3)
  if (SAVE.FLAG) dev.off()
  #
  lengths <- (CI.his-CI.los)
  cat(sprintf("CI lengths by x0s for %s",cats[2,catind]),sep='\n',append=TRUE)
  for (pind in 1:length(ps)) {
    cat(sprintf("p=%4.2f: ",ps[pind]),append=TRUE,sep='')
    cat(sprintf("%4.1f,",lengths[,pind]),append=TRUE)
    cat('\n')
  }
}
cat(sum(filter))

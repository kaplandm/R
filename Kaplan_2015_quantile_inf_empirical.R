# Empirical example(s) using quantile.inf.K11() for Kaplan (2015)
# Files from https://github.com/kaplandm/R
source("quantile_inf_K11.R")
source("quantile_inf.R")

#
# Example E.1: VaR
#
#Google finance: http://www.google.com/finance/historical?cid=26001&startdate=Aug+29%2C+1986&enddate=Aug+31%2C+2012&num=200
#.csv available @ https://drive.google.com/file/d/0B-_LUSJVBv20anFLdTZpQlJfbms/view?usp=sharing
VaR.data.raw <- read.csv("OXY_daily_data_no_holidays.csv")
Y <- VaR.data.raw$DailyLnRet
ALPHA <- 0.1
p <- 0.05;  ns <- c(seq(from=24,to=74,by=5),seq(from=94,to=254,by=20),seq(from=252*2,to=252*10,by=252))
disp <- cbind(ns,lo=NA,hi=NA,H99used=NA,QIlo=NA,QIhi=NA)
#
for (i in 1:length(ns)) {
  disp[i,2:3] <- quantile.inf.K11(X=Y[1:ns[i]],p=p,ALPHA=ALPHA,ONESIDED=0,METHOD.TYPE='single')
  tmp <- quantile.inf(X=Y[1:ns[i]],p=p,ALPHA=ALPHA,ONESIDED=0,METHOD.TYPE='single')
  disp[i,4] <- (tmp$methname=='Hutson')
  disp[i,5:6] <- unlist(tmp$CI)
}
disp
melo <- (disp[,2]<disp[,5]) & (disp[,3]<disp[,6])
Hlo <- (disp[,2]>disp[,5]) & (disp[,3]>disp[,6])
melo
Hlo
sum(melo);  sum(Hlo)
for (i in 1:length(ns)) {
  cat(paste0(c(sprintf("%4d & & [%6.4f, %6.4f] & & [%6.4f, %6.4f] \\\\",
                       disp[i,1],disp[i,2],disp[i,3],disp[i,5],disp[i,6]))),sep="\n")
}

# EOF
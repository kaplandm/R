# Unconditional QTE empirical examples for Goldman and Kaplan (2017), "Nonparametric inference on (conditional) quantile differences and linear combinations, using L-statistics"
# Questions? Comments? kaplandm@missouri.edu

tryCatch(source("quantile_inf.R"),
         error=function(w) { cat("First download the file quantile_inf.R currently available from http://faculty.missouri.edu/~kaplandm/code/quantile_inf.R and save it in the same directory as quantile_inf_np.R (or elsewhere in R's path)."); stop(w) } )

# Input data from Gneezy and List (2006), Table I (p. 1371)
noGiftLib <- c(56,61,58,63,  52,52,51,45,  46,44,52,42,
               45,41,43,38,  41,29,33,25,  38,42,44,46, 
               37,39,38,38,  34,35,32,37,  32,32,28,27, 
               26,30,33,35)
dim(noGiftLib) <- c(4,10)
noGiftLib <- t(noGiftLib)
#
GiftLib  <- c(75,71,60,58,  64,65,63,61,  63,65,59,63, 
              58,40,35,31,  54,42,33,34,  47,35,28,25, 
              42,37,47,39,  37,29,30,30,  25,20,20,22)
dim(GiftLib) <- c(4,9)
GiftLib <- t(GiftLib)
# # output to .csv for ECTJ:
# Lib <- data.frame(GiftTreatmentDummy=c(rep(0,dim(noGiftLib)[1]),rep(1,dim(GiftLib)[1])),
#                   ParticipantIDNumber=1:19,
#                   Period1output=c(noGiftLib[,1],GiftLib[,1]),
#                   Period2output=c(noGiftLib[,2],GiftLib[,2]),
#                   Period3output=c(noGiftLib[,3],GiftLib[,3]),
#                   Period4output=c(noGiftLib[,4],GiftLib[,4]) )
# write.csv(Lib, row.names=FALSE, file="GL2006_TableI_library_data.csv")
#
#
# Table V, p. 1376 of Gneezy and List (2006)
noGiftFund <- c(6,7,  6,21, 20,24, 35,15, 6,25, 8,13, 
                0,4, 41,25, 49,51, 21,14)
dim(noGiftFund) <- c(2,10)
noGiftFund <- t(noGiftFund)
#
GiftFund <- c(35,26, 32,34, 31,20, 14,19, 27,17, 
              42,25, 31,11, 26,3,  15,25, 42,16, 
              77,19, 29,33, 28,26)
dim(GiftFund) <- c(2,13)
GiftFund <- t(GiftFund)
# output to .csv for ECTJ:
Fund <- data.frame(GiftTreatmentDummy=c(rep(0,dim(noGiftFund)[1]),rep(1,dim(GiftFund)[1])),
                  ParticipantIDNumber=1:23,
                  Period1output=c(noGiftFund[,1],GiftFund[,1]),
                  Period2output=c(noGiftFund[,2],GiftFund[,2]) )
write.csv(Fund, row.names=FALSE, file="GL2006_TableV_fundraising_data.csv")

# Period 1 quartile treatment effects for Lib and Fund, 90% nominal 2-sided CI
ps <- c(0.25,0.50,0.75)

# Library, Period 1
disp.spac <- disp.kde  <- cbind(ps,lo=NA,hi=NA)
for (i in 1:length(ps)) {
  disp.spac[i,2:3] <- unlist(quantile.inf(X=list(t=GiftLib[,1],c=noGiftLib[,1]),p=ps[i],METHOD.TYPE='qte', SPACING.FLAG=TRUE, ONESIDED=0, ALPHA=0.10)$CI)
  disp.kde[i,2:3]  <- unlist(quantile.inf(X=list(t=GiftLib[,1],c=noGiftLib[,1]),p=ps[i],METHOD.TYPE='qte', SPACING.FLAG=FALSE, ONESIDED=0, ALPHA=0.10)$CI)
}
# disp.spac
# disp.kde
tmp <- sprintf("(%6.2f&%6.2f) & (%6.2f&%6.2f) & (%6.2f&%6.2f) \\\\",
               disp.kde[1,2],disp.kde[1,3],
               disp.kde[2,2],disp.kde[2,3],
               disp.kde[3,2],disp.kde[3,3])
tmp <- gsub(" ","",tmp)
kde1 <- paste0("1 (kern/MSE)  &",tmp)
#
tmp <- sprintf("(%6.2f&%6.2f) & (%6.2f&%6.2f) & (%6.2f&%6.2f) \\\\",
               disp.spac[1,2],disp.spac[1,3],
               disp.spac[2,2],disp.spac[2,3],
               disp.spac[3,2],disp.spac[3,3])
tmp <- gsub(" ","",tmp)
spc1 <- paste0("2 (spac/bias) &",tmp)

# Library, Period 2
disp.spac <- disp.kde  <- cbind(ps,lo=NA,hi=NA)
for (i in 1:length(ps)) {
  disp.spac[i,2:3] <- unlist(quantile.inf(X=list(t=GiftLib[,2],c=noGiftLib[,2]),p=ps[i],METHOD.TYPE='qte', SPACING.FLAG=TRUE, ONESIDED=0, ALPHA=0.10)$CI)
  disp.kde[i,2:3]  <- unlist(quantile.inf(X=list(t=GiftLib[,2],c=noGiftLib[,2]),p=ps[i],METHOD.TYPE='qte', SPACING.FLAG=FALSE, ONESIDED=0, ALPHA=0.10)$CI)
}
# disp.spac
# disp.kde
tmp <- sprintf("(%6.2f&%6.2f) & (%6.2f&%6.2f) & (%6.2f&%6.2f) \\\\",
               disp.kde[1,2],disp.kde[1,3],
               disp.kde[2,2],disp.kde[2,3],
               disp.kde[3,2],disp.kde[3,3])
tmp <- gsub(" ","",tmp)
kde2 <- paste0("1 (kern/MSE)  &",tmp)
#
tmp <- sprintf("(%6.2f&%6.2f) & (%6.2f&%6.2f) & (%6.2f&%6.2f) \\\\",
               disp.spac[1,2],disp.spac[1,3],
               disp.spac[2,2],disp.spac[2,3],
               disp.spac[3,2],disp.spac[3,3])
tmp <- gsub(" ","",tmp)
spc2 <- paste0("2 (spac/bias) &",tmp)
cat(c("Library",kde1,spc1,kde2,spc2),sep='\n')



# Fundraising, Period 1
disp.spac <- disp.kde  <- cbind(ps,lo=NA,hi=NA)
for (i in 1:length(ps)) {
  disp.spac[i,2:3] <- unlist(quantile.inf(X=list(t=GiftFund[,1],c=noGiftFund[,1]),p=ps[i],METHOD.TYPE='qte', SPACING.FLAG=TRUE, ONESIDED=0, ALPHA=0.10)$CI)
  disp.kde[i,2:3]  <- unlist(quantile.inf(X=list(t=GiftFund[,1],c=noGiftFund[,1]),p=ps[i],METHOD.TYPE='qte', SPACING.FLAG=FALSE, ONESIDED=0, ALPHA=0.10)$CI)
}
# disp.spac
# disp.kde
tmp <- sprintf("(%6.2f&%6.2f) & (%6.2f&%6.2f) & (%6.2f&%6.2f) \\\\",
               disp.kde[1,2],disp.kde[1,3],
               disp.kde[2,2],disp.kde[2,3],
               disp.kde[3,2],disp.kde[3,3])
tmp <- gsub(" ","",tmp)
kde1 <- paste0("1 (kern/MSE)  &",tmp)
#
tmp <- sprintf("(%6.2f&%6.2f) & (%6.2f&%6.2f) & (%6.2f&%6.2f) \\\\",
               disp.spac[1,2],disp.spac[1,3],
               disp.spac[2,2],disp.spac[2,3],
               disp.spac[3,2],disp.spac[3,3])
tmp <- gsub(" ","",tmp)
spc1 <- paste0("2 (spac/bias) &",tmp)

# Fundraising, Period 2
disp.spac <- disp.kde  <- cbind(ps,lo=NA,hi=NA)
for (i in 1:length(ps)) {
  disp.spac[i,2:3] <- unlist(quantile.inf(X=list(t=GiftFund[,2],c=noGiftFund[,2]),p=ps[i],METHOD.TYPE='qte', SPACING.FLAG=TRUE, ONESIDED=0, ALPHA=0.10)$CI)
  disp.kde[i,2:3]  <- unlist(quantile.inf(X=list(t=GiftFund[,2],c=noGiftFund[,2]),p=ps[i],METHOD.TYPE='qte', SPACING.FLAG=FALSE, ONESIDED=0, ALPHA=0.10)$CI)
}
# disp.spac
# disp.kde
tmp <- sprintf("(%6.2f&%6.2f) & (%6.2f&%6.2f) & (%6.2f&%6.2f) \\\\",
               disp.kde[1,2],disp.kde[1,3],
               disp.kde[2,2],disp.kde[2,3],
               disp.kde[3,2],disp.kde[3,3])
tmp <- gsub(" ","",tmp)
kde2 <- paste0("1 (kern/MSE)  &",tmp)
#
tmp <- sprintf("(%6.2f&%6.2f) & (%6.2f&%6.2f) & (%6.2f&%6.2f) \\\\",
               disp.spac[1,2],disp.spac[1,3],
               disp.spac[2,2],disp.spac[2,3],
               disp.spac[3,2],disp.spac[3,3])
tmp <- gsub(" ","",tmp)
spc2 <- paste0("2 (spac/bias) &",tmp)
cat(c("Fundraising",kde1,spc1,kde2,spc2),sep='\n')

# EOF
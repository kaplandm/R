### SOME ANNOTATIONS BELOW; NO CODE CHANGES OTHER THAN VALUES OF I,x,tau,n, AND MAKING SIM NOT AUTO-RUN WHEN SOURCED
### DIFFERENCES WITH PAPER (SEE DETAILS BELOW)---1,3 MAKE NEGLIGIBLE DIFFERENCE:
# 1) RK APPEARS TO BE WRONG? (TOO BIG)
# 2) BASE BANDWIDTH IS NOT YU-JONES BUT KernSmooth::dpill BASED ON Ruppert, Sheather, and Wand (1995), WHICH IS FOR LOCAL LINEAR (MEAN) REGRESSION (AND DOES NOT DEPEND ON x), AND ASSUMES A GAUSSIAN KERNEL.
# 3) this is the (pi/2)^(1/5) from Yu Jones p. 231 ::  ADJUSTMENT IS an[2]=1.095*n^(-1/20) INSTEAD OF n^(-1/20)
### ALSO: NOTE WHEN CI ENDPOINT HAS QUANTILE<0 or >1, JUST USES 0 or 1, SO *WILL* COMPUTE FOR TAIL QUANTILES BUT PROBABLY GIVE BAD PERFORMANCE

#Confidence Interval in Fan&Liu
library('KernSmooth');

# some Auxiliary Functions

# Indicator Function
"Ind"<-function(Y,y){
  if(Y-y>0) a=0 else a=1;
  return(a)
}

# Choose Kernel Functioni
"BK"<-function(u){
  if(abs(u)<1) a=((1-u^2)^2)*(15/16)
  else a=0
  return(a)
}

#Conditional CDF Estimator in Fan & Liu
"FLc"<-function(X,Y,x,y,h){
  # Y is the depedent variable vector, x is the evaluation point
  # Indicator Function
  Indy<-function(Y){
    a=Ind(Y,y);
    return(a)
  }
  #Transformation of Variables
  Fx=ecdf(X);
  U=Fx(X);
  Uh=(U-Fx(x))/h;
  W=sapply(Uh,BK);
  Z=sapply(Y,Indy);
  r=mean(Z*W)/mean(W);
  return(r);
}

#Conditional Quantile Estimator in Fan & Liu
"FLqOLD"<-function(X,Y,x,tau,h){
  if (tau>1) f<-uniroot(function(y) FLc(X,Y,x,y,h)-1, f.lower=min(Y)-3, f.upper=max(Y),lower=min(Y)-1, upper=max(Y)+1, tol=0.00001)
  if (tau<0) f<-uniroot(function(y) FLc(X,Y,x,y,h), f.lower=min(Y)-3, f.upper=max(Y),lower=min(Y)-1, upper=max(Y)+1, tol=0.00001)
  if (tau<=1 & tau>=0) f<-uniroot(function(y) FLc(X,Y,x,y,h)-tau, f.lower=min(Y)-3, f.upper=max(Y),lower=min(Y)-1, upper=max(Y)+1, tol=0.00001)
  return(f$root)
}
FLq <- function(X,Y,x,tau,h,U,u0){
  Uh <- (U-u0)/h;
  W <- sapply(Uh,BK);
  W <- W/sum(W) #normalize sum to 1
  if (tau>1) return(max(Y[W>0]))
  if (tau<0) return(min(Y[W>0]))
  tmp <- order(Y);  Wsort <- W[tmp];  Ysort <- Y[tmp]
  cumW <- cumsum(Wsort)
  ind1 <- max(which(cumW<=tau))
  if (is.infinite(ind1)) return(Ysort[1]) else return(Ysort[ind1])
}

FL.CI.fn <- function(X,Y,x0,tau,ALPHA,UNIF=FALSE,ORIGINAL=TRUE) {
  RK <- 5/7;  
  n <- length(Y)
  z <- qnorm(1-ALPHA/2)
  FLdelta <- 1/5+1/20
  if (UNIF) {
    intKsq <- 5/7;  intKpsq <- 15/7  #mu2K <- 1/7;  
    z <- (log(2)-log(abs(log(1-ALPHA))))/sqrt(2*FLdelta*log(n)) + 
      sqrt(2*FLdelta*log(n)) + log(intKpsq/(4*pi*intKsq)) / sqrt(2*FLdelta*log(n)) #FL eqn (20)
  }
  Fx <- ecdf(X);  U <- Fx(X)
  h2 <- dpill(U,Y)*(pi/2)^(1/5)*n^(1/5-FLdelta)
  if (!ORIGINAL) {
    RKb <- 5/7; mu2Kb <- 1/7; RKg <- 1/(2*sqrt(pi)); mu2Kg <- 1
    h2 <- h2 * ((RKb/mu2Kb^2)/(RKg/mu2Kg^2))^(1/5) #dpill presumes Gaussian (not bisquare) kernel, so need to adjust
  }
  h2 <- h2*(2*tau*(1-tau)/(pi*dnorm(qnorm(tau))^2))^(1/5) #Yu and Jones (1998, bottom left of p. 231)
  Snp <- sqrt(tau*(1-tau)*RK/(n*h2))
  return(c(FLq(X,Y,x0,tau-z*Snp,h2,U,Fx(x0)), FLq(X,Y,x0,tau+z*Snp,h2,U,Fx(x0))))
}

if (!exists("RUNSIM")) RUNSIM <- FALSE

if (RUNSIM) {
  starttime <- Sys.time()
  #Simulation Results
  # Set Parameter Values
  n=500;                                #Sample Size
  I=100;                                 #Iteration Numbers
  tau=0.5;                              #Quantile Level
  x=0;                                 #Evaluation Point
  alpha=0.05;                            #Type I Error of CI
  
  an=c(1.13*n^(-1/20),1.095*n^(-1/20));       #Yu-Jones Adjustment in Bandwidth selection
  z=qnorm(1-alpha/2);
  RK=(15/8)*(44/45-4/7);               #Roughness of the Bisquare Kernel
  ### THIS SEEMS WRONG--SHOULD BE 5/7: http://www.wolframalpha.com/input/?i=integral+of+%28%2815%2F16%29*%281-u%5E2%29%5E2%29%5E2+from+u%3D-1+to+1
  RK.adj <- 5/7 ### ADDED
  
  for (tau in c(0.25,0.5,0.75)) {
    for (x in c(0,0.75,1.5)) {
      
      CIU=rep(0,I);                        #Null Vector to store upper/lower bounds of CI
      CIL=rep(0,I);
      
      CIU.adj <- CIL.adj <- rep.int(0,I)
      CIU.chk <- CIL.chk <- rep.int(0,I)
      
      set.seed(123);
      for (j in 1:I){
        X=rnorm(n);
        Z=rnorm(n);
        Y=2.5+sin(2*X)+2*exp(-16*(X)^2)+0.5*Z;
        
        Fx=ecdf(X);
        U=Fx(X);  
        
        
        h2=dpill(U,Y)*an[2];   ### dpill ASSUMES GAUSSIAN KERNEL, BUT USING BISQUARE KERNEL
        ### ALSO: PAPER SAYS n^(-1/20) ADJUSTMENT, BUT an[2]=1.095*n^(-1/20)
        h2.adj <- h2/an[2]*n^(-1/20) ### ADDED
        
        #Asymptotic Standard Deviations  
        
        Snp=sqrt(((tau*(1-tau))*RK)/(n*h2));  
        
        Snp.adj <- sqrt(tau*(1-tau)*RK.adj/(n*h2.adj)) ### ADDED
        
        #Calculate Lower/Upper Bounds of CI in Fan&Liu
#         CIL[j]=FLqOLD(X,Y,x,tau-z*Snp,h2);
#         CIU[j]=FLqOLD(X,Y,x,tau+z*Snp,h2);  
#         
#         CIL.adj[j] <- FLqOLD(X,Y,x,tau-z*Snp.adj,h2.adj)
#         CIU.adj[j] <- FLqOLD(X,Y,x,tau+z*Snp.adj,h2.adj)
        
        tmp <- FL.CI.fn(X,Y,x,tau,alpha)
        CIL.chk[j] <- tmp[1];  CIU.chk[j] <- tmp[2]
      }
      
      ### ADDED BY KAPLAN
      y0 <- 2.5+sin(2*x)+2*exp(-16*x^2) + 0.5*qnorm(tau)
      cat(sprintf("tau=%g,x=%g,n=%d,alpha=%g,I=%d\n",tau,x,n,alpha,I))
      cat(sprintf("CP=%6.4f,CP.adj=%6.4f,CP.chk=%6.4f; %g seconds elapsed\n",
                  1-mean(CIU<y0 | CIL>y0),1-mean(CIU.adj<y0 | CIL.adj>y0),1-mean(CIU.chk<y0 | CIL.chk>y0),
                  Sys.time()-starttime))
    }
  }
}

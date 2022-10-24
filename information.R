### this is the basic file defining the s.Cij function.
### The inputs should be:
### (1) number of treatment "t". some check should be made here.
### (2) number of period "p". some check should be made here.
### (3) corvariance matrix "sigma".
### (4) model index, denoted by "model" with model=0 as NULL:
# "0" for crossover design      # "1" for interference model
### (5) the effect to estimate, denoted by "effect" with effect=0 as NULL:
# "0" for direct effect         # "1" for total effect
### (6) the flag of circular/non-circular, denoted by "flag_circular" with effect=0 as NULL:
# "0" for non-circular          # "1" for circular
library(nloptr)
library(partitions)
library(MASS)
library(Matrix)
library(psych)
library(plyr)
library(permute)
library(combinat)
library(lpSolve)
s.Cij<-function(s,t,p,model=1,effect=1,flag_circular=0,sigma=diag(p)){
  if(model==1 || model==0)
  {
    ident=diag(t) # this is only used for constructional purpose
    sigma=diag(p) # covariance matrix
    #for (i in 1:p){
     # for(j in 1:p){
      #    sigma[i,j]=0.2^abs(i-j)
      #}
    #}
    smid=s[2:(length(s)-1)]
    sleft=s[1:(length(s)-2)]
    sright=s[3:length(s)]
    sigma_inverse=solve(sigma)
    v=sigma_inverse
    tb=v-v%*%matrix(1,p,p)%*%v/sum(v)
    bt=diag(1,t)-matrix(1/t,t,t)
    B.tilde=sigma_inverse-sigma_inverse%*%((rep(1,p))%*%t(rep(1,p)))%*%
      sigma_inverse/as.numeric(t(rep(1,p))%*%sigma_inverse%*%rep(1,p))
    tu=ident[smid,]
    lu=ident[sleft,]
    ru=ident[sright,]
    if(model==0){
      if(model==0&&effect==0){fu=lu}
      if(model==0&&effect==1){ld=lu-tu;fu=ld}
      coefb1=tr(bt%*%t(tu)%*%tb%*%tu%*%bt);coefb2=tr(bt%*%t(tu)%*%tb%*%fu%*%bt);coefb3=tr(bt%*%t(fu)%*%tb%*%fu%*%bt)
      C00=t(tu)%*%B.tilde%*%tu;C01=t(tu)%*%B.tilde%*%ld;C11=t(ld)%*%B.tilde%*%ld
      return(list(C00=C00,C01=C01,C11=C11,coefb1=coefb1,coefb2=coefb2,coefb3=coefb3))
    }
    if(model==1){
      if(model==1&&effect==0){fu=lu+ru}
      if(model==1&&effect==1){lu=lu-tu;ru=ru-tu;fu=lu+ru}
      coefb1=tr(bt%*%t(tu)%*%tb%*%tu%*%bt);coefb2=tr(bt%*%t(tu)%*%tb%*%fu%*%bt);coefb3=tr(bt%*%t(fu)%*%tb%*%fu%*%bt)
      C00=t(tu)%*%B.tilde%*%tu;C01=t(tu)%*%B.tilde%*%lu;C02=t(tu)%*%B.tilde%*%ru
      C11=t(lu)%*%B.tilde%*%lu;C12=t(lu)%*%B.tilde%*%ru;C22=t(ru)%*%B.tilde%*%ru
      return(list(C00=C00,C01=C01,C02=C02,C11=C11,C12=C12,C22=C22,coefb1=coefb1,coefb2=coefb2,coefb3=coefb3))
    } 
  }
  else{print("A self-defined model is used. Make sure you define the information matrix of each given sequence correctly.")}
}

### This program addresses the optimal sequence to be used for the final design. The inputs should be 
### (1) number of treatment "t". some check should be made here.
### (2) number of period "p". some check should be made here.
### (3) corvariance matrix "sigma".
### (4) model index, denoted by "model" with model=0 as NULL:
# "0" for crossover design      # "1" for interference model
### (5) the effect to estimate, denoted by "effect" with effect=0 as NULL:
# "0" for direct effect         # "1" for total effect
### (6) the flag of circular/non-circular, denoted by "flag_circular" with effect=0 as NULL:
# "0" for non-circular          # "1" for circular
### Packages to be used include:
source("information.R")
library(nloptr)
library(partitions)
library(MASS)
library(Matrix)
library(psych)
library(plyr)
library(permute)
library(combinat)
library(lpSolve)
#### define function "find_good_sequence" as follows

find_good_sequence<-function(p,t,sigma=diag(p),model=0,effect=0,flag_circular=0){
  error=1e-2
  sum.parts=restrictedparts(p+2,t)
  s.block=as.matrix(setparts(sum.parts))[,-1] 
  n.s=ncol(s.block)
  coefb=matrix(0,n.s,3)
  for(i in 1:n.s){
    temp=s.Cij(s.block[,i],t,p,model,effect,flag_circular)
    coefb[i,1]=temp$coefb1;coefb[i,2]=temp$coefb2;coefb[i,3]=temp$coefb3
  }
  maxfunc<-function(x){
    return(max(    coefb[,1]+2*coefb[,2]*x+coefb[,3]*x^2    ))
  }
  ymax=nlm(maxfunc,1)$minimum
  xach=nlm(maxfunc,1)$estimate
  indexset=NULL
  for(i in 1:n.s){
    temp=coefb[i,1]+2*coefb[i,2]*xach+coefb[i,3]*xach^2
    if(temp>(ymax-error)){indexset=append(indexset,i)}
  }
  final.s.block=s.block[,indexset]
    return(list(final.s.block=final.s.block,zstar=xach,ymax=ymax))
}

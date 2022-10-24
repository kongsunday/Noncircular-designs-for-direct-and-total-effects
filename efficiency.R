### This program addresses the exact design efficiency.
### The small subgroup of sequences to be used is given in the file "part_1.find_opt_sequences.R" and "part_2.find_exact_opt_design.R"
source("information.R")
source("sequence.R")
### The integer programming package to be used in this file is LpSolve.
library(lpSolve)
library(permute)
library(combinat)

effi_cal<-function(nn,designmatrix,opt.type="T",model=1,effect=1,flag_circular=0){
  p=nrow(designmatrix)-2
  n=ncol(designmatrix)
  t=max(as.vector(designmatrix))
  design.base=find_good_sequence(p,t,model = model,effect = effect,flag_circular = flag_circular)
  # The following functions defines the C_ij matrix.
  if(model==1){
    C_xi00=matrix(0,t,t);C_xi01=matrix(0,t,t);C_xi02=matrix(0,t,t);C_xi11=matrix(0,t,t);C_xi12=matrix(0,t,t);C_xi22=matrix(0,t,t)
    for(i in 1:n){
      Clist=s.Cij(designmatrix[,i],t,p,model=model,effect=effect,flag_circular = flag_circular)
      C_xi00=C_xi00+Clist$C00;C_xi01=C_xi01+Clist$C01;C_xi02=C_xi02+Clist$C02
      C_xi11=C_xi11+Clist$C11;C_xi12=C_xi12+Clist$C12;C_xi22=C_xi22+Clist$C22 }
    FINALC=C_xi00-cbind(C_xi01,C_xi02)%*%ginv(rbind(cbind(C_xi11,C_xi12),cbind(t(C_xi12),C_xi22)))%*%t(cbind(C_xi01,C_xi02))
  }
  if(model==0){
    C_xi00=matrix(0,t,t);  C_xi01=matrix(0,t,t);  C_xi11=matrix(0,t,t)
    for(i in 1:n){
      Clist=s.Cij(designmatrix[,i],t,p,model=model,effect=effect,flag_circular = flag_circular)
      C_xi00=C_xi00+Clist$C00;C_xi01=C_xi01+Clist$C01;C_xi11=C_xi11+Clist$C11 }
    FINALC=C_xi00-C_xi01%*%ginv(C_xi11)%*%t(C_xi01)
  }
  eig=(sort(eigen(FINALC)$values)[2:t])
  efficiency.A=(t-1)^2/nn/design.base$ymax/sum(  1/eig  )
  efficiency.D=(t-1)/nn/design.base$ymax*(prod(eig))^(1/(t-1))
  efficiency.E=(t-1)*min(eig)/nn/design.base$ymax
  efficiency.T=sum(eig)/nn/design.base$ymax
  if(opt.type=="E"){return(efficiency.E)}
  if(opt.type=="A"){return(efficiency.A)}
  if(opt.type=="D"){return(efficiency.D)}
  if(opt.type=="T"){return(efficiency.T)}
  if(opt.type=="ALL"){return(c(efficiency.E,efficiency.A,efficiency.D,efficiency.T))}
}

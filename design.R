### This program addresses the exact optimal design given the sequences to be used.
### The small subgroup of sequences to be used is given in the file "part_1.find_opt_sequences.R"
source("information.R")
source("sequence.R")
source("efficiency.R")
### The integer programming package to be used in this file is LpSolve.
library(lpSolve)
library(permute)
library(combinat)

### The function for finding optimal exact designs
find_good_design<-function(n=20,p=4,t=3,model=1,effect=1,flag_circular=0){
  ### the searching process
# covariance matrix
  bt=diag(1,t)-matrix(1/t,t,t)
  design.base=find_good_sequence(p,t,model=model,effect=effect,flag_circular = flag_circular)##more general function?
  num.col=ncol(design.base$final.s.block)
  if(length(num.col)==0){num.col=1}
  design.base.matrix=design.base$final.s.block
  temp.permu=matrix(unlist(permn(t)),nrow=t)
  full.Tau=NULL
  for(i in 1:num.col)
  {
    if(num.col==1){full.Tau=temp.permu[design.base$final.s.block,];break}
    full.Tau=cbind(full.Tau,temp.permu[design.base$final.s.block[,i],])
  }
  num.mat=length(s.Cij(full.Tau[,1],t,p,model = model,effect = effect,flag_circular = flag_circular))-3
  if(num.mat==3){
    all.matrix=apply(full.Tau,2,function(x){temp=s.Cij(x,t,p,model = model,effect = effect,flag_circular = flag_circular);
      ret1=temp[[2]]%*%bt;ret2=t(temp[[2]]);ret3=temp[[3]]%*%bt;
      return(c((temp[[1]]+design.base$zstar*ret1),(ret2+design.base$zstar*ret3)))})  
    final.restrict.1=cbind(    cbind(all.matrix[1:t^2,],-diag(t^2))  ,  matrix(0,t^2,t^2)    )
    final.restrict.2=cbind(    all.matrix[(t^2+1):(2*t^2),]  ,  matrix(0,t^2,t^2)  ,  -diag(t^2)    )
    final.restrict.3=cbind(    diag(rep(1,ncol(all.matrix)))  ,  matrix(0,ncol(all.matrix),2*t^2)    )
    final.restrict.t=rbind( final.restrict.1 , final.restrict.2 , c(rep(1,num.col*prod(1:t)),rep(0,2*t^2)) , final.restrict.3)
    restrict.vector=c(    n*as.vector(design.base$ymax*bt/(t-1)) , rep(0,t^2) , n , rep(0,ncol(all.matrix))    )
    Q.matrix=diag(   c(rep(0,(ncol(final.restrict.t)-2*t^2)),rep(1,t^2),rep(1,t^2))  )
    obj=rep(0,ncol(final.restrict.t))
    sense=c(  rep("=",(nrow(final.restrict.t)-ncol(all.matrix))) , rep(">=",ncol(all.matrix))   )
    v.type=c(   rep("I",(ncol(final.restrict.t)-2*t^2)) , rep("C",2*t^2)   )
    
    test<-list()
    
    test$Q<-Q.matrix    ### Q should be semi-positive definite
    test$obj<-obj
    
    test$A<-final.restrict.t
    test$rhs<-restrict.vector
    test$sense<-sense
    
    
    test$modelsense <- "min"
    test$vtype<-v.type
    test$lb=rep(-10,ncol(final.restrict.t))
    test$ub=rep(20,ncol(final.restrict.t))
    
    param=list(timeLimit = 1000)
    result=gurobi(test,param)
  }
  if(num.mat==6){
    all.matrix=apply(full.Tau,2,function(x){temp=s.Cij(x,t,p,model = model,effect = effect,flag_circular = flag_circular);
      ret1=(temp[[2]]+temp[[3]])%*%bt;ret2=rbind(t(temp[[2]]),t(temp[[3]]));ret3=rbind((temp[[4]]+temp[[5]])%*%bt,(t(temp[[5]])+temp[[6]])%*%bt);
      return(c((as.vector(temp[[1]])+design.base$zstar*as.vector(ret1)),(as.vector(ret2)+design.base$zstar*as.vector(ret3))))})  
    final.restrict.1=cbind(    cbind(all.matrix[1:t^2,],-diag(t^2))  ,  matrix(0,t^2,2*t^2)    )
    final.restrict.2=cbind(    all.matrix[(t^2+1):(3*t^2),]  ,  matrix(0,2*t^2,t^2)  ,  -diag(2*t^2)   )
    final.restrict.3=cbind(    diag(rep(1,ncol(all.matrix)))  ,  matrix(0,ncol(all.matrix),3*t^2)    )
    final.restrict.t=rbind( final.restrict.1 , final.restrict.2 , c(rep(1,num.col*prod(1:t)),rep(0,3*t^2)) , final.restrict.3)
    restrict.vector=c(    n*as.vector(design.base$ymax*bt/(t-1)) , rep(0,2*t^2) , n , rep(0,ncol(all.matrix))    )
    Q.matrix=diag(   c(rep(0,(ncol(final.restrict.t)-3*t^2)),rep(1,t^2),rep(1,2*t^2))  )
    obj=rep(0,ncol(final.restrict.t))
    sense=c(  rep("=",(nrow(final.restrict.t)-ncol(all.matrix))) , rep(">=",ncol(all.matrix))   )
    v.type=c(   rep("I",(ncol(final.restrict.t)-3*t^2)) , rep("C",3*t^2)   )
    
    test<-list()
    
    test$Q<-Q.matrix    ### Q should be semi-positive definite
    test$obj<-obj
    
    test$A<-final.restrict.t
    test$rhs<-restrict.vector
    test$sense<-sense
    
    
    test$modelsense <- "min"
    test$vtype<-v.type
    test$lb=rep(-10,ncol(final.restrict.t))
    test$ub=rep(20,ncol(final.restrict.t))
    
    param=list(timeLimit = 1000)
    result=gurobi(test,param)
  }
  index=as.integer(result$x[1:(num.col*prod(1:t))]+0.1)
  final.design=full.Tau[,rep((1:length(index)),index)]
  ef=effi_cal(n,final.design,opt.type = "ALL",model=model,effect=effect,flag_circular = flag_circular)
  print(num.mat)
  return(list(final.design=final.design,ef=ef))
}


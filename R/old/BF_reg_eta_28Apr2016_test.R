# revised based on BF_25Feb216.R

BF_reg_eta = function(Data,par.set,KK,nsites){
  
  eta.star.reg = par.set[1] 
  k.star.reg = par.set[2] 
  eta1.star.reg = par.set[3] 
  k1.star.reg = par.set[4] 
  eta0.star.reg = par.set[5] 
  k0.star.reg = par.set[6] 
  
  bfexch_case = function(theta,datapar){
    x = datapar$data
    n = nsites
    K = datapar$K
    eta.pos = exp(theta)/(1+exp(theta))
    N = length(x)
    alpha.star = eta1.star.reg*k1.star.reg
    beta.star = k1.star.reg*(1-eta1.star.reg)
    z= 0*theta;
    for(i in 1:N)
      z=z+lbeta(K*eta.pos+x[i],K*(1-eta.pos)+n-x[i])
    z=z-N*lbeta(K*eta.pos,K*(1-eta.pos))+log(eta.pos*(1-eta.pos))
    z = z+dbeta(eta.pos,alpha.star,beta.star,log=TRUE)
    return(z)
  }
  
  bfexch_control = function(theta,datapar){
    x = datapar$data
    n = nsites
    K = datapar$K
    eta.pos = exp(theta)/(1+exp(theta))
    N = length(x) 
    alpha.star = eta0.star.reg*k0.star.reg
    beta.star = k0.star.reg*(1-eta0.star.reg)
    z= 0*theta;
    for(i in 1:N)
      z=z+lbeta(K*eta.pos+x[i],K*(1-eta.pos)+n-x[i])
    z=z-N*lbeta(K*eta.pos,K*(1-eta.pos))+log(eta.pos*(1-eta.pos))
    z = z+dbeta(eta.pos,alpha.star,beta.star,log=TRUE)
    return(z)
  }
  
  bfexch_total = function(theta,datapar){
    x = datapar$data
    n = nsites
    K = datapar$K
    eta.pos = exp(theta)/(1+exp(theta))
    N = length(x)
    alpha.star = eta.star.reg*k.star.reg
    beta.star = k.star.reg*(1-eta.star.reg)
    z= 0*theta;
    for(i in 1:N)
      z=z+lbeta(K*eta.pos+x[i],K*(1-eta.pos)+n-x[i])
    z=z-N*lbeta(K*eta.pos,K*(1-eta.pos))+log(eta.pos*(1-eta.pos))
    z = z+dbeta(eta.pos,alpha.star,beta.star,log=TRUE)
    return(z)
  }
  
#   library(LearnBayes)
  
  # BF function based on non-informative prior, uniform prior
  BF_function = function(K0,obs.data){
    N = nrow(obs.data)
    cases = which(obs.data[,2] == 1)
    controls = which(obs.data[,2] == 0)
    
    #cat(cases,"\n")
    #cat(controls,"\n")
    #cat("par=",par.set,"\n")
    
    
    logm_total = laplace(bfexch_total,0,list(data=obs.data[,1],K=K0))$int
    logm_case = laplace(bfexch_case,0,list(data=obs.data[cases,1],K=K0))$int
    #logm_case
    logm_control = laplace(bfexch_control,0,list(data=obs.data[controls,1],K=K0))$int
    #logm_control
    print(c(logm_case , logm_control , logm_total))
    return(exp(logm_case + logm_control - logm_total))
  }
  
  return(BF_function(KK,Data))
}
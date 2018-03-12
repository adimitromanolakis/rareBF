# revised based on BF_25Feb216.R


#' calculate BF with beta prior
#'
#' @param obs.data N*2 data matrix, N is the sample size, 1st column is sum of rare variants for each individual, 2nd column is the disease status
#'
#' @return Vector of 2 elements. 1st: precision parameter of beta distribution; 2nd: BF
#' @export
#' 
# rename function to BFbeta
BFbeta = function(obs.data)  {
  
  
  
  #print(str(Data))
  
  #obs.data = Data
  
  bfexch = function(par.vec,data.x){
    
    theta = par.vec[1]
    K = par.vec[2]
    x = data.x
    n = length(causal.pool)
    eta.pos = exp(theta)/(1+exp(theta))
    N = length(x)
    # alpha.star = eta.star.reg*k.star.reg
    # beta.star = k.star.reg*(1-eta.star.reg)
    z= 0
    for(i in 1:N)
      z=z+lbeta(K*eta.pos+x[i],K*(1-eta.pos)+n-x[i])
    z=z-N*lbeta(K*eta.pos,K*(1-eta.pos))
    # z = z+dbeta(eta.pos,alpha.star,beta.star,log=TRUE)
    return(z)
  }
  
  # beta.root = optim(c(0,500),bfexch,data.x=obs.data[,1],method="L-BFGS-B",lower=c(-20,10),upper=c(0.0001,100000),control=list(fnscale=-1))$par
  beta.root = optim(c(0,200),bfexch,data.x=obs.data[,1],method="L-BFGS-B",lower=c(-20,10),upper=c(0,100000),control=list(fnscale=-1))$par
  
  KK = beta.root[2]
  
  bfexch_all = function(theta,datapar){
    x = datapar$data
    n = length(causal.pool)
    K = datapar$K
    eta.pos = exp(theta)/(1+exp(theta))
    N = length(x)
    # alpha.star = eta.star.reg*k.star.reg
    # beta.star = k.star.reg*(1-eta.star.reg)
    z= 0
    for(i in 1:N)
      z=z+lbeta(K*eta.pos+x[i],K*(1-eta.pos)+n-x[i])
    z=z-N*lbeta(K*eta.pos,K*(1-eta.pos))
    # z = z+dbeta(eta.pos,alpha.star,beta.star,log=TRUE)
    return(z)
  }
  

    
    # N = nrow(obs.data)
    logm_total = optim(0,bfexch_all,gr=NULL,list(data=obs.data[,1],K=KK),hessian = T,control=list(fnscale=-1),method="Brent",lower=-100,upper = 0)
    eta_total = exp(logm_total$par)/(1+exp(logm_total$par))
    var_total = -solve(logm_total$hessian) * (exp(logm_total$par)/(1+exp(logm_total$par))^2)^2
  
    logm_case = optim(0,bfexch_all,gr=NULL,list(data=obs.data[1:sum(obs.data[,2]),1],K=KK),hessian = T,control=list(fnscale=-1),method="Brent",lower=-100,upper = 0)
    eta_case = exp(logm_case$par)/(1+exp(logm_case$par))
    var_case = -solve(logm_case$hessian) * (exp(logm_case$par)/(1+exp(logm_case$par))^2)^2
    
    logm_control = optim(0,bfexch_all,gr=NULL,list(data=obs.data[(sum(obs.data[,2])+1):nrow(obs.data),1],K=KK),hessian = T,control=list(fnscale=-1),method="Brent",lower=-100,upper = 0)
    eta_control = exp(logm_control$par)/(1+exp(logm_control$par))
    var_control = -solve(logm_control$hessian) * (exp(logm_control$par)/(1+exp(logm_control$par))^2)^2
    # print(c(eta_total,var_total,eta_case,var_case,eta_control,var_control))

    k.star.reg = eta_total/var_total
    alpha.star = eta_total*k.star.reg
    beta.star = k.star.reg*(1-eta_total)
    I_total = log(var_total)/2 + (log(2*pi)/2)+ logm_total$value + dbeta(eta_total,alpha.star,beta.star,log=TRUE)
    # I_total_1 = bfexch_case(logm_total$par,list(data=obs.data[1:sum(obs.data[,2]),1],K=K0))
    # I_total_0 = bfexch_control(logm_total$par,list(data=obs.data[(sum(obs.data[,2])+1):nrow(obs.data),1],K=K0))
    #   
    k1.star.reg = eta_case/var_case
    alpha.star = eta_case*k1.star.reg
    beta.star = k1.star.reg*(1-eta_case)
    I_case = log(var_case)/2 + (log(2*pi)/2)+ logm_case$value + dbeta(eta_case,alpha.star,beta.star,log=TRUE)
    # 
    k0.star.reg = eta_control/var_control
    alpha.star = eta_control*k0.star.reg
    beta.star = k0.star.reg*(1-eta_control)
    I_control = log(var_control)/2 + (log(2*pi)/2)+ logm_control$value +dbeta(eta_control,alpha.star,beta.star,log=TRUE)
    
  #return(c(KK,eta_total,eta_case,eta_control,logm_total$value,logm_case$value,logm_control$value,exp(I_case+I_control-I_total)))
    # return(c(logm_case$value-I_total_1,logm_control$value-I_total_0))
    
    
    vec = c(KK,exp(I_case+I_control-I_total))
    
    names(vec) = c("KK","BF")
    return(vec)
    
    
    
}

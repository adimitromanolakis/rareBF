# # The BF is based on minor allele counts 
# mixed prior for p, compare both w0 and eta

BF_mix_both = function(Data,par.set){
  
  eta.star=par.set[1]
  k.star=par.set[2]
  eta1.star=par.set[3]
  k1.star=par.set[4]
  eta0.star=par.set[5]
  k0.star=par.set[6]
  eta.tilde=par.set[7]
  k.tilde=par.set[8]
  eta1.tilde=par.set[9]
  k1.tilde=par.set[10]
  eta0.tilde=par.set[11]
  k0.tilde=par.set[12]
  K=par.set[13]
    
  bfexch_case = function(theta,datapar){
    x = sort(datapar$data)
    n = length(causal.pool)
    m0 = datapar$m0
    j0 = datapar$j0
    eta.pos = exp(theta)/(1+exp(theta))
    N = length(x)
    alpha.star = eta1.star*k1.star
    beta.star = k1.star*(1-eta1.star)
    z= 0*theta
    z = z + (m0-j0)*lbeta(K*eta.pos,K*(1-eta.pos)+n) - (N-j0)*lbeta(K*eta.pos,K*(1-eta.pos)) 
    for(s in (m0+1):N)
      z=z+lbeta(K*eta.pos+x[s],K*(1-eta.pos)+n-x[s])
    z=z+log(eta.pos*(1-eta.pos))
    z = z+dbeta(eta.pos,alpha.star,beta.star,log=TRUE)
    return(z)
    #return(c(m0,N))
  }
  
  bfexch_control = function(theta,datapar){
    x = sort(datapar$data)
    n = length(causal.pool)
    m0 = datapar$m0
    j0 = datapar$j0
    eta.pos = exp(theta)/(1+exp(theta))
    N = length(x) 
    alpha.star = eta0.star*k0.star
    beta.star = k0.star*(1-eta0.star)
    z= 0*theta;
    z = z + (m0-j0)*lbeta(K*eta.pos,K*(1-eta.pos)+n) - (N-j0)*lbeta(K*eta.pos,K*(1-eta.pos)) 
    for(s in (m0+1):N)
      z=z+lbeta(K*eta.pos+x[s],K*(1-eta.pos)+n-x[s])
    z=z+log(eta.pos*(1-eta.pos))
    z = z+dbeta(eta.pos,alpha.star,beta.star,log=TRUE)
    return(z)
  }
  
  bfexch_total = function(theta,datapar){
    x = sort(datapar$data)
    n = length(causal.pool)
    
    m0 = datapar$m0
    j0 = datapar$j0
    eta.pos = exp(theta)/(1+exp(theta))
    N = length(x)
    alpha.star = eta.star*k.star
    beta.star = k.star*(1-eta.star)
    z= 0*theta;
    z = z + (m0-j0)*lbeta(K*eta.pos,K*(1-eta.pos)+n) - (N-j0)*lbeta(K*eta.pos,K*(1-eta.pos)) 
    for(s in (m0+1):N)
      z=z+lbeta(K*eta.pos+x[s],K*(1-eta.pos)+n-x[s])
    z=z+log(eta.pos*(1-eta.pos))
    z = z+dbeta(eta.pos,alpha.star,beta.star,log=TRUE)
    return(z)
  }
  
#   library(LearnBayes)
  
  T.const = function(m.pos,DATA){
    t.const = rep(NA,m.pos+1)
    for(l in 0:m.pos){
      
      t.const[l+1] = lchoose(m.pos,l) + lbeta(l+eta.tilde*k.tilde,nrow(DATA)-l+k.tilde*(1-eta.tilde)) + I_fun(m.pos,l,"total",DATA)
    }
    #   return(list(mean(t.const),pos.indx = which(t.const>mean(t.const))))
    return(mean(t.const))
  }
  
  sum_function = function(m.pos,type,T_total,DATA){
    S = 0
    if(type=="total"){
      for(l in 0:m.pos){
        
        S = S + exp(lchoose(m.pos,l) + lbeta(l+eta.tilde*k.tilde,nrow(DATA)-l+k.tilde*(1-eta.tilde)) + I_fun(m.pos,l,type,DATA)-T_total)
        print(lchoose(m.pos,l) + lbeta(l+eta.tilde*k.tilde,nrow(DATA)-l+k.tilde*(1-eta.tilde)) + I_fun(m.pos,l,type,DATA)-T_total)
      }
      return(log(S) - lbeta(eta.tilde*k.tilde,k.tilde*(1-eta.tilde)))
    }else if(type == "case"){
      for(l in 0:m.pos){
        
        S = S + exp(lchoose(m.pos,l) + lbeta(l+eta1.tilde*k1.tilde,nrow(DATA)-l+k1.tilde*(1-eta1.tilde)) + I_fun(m.pos,l,type,DATA)-(T_total/2))
        print(lchoose(m.pos,l) + lbeta(l+eta1.tilde*k1.tilde,nrow(DATA)-l+k1.tilde*(1-eta1.tilde)) + I_fun(m.pos,l,type,DATA)-(T_total/2))
      }
      return(log(S) - lbeta(eta1.tilde*k1.tilde,k1.tilde*(1-eta1.tilde)))
    }
    else{
      for(l in 0:m.pos){
        
        S = S + exp(lchoose(m.pos,l) + lbeta(l+eta0.tilde*k0.tilde,nrow(DATA)-l+k0.tilde*(1-eta0.tilde)) + I_fun(m.pos,l,type,DATA)-(T_total/2))
        print(lchoose(m.pos,l) + lbeta(l+eta0.tilde*k0.tilde,nrow(DATA)-l+k0.tilde*(1-eta0.tilde)) + I_fun(m.pos,l,type,DATA)-(T_total/2))
      }
      return(log(S) - lbeta(eta0.tilde*k0.tilde,k0.tilde*(1-eta0.tilde)))
    }
  }
  
  I_fun = function(mm.pos,j.pos,type,dataset){
    if(type == "total")
      logm = laplace(bfexch_total,0,list(data=dataset[,1],m0 = mm.pos,j0 = j.pos))$int
    else if(type == "case")
      logm = laplace(bfexch_case,0,list(data=dataset[,1],m0 = mm.pos,j0 = j.pos))$int
    else logm = laplace(bfexch_control,0,list(data=dataset[,1],m0 = mm.pos,j0 = j.pos))$int
    return(logm)
  }
  
  log_BF = function(obs.data){
    
    obs.data1 = obs.data[1:sum(obs.data[,2]),]
    obs.data0 = obs.data[-(1:sum(obs.data[,2])),]
    
    m.total = sum(obs.data[,1]==0)
    m.case = sum(obs.data1[,1]==0)
    m.control = sum(obs.data0[,1]==0)
    
    
    TT_total = T.const(m.total,obs.data)
    
    sum_case = sum_function(m.case, "case",TT_total,obs.data1)
    # if I don't add obs.data as argument of the function, the result can not be updated for each simulated dataset
    sum_control = sum_function(m.control, "control",TT_total,obs.data0 )
    sum_total = sum_function(m.total, "total",TT_total,obs.data)
    print(c(sum_case,sum_control,sum_total))
    
    return(sum_case+sum_control-sum_total)
  }
  return(exp(log_BF(Data)))
}



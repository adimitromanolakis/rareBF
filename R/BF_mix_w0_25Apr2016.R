# # mixed prior for p, compare w0 between cases and controls

BF_mix_w0 = function(Data,par.set){
  
  z.function = function(eta.pos,k.pos){
    return(beta(eta.pos*k.pos,(length(causal.pool)+k.pos*(1-eta.pos)))/beta(eta.pos*k.pos,(k.pos*(1-eta.pos))))
  }
  
  I_function = function(m.pos,eta.s,k.s,N,tt,Z){
    S = 0
    for(l in 0:m.pos){
      S = S + exp(lchoose(m.pos,l) + (m.pos-l)*log(Z) + lbeta((l + eta.s*k.s), (-l+N+k.s*(1-eta.s))) - tt)
    }
    return(log(S))
  }
  
  t_function = function(m.pos,eta.s,k.s,N,Z){
    t_sum = rep(NA,m.pos+1)
    for(l in 0:m.pos){
      t_sum[l+1] = lchoose(m.pos,l) + (m.pos-l)*log(Z) + lbeta((l + eta.s*k.s), (-l+N+k.s*(1-eta.s))) 
    }
    return(mean(t_sum))
  }
    
  C_function = function(eta.s,k.s){
    S = 0
    S = S - lbeta((eta.s*k.s),(k.s*(1-eta.s)))
    return(S) 
  }
  
  eta.tilde = par.set[1]
  k.tilde = par.set[2]
  eta1.tilde = par.set[3]
  k1.tilde = par.set[4]
  eta0.tilde = par.set[5]
  k0.tilde = par.set[6]
  eta.par = par.set[7]
  k.par = par.set[8]
  
  
  log_BF = function(obs.data){
    obs.data1 = obs.data[1:sum(obs.data[,2]),]
    obs.data0 = obs.data[-(1:sum(obs.data[,2])),]
    
    
    m.total = sum(obs.data[,1]==0)
    m.case = sum(obs.data1[,1]==0)
    m.control = sum(obs.data0[,1]==0)
    
    C_case = C_function(eta1.tilde,k1.tilde)
    C_control = C_function(eta0.tilde,k0.tilde)
    C_total = C_function(eta.tilde,k.tilde)
    
    t = t_function(m.total,eta.tilde,k.tilde,nrow(obs.data),z.function(eta.par,k.par))/2
    
    I_case = I_function(m.case, eta1.tilde, k1.tilde, nrow(obs.data1), t,z.function(eta.par,k.par))
    I_control = I_function(m.control,eta0.tilde,k0.tilde,nrow(obs.data0),t,z.function(eta.par,k.par))
    I_total = I_function(m.total,eta.tilde,k.tilde,nrow(obs.data),2*t,z.function(eta.par,k.par))
    
    return(C_case + C_control -C_total + I_case + I_control - I_total)
  }
  return(exp(log_BF(Data)))
         
}

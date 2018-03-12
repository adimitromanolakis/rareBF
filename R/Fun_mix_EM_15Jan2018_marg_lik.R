# revised based on BF_25Feb216.R

## Rename to BFmixture


#' calculate BF with mixture prior.
#'
#' @param obs.data N*2 data matrix, N is the sample size, 1st column is sum of rare variants for each individual, 2nd column is the disease status
#' @param nvariants number of sites of the region
#' @param low.bound lower bound
#'
#' @return Vector of 3 elements. 1st, probability of p=0; 2nd, precision parameter of beta distribution,3rd: BF
#' 
#' @examples 
#' 
#'nsites = 500
#'ncase = 200
#'ncontrol = 200
#'maf = 0.005
#'rexp_param = 2
#'
#'model1 = function(i) {
#'  
#'  p1 = rbinom( ncontrol  , nsites , maf   )
#'  p2 = rbinom( ncase     , nsites , maf )  + H1 * round(rexp(ncase,rexp_param))
  
#'  disease_status = c (  rep(0,ncontrol),rep(1,ncase)  )
#'  data = data.frame(c(p1,p2), disease_status)
  
#'  bf = BFmixture(data, nsites)[3]
#'  bf
#'}
#'
#'nrep = 20
#'
#'H1 = 0
#'bf1 = ( replicate(nrep,model1() ) )
#'
#'H1 = 1
#'bf2 = ( replicate(nrep,model1() ) )
#'
#'col = c( rep("blue",nrep) , rep("red", nrep))
#'
#'
#'plot(log(c(bf1,bf2)),col = col ,pch=20,cex=1.5)
#' 
#' @export
#' 
BFmixture = function(obs.data, nvariants, low.bound = 10^(-2) ) {
  # 
  converge.ind = 0 # if converge.ind is equal to 1, meaning EM estimate for w0, eta and k doesn't converge
  obs.data1 = obs.data[1:sum(obs.data[,2]),]
  obs.data2 = obs.data[-(1:sum(obs.data[,2])),]
  
  m = sum(obs.data[,1]==0)
  m1 = sum(obs.data1[,1]==0)
  m2 = sum(obs.data2[,1]==0)
  
  x = obs.data[,1]
  x = sort(x)
  N = length(x)
  
  x1 = obs.data1[,1]
  x1 = sort(x1)
  N1 = length(x1)
  
  x2 = obs.data2[,1]
  x2 = sort(x2)
  N2 = length(x2)
  
  obs.data$p.hat = obs.data[,1]/nvariants
  w0.initial = sum(obs.data[,1]==0)/N
  # eta.initial = mean(c(obs.data$p.hat[which(obs.data[,1]!=0)],rep(0,m-(w0.initial*N))))
  # k.initial = eta.initial*(1-eta.initial)/var(c(obs.data$p.hat[which(obs.data[,1]!=0)],rep(0,m-(w0.initial*N))))-1
  
  eta.initial = mean(obs.data$p.hat[which(obs.data[,1]!=0)])
  if(var(obs.data$p.hat[which(obs.data[,1]!=0)]) == 0){
    return(rep(NA,9))
  }else{
    k.initial = eta.initial*(1-eta.initial)/var(obs.data$p.hat[which(obs.data[,1]!=0)])-1
    
    # like_all_eta_k = function(par.est,w0){
    #   eta.pos = exp(par.est[1])/(1+exp(par.est[1]))
    #   k.pos = par.est[2]
    # 
    #   z=0
    #   # z = z + m.pos * log(w0 + (1-w0)*beta(eta.pos*k.pos,(k.pos*(1-eta.pos)+n))/beta(eta.pos*k.pos,(k.pos*(1-eta.pos))))
    #   z = z + m * log(w0 + (1-w0)*exp(lbeta(eta.pos*k.pos,(k.pos*(1-eta.pos)+nvariants))-lbeta(eta.pos*k.pos,(k.pos*(1-eta.pos)))))+log(eta.pos*(1-eta.pos))
    #   for(l in (m+1):N){
    #     z = z + log(1-w0) + lbeta((eta.pos*k.pos+x[l]),(k.pos*(1-eta.pos)+nvariants-x[l])) - lbeta(eta.pos*k.pos,k.pos*(1-eta.pos))
    #   }
    #   return(z)
    # }
    # 
    # eta.par.est = optim(c(log(eta.initial/(1-eta.initial)),k.initial),like_all_eta_k,w0.initial,gr=NULL,control=list(fnscale=-1))$par
    # eta.initial = exp(eta.par.est[1])/(1+exp(eta.par.est[1]))
    # k.initial = eta.par.est[2] # from what we observed, k.hat fluctuate very much across the simulated datasets
    
    #library(rootSolve)
    count = 1
    par.set = matrix(NA,nrow=200,ncol=3)
    par.set[1,1] = w0.initial
    par.set[1,2] = eta.initial
    par.set[1,3] = k.initial
    delta.cut = 0.01
    delta.par = c(1,1,1)
    
    # n = nvariants
    while(max(delta.par)>delta.cut){
      count = count+1
      # for(count in 2:100){
      exp.z = (1-par.set[count-1,1])*exp(lbeta(par.set[count-1,2]*par.set[count-1,3],(par.set[count-1,3]*(1-par.set[count-1,2])+nvariants))-lbeta(par.set[count-1,2]*par.set[count-1,3],par.set[count-1,3]*(1-par.set[count-1,2])))
      exp.z = exp.z/(par.set[count-1,1]+exp.z)
      m.pi =  m*(1-exp.z)
      par.set[count,1] = m.pi/N
      
      loglike_der = function(beta.par){
        eta = exp(beta.par[1])/(1+exp(beta.par[1]))
        # eta = beta.par[1]
        k = beta.par[2]
        z = (m-m.pi)*(lbeta(k*eta,k*(1-eta)+nvariants)-lbeta(k*eta,k*(1-eta))) - (N-m)*lbeta(k*eta,k*(1-eta))
        for(i in (m+1):N){
          z = z+ lbeta(k*eta+x[i],k*(1-eta)+nvariants-x[i])
        }
        return(z)
      }
      #beta.root = optim(c(par.set[1,2],par.set[1,3]),loglike_der,method="L-BFGS-B",lower=c(0.0001,10),upper=c(0.5,100000),control=list(fnscale=-1))$par
      beta.root = optim(c(log(par.set[count-1,2]/(1-par.set[count-1,2])),par.set[count-1,3]),loglike_der,method="L-BFGS-B",lower=c(-10,10),upper=c(0,100000),control=list(fnscale=-1))$par
      # beta.root = optim(c(log(eta.initial/(1-eta.initial)),k.initial),loglike_der,method="L-BFGS-B",lower=c(0.0001,10),upper=c(0.5,100000),control=list(fnscale=-1))$par
      
      par.set[count,2] = exp(beta.root[1])/(1+exp(beta.root[1]))
      # par.set[count,2] = beta.root[1]
      par.set[count,3] = beta.root[2]
      
      # print(c(count,par.set[count,]))
      
      delta.par = abs(par.set[count,]-par.set[count-1,])/par.set[count-1,]
      
      if(count == nrow(par.set)| par.set[count,1]<=low.bound){
        delta.par = c(0.01,0.01,0.01)
        converge.ind=1
      }
    }
    
    w0.fix = par.set[count,1]
    k.fix = par.set[count,3]
    
    # given w0 and k, estimate eta and its std. error using univariate EM algorithm
    Eta_EM = function(sample_all, data.x, sample_m){
      par.eta = rep(NA,50)
      par.eta[1] = par.set[count,2]
      delta.cut = 0.001
      delta.par = 1
      
      count0 = 1
      while(delta.par>delta.cut){
        count0 = count0+1
        # for(count in 2:100){
        exp.z = (1-w0.fix)*exp(lbeta(par.eta[count0-1]*k.fix,(k.fix*(1-par.eta[count0-1])+nvariants))-lbeta(par.eta[count0-1]*k.fix,k.fix*(1-par.eta[count0-1])))
        exp.z = w0.fix/(w0.fix+exp.z)
        m.pi =  sample_m*exp.z
        
        # loglike_der_eta = function(xx){
        #   eta = xx
        #   k = k.fix
        #   
        #   # w.r.t. eta
        #   zz = (sample_m-m.pi)*(-digamma(k*(1-eta)+nvariants) + digamma(k*(1-eta))) - (sample_all-sample_m)*(digamma(k*eta) - digamma(k*(1-eta)))
        #   for(i in (sample_m+1):sample_all){
        #     zz = zz + digamma(k*eta+data.x[i]) - digamma(k*(1-eta)+nvariants-data.x[i])
        #   }
        #   return(zz)
        # }
        # eta.root = uniroot(f=loglike_der_eta,c(0.0001,0.5))$root
        
        loglike_der_eta = function(xx){
          eta = exp(xx)/(1+exp(xx))
          k = k.fix
          z = (sample_m-m.pi)*(lbeta(k*eta,k*(1-eta)+nvariants)-lbeta(k*eta,k*(1-eta))) - (sample_all-sample_m)*lbeta(k*eta,k*(1-eta))
          for(i in (sample_m+1):sample_all){
            z = z+ lbeta(k*eta+data.x[i],k*(1-eta)+nvariants-data.x[i])
          }
          return(z)
        }
        
        eta.root = optim(log(par.eta[count0-1]/(1-par.eta[count0-1])),loglike_der_eta,method="Brent",lower=-10,upper=0,control=list(fnscale=-1))$par
        
        par.eta[count0] = exp(eta.root)/(1+exp(eta.root))
        # print(c(count0,par.eta[count0]))
        delta.par = abs(par.eta[count0]-par.eta[count0-1])/par.eta[count0-1]
      }
      
      eta.em = par.eta[count0]
      # Louis formula
      
      exp_z = (1-w0.fix)*exp(lbeta(eta.em*k.fix,(k.fix*(1-eta.em)+nvariants))-lbeta(eta.em*k.fix,k.fix*(1-eta.em)))
      exp_z = w0.fix/(exp_z+w0.fix)
      m_pi =  sample_m*exp_z
      
      i_com = (sample_m-m_pi)*k.fix^2*( trigamma(k.fix*(1-eta.em)+nvariants) - trigamma(k.fix*(1-eta.em))) - (sample_all-sample_m)*k.fix^2*(trigamma(k.fix*eta.em) + trigamma(k.fix*(1-eta.em)))
      for(i in (sample_m+1):sample_all){
        i_com = i_com + k.fix^2*(trigamma(k.fix*eta.em+data.x[i]) + trigamma(k.fix*(1-eta.em)+nvariants- data.x[i]))
      }
      i_com = -i_com
      i_miss = (k.fix*( - digamma(k.fix*(1-eta.em)+nvariants) + digamma(k.fix*(1-eta.em))))^2 * sample_m * (exp_z - exp_z^2)
      
      return(c(eta.em,  1/(i_com-i_miss)))
    }
    
    est_total = Eta_EM(N,x,m)
    est_case = Eta_EM(N1,x1,m1)
    est_control = Eta_EM(N2,x2,m2)
    
    loglike_obs = function(sample_all,data.x,sample_m,eta){
      z = sample_m*log(w0.fix+(1-w0.fix)*exp(lbeta(eta*k.fix,k.fix*(1-eta)+nvariants)-lbeta(eta*k.fix,k.fix*(1-eta))))
      for(i in (sample_m+1):sample_all){
        z = z + log(1-w0.fix) + lchoose(nvariants,data.x[i]) + lbeta(eta*k.fix+data.x[i],k.fix*(1-eta)+nvariants-data.x[i]) - lbeta(eta*k.fix,k.fix*(1-eta))
      }
      return(z)
    }
    
    k.star.reg = est_total[1]/est_total[2]
    alpha.star = est_total[1]*k.star.reg
    beta.star = k.star.reg*(1-est_total[1])
    I_total = log(est_total[2])/2 + (log(2*pi)/2)+ loglike_obs(N,x,m,est_total[1]) + dbeta(est_total[1],alpha.star,beta.star,log=TRUE)
    
    #   
    k1.star.reg =  est_case[1]/est_case[2]
    alpha.star = est_case[1]*k1.star.reg
    beta.star = k1.star.reg*(1-est_case[1])
    I_case = log(est_case[2])/2 + (log(2*pi)/2)+ loglike_obs(N1,x1,m1,est_case[1]) + dbeta(est_case[1],alpha.star,beta.star,log=TRUE)
    
    # 
    k0.star.reg = est_control[1]/est_control[2]
    alpha.star = est_control[1]*k0.star.reg
    beta.star = k0.star.reg*(1-est_control[1])
    I_control = log(est_control[2])/2 + (log(2*pi)/2)+ loglike_obs(N2,x2,m2,est_control[1]) +dbeta(est_control[1],alpha.star,beta.star,log=TRUE)
    
    # return(c(w0.fix,est_total,est_case,est_control,loglike_obs(N,x,m,est_total[1]),loglike_obs(N1,x1,m1,est_case[1]),loglike_obs(N2,x2,m2,est_control[1]),exp(I_case+I_control-I_total)))
    #return(c(w0.fix,k.fix,est_total[1],est_case[1],est_control[1],loglike_obs(N,x,m,est_total[1]),loglike_obs(N1,x1,m1,est_case[1]),loglike_obs(N2,x2,m2,est_control[1]),exp(I_case+I_control-I_total)))
    
    vec = c(w0.fix,k.fix,exp(I_case+I_control-I_total))
    
    names(vec) = c("prob.p.eq.0", "precision","BF")
    
    return( vec )
    
  }
  
}
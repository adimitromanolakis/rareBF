
library(sim1000G)
library(rareBF)


# total sample size, assuming equal sample sizes for cases and controls
ss = 500 

# number of sites within the region
nvariants = 72 


# number of causal variants under H1
p = 8 



filename = "region-chr4-93-TMEM156.vcf.gz"




convert.gt = function(xx) {
  
  xx = ifelse(xx==2,1,xx)
  
  return(xx)
}


# function to simulate data under null hypothesis

data_sim_H0 = function(vcf_sub,seed.num)  {

    startSimulation(vcf_sub, totalNumberOfIndividuals = sample.size) 
  
    SIM$reset()  
    
    set.seed(seed.num)
    
    id = generateUnrelatedIndividuals(sample.size)
    gt = retrieveGenotypes(id)
    gt = matrix(unlist(lapply(gt,convert.gt)),sample.size,nvariants)
    
    freq = apply(gt,2,sum)/(2*nrow(gt))
    
    if(sum(freq>=0.5)>0) {
      
      for(jj in which(freq>=0.5)) {
        gt[,jj] = 1-gt[,jj]
      }
    }
    
    
    if(sum(freq>=0.05)>0)  {
      gt = gt[,-which(freq>=0.05)]
    }
  
    set.seed(seed.num)
    genocase = sample(sample.size,ss)
    
    
    return(list(geno = gt[genocase,],nvar = ncol(gt)))
    
}

# function to simulate data under alternative hypothesis
data_sim = function(vcf_sub,seed.num) {
  
  startSimulation(vcf_sub, totalNumberOfIndividuals = sample.size) 
  SIM$reset() 
  
  set.seed(seed.num)
  
  id = generateUnrelatedIndividuals(sample.size)
  gt = retrieveGenotypes(id)
  gt = matrix(unlist(lapply(gt,convert.gt)),sample.size,nvariants)
  
  freq = apply(gt,2,sum)/(2*nrow(gt))
  
  if(sum(freq>=0.5)>0)  {
    for(jj in which(freq>=0.5))  {
      gt[,jj] = 1-gt[,jj]
    }
  }
  
  
  if(sum(freq>=0.05)>0)  {
    gt = gt[,-which(freq>=0.05)]
  }
  
  
  freq = apply(gt,2,sum)/(2*nrow(gt))
  causal = sample(setdiff(1:ncol(gt),which(freq==0)),p)
  
  beta.sign = rep(1,p)
  c.value = 0.402
  beta.abs = c.value*abs(log10(freq[causal]))
  beta.val = beta.sign*beta.abs
  x.bar = apply(gt[,causal],2,mean)
  x.bar = as.matrix(x.bar)
  beta.val = t(as.matrix(beta.val))
  beta0 = 0-beta.val %*% x.bar
  
  eta = beta.val %*% t(gt[,causal])
  eta = as.vector(eta) + rep(beta0,nrow(gt))
  prob = exp(eta)/(1+exp(eta))
  
  genocase = rep(NA, sample.size)
  for(i in 1:sample.size){
    genocase[i] = rbinom(1, 1, prob[i])
  }
  case.idx = sample(which(genocase==1),ss/2)
  control.idx = sample(which(genocase==0),ss/2)
  return(geno = gt[c(case.idx,control.idx),])
  
}



source("subset_vcf_11Jan2018.R")


result.reg= NULL
result.mix= NULL
result.mix.both = NULL


## Download and use full chromosome genetic map
downloadGeneticMap(4)
readGeneticMap(4)

  vcf_file = filename
  vcf = readVCF( vcf_file, maxNumberOfVariants = 1e6 , min_maf = 1e-6 ,max_maf = 0.01)
  
  print(vcf_file)
  for(kk in 1:10){
    set.seed(kk)
    
    vcf_sub = subset_vcf(kk)
    sample.size=ss*1.5
    
    simData0 = data_sim_H0(vcf_sub,kk)
    Z0.skat = simData0$geno
    nsites = simData0$nvar
    Z1.skat = data_sim(vcf_sub,kk)
    dataset_0 = cbind(apply(Z0.skat,1,sum),c(rep(1,ss/2),rep(0,ss/2)))
    dataset_0 = data.frame(dataset_0)
    dataset_1 = cbind(apply(Z1.skat,1,sum),c(rep(1,ss/2),rep(0,ss/2)))
    dataset_1 = data.frame(dataset_1)
    
    
    
    result.mix = rbind(result.mix,c(kk,BFmixture(dataset_0,nsites),BFmixture(dataset_1,nsites)))
    causal.pool = 1:nsites
    par.theory0 = BFbeta(dataset_0)
    par.theory1 = BFbeta(dataset_1)
    result.reg = rbind(result.reg,c(kk,par.theory0,par.theory1))
    result.mix = data.frame(result.mix)
    result.reg = data.frame(result.reg)
    
    names(result.mix) = c("seed_number","w0_0","K0","BF0","w0_1","K1","BF1")
    names(result.reg) = c("seed_number","K0","BF0","K1","BF1")
  }
  
  
  

cat("results from mixture prior:\n")
result.mix

cat("results from regular prior:\n")
result.reg

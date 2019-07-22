# In this code, first simulate data using Sim1000G package, then run BF analysis with KS test p-value as prior

region = "100k" #gene size
ss = 1000 #sample size
jj = 1 # VCF file index
protect.prop = 0.1 # proportion of protective variants in the gene
Rep  = 1 # simulation replicate number

setwd("/Users/Xiong/Dropbox/BF_papers/Biometrics/Reply/Code/")
filename = list.files("./","*\\.vcf.gz")
if(region=="15k"){
  nvariants = 22
  p = 11
}else if(region=="30k"){
  nvariants = 45
  p=15
}else if(region=="50k"){

  nvariants = 72
  p=18
}else if(region=="100k"){
  nvariants = 145
  p=29
}


library(sim1000G)
## Download and use full chromosome genetic map
downloadGeneticMap(4,dir="/Users/Xiong/Dropbox/BF_papers/Biometrics/Reply/Code/")
#readGeneticMap(4,dir="/home/briollaislab/jxu/BayesFactor/sim1000G/mapfile")

convert.gt = function(xx){
  xx = ifelse(xx==2,1,xx)
  return(xx)
}

# simulate gene under H0
data_sim_H0 = function(vcf_sub,seed.num){
  set.seed(seed.num)
  startSimulation(vcf_sub, totalNumberOfIndividuals = sample.size) #what does this step do, usually how can I choose number of individuals
  SIM$reset() #what does this step do 
  

  id = generateUnrelatedIndividuals(sample.size)
  # print(id)
  gt = retrieveGenotypes(id)
  gt = matrix(unlist(lapply(gt,convert.gt)),sample.size,nvariants)
  freq = apply(gt,2,sum)/(2*nrow(gt))

  if(sum(freq>=0.05)>0){
    gt = gt[,-which(freq>=0.05)]
  }
  set.seed(seed.num)
  genocase = sample(sample.size,ss)
  return(list(geno = gt[genocase,],nvar = ncol(gt)))
}

# simulate gene under H1
data_sim = function(vcf_sub,seed.num){
  set.seed(seed.num)
  startSimulation(vcf_sub, totalNumberOfIndividuals = sample.size) 
  SIM$reset() 
  

  id = generateUnrelatedIndividuals(sample.size)
  # print(id)
  gt = retrieveGenotypes(id)
  gt = matrix(unlist(lapply(gt,convert.gt)),sample.size,nvariants)
  freq = apply(gt,2,sum)/(2*nrow(gt))

  if(sum(freq>=0.05)>0){
    gt = gt[,-which(freq>=0.05)]
  }
  freq = apply(gt,2,sum)/(2*nrow(gt))
  causal = sample(setdiff(1:ncol(gt),which(freq==0)),p)
  
  beta.sign = rep(1,p)
  pos.protect = sample(p,floor(p*protect.prop))
  beta.sign[pos.protect] = (-1)
  
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
  # set.seed(seed.num)
  for(i in 1:sample.size){
    genocase[i] = rbinom(1, 1, prob[i])
  }
  case.idx = sample(which(genocase==1),ss/2)
  control.idx = sample(which(genocase==0),ss/2)
  return(geno = gt[c(case.idx,control.idx),])
  
}
source("Fun_reg_MLE_K_03Jan2018_marg_lik_unify.R")
source("Fun_mix_EM_15Jan2018_marg_lik_unify.R")
source("subset_vcf_11Jan2018.R")

result.reg= NULL
result.mix= NULL

  print(filename[jj])
  vcf_file = paste(filename[jj]) 
  vcf = readVCF( vcf_file, maxNumberOfVariants = 1e6 , min_maf = 1e-6 ,max_maf = 0.01)
  
  print(vcf_file)

  p_res = NULL
  
  for(kk in (Rep*10+1:10)){
    set.seed(kk)
    
    vcf_sub = subset_vcf(kk)
    sample.size=ss*1.5
    
    # under H0
    simData0 = data_sim_H0(vcf_sub,kk)
    Z0.skat = simData0$geno
    nsites = simData0$nvar
    Z1.skat = data_sim(vcf_sub,kk)
    maf.vec = (apply(Z0.skat,2,sum)+1)/(2*ss) # MAF estimate 
    
  
    dataset_1 = cbind(apply(Z1.skat,1,sum),c(rep(1,ss/2),rep(0,ss/2)))
    dataset_1 = data.frame(dataset_1)
    

    ##################### use adjusted gneotype to calculate p-value for each variant
    maf.std = function(geno.vec,maf){
       return(floor(geno.vec * min(0.01/maf,10)))
    } 
    


    Z1.std = NULL
    for(v in 1:ncol(Z1.skat)){
        Z1.std = cbind(Z1.std,maf.std(Z1.skat[,v],maf.vec[v]))
}



    p_locus_function = function(c1,c2){
      if((c1+c2)<5) return(NA)
      else return(2*pnorm(abs(c1-c2)/(c1+c2)^0.5,lower.tail = F))
    }
    
    
    sum.Z1.case = apply(Z1.std[1:(ss/2),],2,sum)
    sum.Z1.control = apply(Z1.std[(ss/2+1):ss,],2,sum)
    
    p_z1 = mapply(p_locus_function,sum.Z1.case,sum.Z1.control)
    p_z1[which(maf.vec<0.001)]=NA

    p_res = c(p_res,p_z1)

    
    causal.pool = 1:nsites
    BFbeta1 = 2*log(BFbeta(dataset_1)[2])
    BFmix1 = 2*log(BFmixture(dataset_1,nsites)[3])

    result.reg = rbind(result.reg,c(jj,kk,nsites,BFbeta1))
    result.mix = rbind(result.mix,c(jj,kk,nsites,BFmix1))
    
   }
print(gc())
p_res = data.frame(p_res)
names(p_res)="Data1"
result.reg = data.frame(result.reg)
names(result.reg) = c("VCF_file","seed_number","sites","BF_1")
result.mix = data.frame(result.mix)
names(result.mix) = c("VCF_file","seed_number","sites","BF_1")


library(dgof)
p.vec=NULL
p_H0 = p_res[,1]
p_H0 = p_H0[is.na(p_H0)==F]
for(rr in 1:10){
  #      for(rr in which(BFresult$TrueData==1)){
  pv.gene = p_res[((rr-1)*nvariants+1):(rr*nvariants),1]
  pv.gene = pv.gene[is.na(pv.gene)==F]
  PV = ifelse(length(pv.gene)<5,1,ks.test(pv.gene, ecdf(p_H0),alternative = "greater")$p.value)
  p.vec = c(p.vec,PV)
}
result.reg$BF_inf = result.reg$BF_1 - 2*log(p.vec)
result.mix$BF_inf = result.mix$BF_1 - 2*log(p.vec)


#print("Loaded R-bf-permute")
#library(LearnBayes)



#' Wrapper for Bayes factor methods
#'
#' @param snps Variants
#'
#' @return BF
run_BF = function(snps, pheno, method, permuteSamples, KK,  verbose = F) {
  
  # Hyper parameters
  # eta.star, k.star, eta1.star, k1.star, eta0.star, k0.star,w0.par,k.par
  Eta.par = c(8.007217e-03, 2.348949e+04, 9.528792e-03, 2.319012e+04, 8.007217e-03, 2.348949e+04, 7.150010e-01, 1.225571e+04)
  Eta.par.reg = c(2.022313e-03, 4.496896e+04, 2.632804e-03, 3.791696e+04, 2.022313e-03, 4.496896e+04)
  # eta.star, k.star, eta1.star, k1.star, eta0.star, k0.star, eta.tilde, k.tilde, eta1.tilde, k1.tilde, eta0.tilde, k0.tilde,k.par
  Both.par=c(8.007217e-03, 2.348949e+04, 9.528792e-03, 2.319012e+04, 8.007217e-03, 2.348949e+04, 7.357986e-01, 
             2.448824e+02, 6.914371e-01, 1.610376e+02, 7.357986e-01, 2.448824e+02, 1.225571e+04)
  # eta.tilde, k.tilde, eta1.tilde, k1.tilde, eta0.tilde, k0.tilde,eta.par,k.par
  W.par = c(7.342742e-01,1.979900e+02,6.917673e-01,2.075935e+02,7.342742e-01,1.979900e+02,8.741540e-03,1.262786e+04)
  
  
  
  # snps : genotypes coded as NA/0/1/2
  # pheno : phenotypes
  
  
  
  # Recode as 0: no alleles 
  #snps = ( snps > 0 ) ^ 1     
  causal.pool = 1:nrow(snps)
  
  variants_per_individual = apply(snps,2, function(x) sum(x>0,na.rm=T))
  
  # Number of non missing sites per individual
  non_missing_sites = apply(snps,2,function(x) sum(!is.na(x)))
  
  
  if(verbose) cat("nvarlength = ", length(variants_per_individual), "pheno length = ", length(pheno), "\n")
  
  input_data = data.frame(variants=variants_per_individual, pheno = pheno )
  
  if(permuteSamples > 0) input_data$pheno = sample(input_data$pheno)

  
  if(1) {   ## computation of Eta.par per gene
    

    # Vector of invididual p.hat
    input_data$p.hat = variants_per_individual / non_missing_sites 
    
    #input_data$cc = input_data$pheno
    #input_data$geno_sum = input_data$variants
    
    
    p0.hat = mean(input_data$p.hat[which(pheno==0)])
    p1.hat = mean(input_data$p.hat[which(pheno==1)])
    Eta.par.reg[1] = Eta.par.reg[5] = p0.hat
    Eta.par.reg[3] = p1.hat
    
    new.geno_nonzero = input_data[which(variants_per_individual > 0),]
    
    case_nonzero.index = which(new.geno_nonzero$pheno == 1)
    control_nonzero.index = which(new.geno_nonzero$pheno == 0)
    
    p0.hat_nonzero = mean(new.geno_nonzero$p.hat[control_nonzero.index])
    p1.hat_nonzero = mean(new.geno_nonzero$p.hat[case_nonzero.index])
    
    # Eta.par[7] = 1-nrow(new.geno_nonzero)/nrow(new.geno) # Aug 2nd 2016 for TP63
    Eta.par[1] = min(p0.hat_nonzero,p1.hat_nonzero)
    Eta.par[5] = p0.hat_nonzero
    Eta.par[3] = p1.hat_nonzero
    Both.par[1] = Both.par[5] = p0.hat_nonzero
    Both.par[3] = p1.hat_nonzero
    w0.prop1 = 1-length(case_nonzero.index)/sum(input_data$pheno)
    w0.prop0 = 1-length(control_nonzero.index)/(nrow(input_data)-sum(input_data$pheno))
    W.par[1] = W.par[5] = w0.prop0
    W.par[3] = w0.prop1 
    
  }
  
  
  dataset = data.frame(variants=input_data$variants , pheno=input_data$pheno  )
  
  #sample.order = order(dataset$pheno)
  #sample.order
  #dataset = dataset[sample.order ,]
  #snps = snps[, sample.order ]
  
  dataset$sites.num = non_missing_sites
  
  

  
  
  if(verbose) cat("#### Dataset:\n")
  if(verbose) print(dataset)
    
  if(verbose) cat("KK=",KK, "\n");
  if(verbose) cat("params = ", Eta.par.reg , "\n" )
  
  
  BayesFactor = NA
  
  
  if(method == "reg_eta_miss") {
    BayesFactor = try( BF_reg_eta_miss(dataset, Eta.par.reg, KK) )
  }
  
  if(method == "mix_eta") {
    BayesFactor = try( BF_mix_eta(dataset, Eta.par,  nrow(snps)  ) )
  }
 
  
  # BF.mix.eta = 1
  
  if( "try-error" %in% class(BayesFactor) ) BayesFactor = NA
  
  
  return(BayesFactor)
  
  
}



library(parallel)


if(0) {
  system.time ( t<- sapply(1:1000, function(x) ComputeBF(snps,group12,1)) )
  system.time ( t<-mclapply(mc.cores=6,1:100, function(x) ComputeBF(snps,group12,0)) )
}

parapply = function(n, f) { mclapply(mc.cores=mc.cores ,n, f) }
#parapply = function(n, f) { lapply(n, f) }

adaptivePermutation = function(f)  {
  
  n = 0
  b = 0
  
  b = b + sum ( unlist( parapply(1:10, f ) ) )
  n = n + 10
  
  if(b > 3) return(b/n)
  
  permToRun = 40
  
  #while( b < 10 && n < 500000) {
  while( b < 20 && n < 3000000) {
    b = b + sum ( unlist( parapply(1:permToRun, f ) ) )
    n = n + permToRun
    permToRun = round(1.3*permToRun)
  }
  
  
  return ( c(b/n,b,n) )
  
}

if(0) {
  f = function(t=1) { rbinom(1,1,0.0001) > 0 }
  mean(sapply(1:2000,f))
  mc.cores = 8
  print(system.time( adaptivePermutation( f ) ))
  parapply = function(n, f) { lapply(n, f) }
  print(system.time( adaptivePermutation( f ) ))
}



#print("Loaded R-bf-permute")
#library(LearnBayes)
# change
#' Return Hyper parameters
#'
#' @param variants_per_individual number of variants per individual (n)
#' @param non_missing_sites numer of non missing sites per individual
#' @param pheno phenotype in 0,1
#' @param method Which method to use
#'
#' @return BF
#' @seealso \code{\link{BF}}
#' @export
compute_hyper_parameters = function(variants_per_individual, non_missing_sites, pheno, method, verbose = F) {
  
  # Initial Hyper parameters
  # eta.star, k.star, eta1.star, k1.star, eta0.star, k0.star,w0.par,k.par
#   Eta.par = c(8.007217e-03, 2.348949e+04, 9.528792e-03, 2.319012e+04, 8.007217e-03, 2.348949e+04, 7.150010e-01, 1.225571e+04)
  # Eta.par.reg = c(2.022313e-03, 4.496896e+04, 2.632804e-03, 3.791696e+04, 2.022313e-03, 4.496896e+04)
  # eta.star, k.star, eta1.star, k1.star, eta0.star, k0.star, eta.tilde, k.tilde, eta1.tilde, k1.tilde, eta0.tilde, k0.tilde,k.par
#   Both.par=c(8.007217e-03, 2.348949e+04, 9.528792e-03, 2.319012e+04, 8.007217e-03, 2.348949e+04, 7.357986e-01, 
#              2.448824e+02, 6.914371e-01, 1.610376e+02, 7.357986e-01, 2.448824e+02, 1.225571e+04)
  # eta.tilde, k.tilde, eta1.tilde, k1.tilde, eta0.tilde, k0.tilde,eta.par,k.par
#   W.par = c(7.342742e-01,1.979900e+02,6.917673e-01,2.075935e+02,7.342742e-01,1.979900e+02,8.741540e-03,1.262786e+04)
  Eta.par.reg = rep(NA,6)
  # New initial hyper parameters are computed based on new hyperparameter function "Fun_all_par_H0_13Nov2016.R" and 100kb simulated data "_06Dec2016" (11Dec2016)
  Eta.par = c(6.948017e-03, 3.734636e+04, 6.938008e-03, 5.222202e+04, 6.948017e-03, 3.734636e+04, 4.387904e-01, 9.480634e+03)
  Both.par = c(6.948017e-03, 3.734636e+04, 6.938008e-03, 5.222202e+04, 6.948017e-03, 3.734636e+04, 3.603100e-01, 5.657617e+01, 3.590419e-01, 7.986061e+01, 3.603100e-01, 5.657617e+01, 9.480634e+03)
  W.par = c(3.603100e-01, 5.657617e+01, 3.590419e-01, 7.986061e+01, 3.603100e-01, 5.657617e+01, 7.087595e-03, 9.480634e+03)
    
    # Vector of invididual p.hat
    p.hat = variants_per_individual / non_missing_sites 
  
    
    p0.hat = mean(p.hat[ which(pheno==0) ])
    p1.hat = mean(p.hat[ which(pheno==1) ])
    
    if(method == "reg_eta_miss" || method == "reg_eta") {
  
      # Eta.par.reg[1] = mean(p.hat)
      # Eta.par.reg[5] = p0.hat
      # Eta.par.reg[3] = p1.hat
      # return (Eta.par.reg);
    
      Eta.par.reg[1] = mean(p.hat)
      Eta.par.reg[5] = p0.hat
      Eta.par.reg[3] = p1.hat
      Eta.par.reg[2] = length(non_missing_sites) *length(pheno)
      Eta.par.reg[4] = length(non_missing_sites) *sum(pheno==1)
      Eta.par.reg[6] = length(non_missing_sites) *sum(pheno==0)
      return (Eta.par.reg);
    }
    
    case_nonzero.index    = which(pheno == 1 &  variants_per_individual > 0 )
    control_nonzero.index = which(pheno == 0 &  variants_per_individual > 0 )
    
    p0.hat_nonzero = mean( p.hat[ control_nonzero.index ] )
    p1.hat_nonzero = mean( p.hat[ case_nonzero.index    ] )
    
    if(method == "mix_eta") {

        # Eta.par[7] = 1-nrow(new.geno_nonzero)/nrow(new.geno) # Aug 2nd 2016 for TP63
        Eta.par[1] = min(p0.hat_nonzero,p1.hat_nonzero)
        Eta.par[5] = p0.hat_nonzero
        
        #Eta.par[1] = Eta.par[5] = p0.hat_nonzero # Old method
        
        
        Eta.par[3] = p1.hat_nonzero
        return ( Eta.par )
    }
    
    
    if(method == "mix_both") {
      
      Both.par[1] = Both.par[5] = p0.hat_nonzero
      Both.par[3] = p1.hat_nonzero
      return ( Both.par )
      
    }
    
    if(method == "mix_w0") {

      w0.prop1 = 1-length(case_nonzero.index)/sum(pheno)
      w0.prop0 = 1-length(control_nonzero.index)/(length(pheno)-sum(pheno))
      
      W.par[1] = W.par[5] = w0.prop0
      W.par[3] = w0.prop1 
      return(W.par)
      
    }
    
    

  
  
  return( list(Eta.par = Eta.par, Eta.par.reg = Eta.par.reg, Both.par= Both.par, W.par = W.par))
  
}


#' Wrapper for Bayes factor methods
#'
#' @param variants_per_individual vector
#' @param non_missing_sites vector
#' @param pheno vector
#'
#' @return BF
run_BF = function(variants_per_individual, non_missing_sites, pheno, method, permuteSamples, KK,  hyper, verbose = F) {

  
  
  
  
  
  if(verbose) cat("nvarlength = ", length(variants_per_individual), "pheno length = ", length(pheno), "\n")
  
  if(0) {   ## computation of Eta.par per gene
      
    # Hyper parameters
    # eta.star, k.star, eta1.star, k1.star, eta0.star, k0.star,w0.par,k.par
    Eta.par = c(8.007217e-03, 2.348949e+04, 9.528792e-03, 2.319012e+04, 8.007217e-03, 2.348949e+04, 7.150010e-01, 1.225571e+04)
    #Eta.par.reg = c(2.022313e-03, 4.496896e+04, 2.632804e-03, 3.791696e+04, 2.022313e-03, 4.496896e+04)
    # eta.star, k.star, eta1.star, k1.star, eta0.star, k0.star, eta.tilde, k.tilde, eta1.tilde, k1.tilde, eta0.tilde, k0.tilde,k.par
    Both.par=c(8.007217e-03, 2.348949e+04, 9.528792e-03, 2.319012e+04, 8.007217e-03, 2.348949e+04, 7.357986e-01, 
               2.448824e+02, 6.914371e-01, 1.610376e+02, 7.357986e-01, 2.448824e+02, 1.225571e+04)
    # eta.tilde, k.tilde, eta1.tilde, k1.tilde, eta0.tilde, k0.tilde,eta.par,k.par
    W.par = c(7.342742e-01,1.979900e+02,6.917673e-01,2.075935e+02,7.342742e-01,1.979900e+02,8.741540e-03,1.262786e+04)

  
    input_data = data.frame(variants=variants_per_individual, pheno = pheno )

    # Vector of invididual p.hat
    input_data$p.hat = variants_per_individual / non_missing_sites
    
    #input_data$cc = input_data$pheno
    #input_data$geno_sum = input_data$variants
    
    
    p0.hat = mean(input_data$p.hat[which(pheno==0)])
    p1.hat = mean(input_data$p.hat[which(pheno==1)])
    Eta.par.reg[1] = mean(input_data$p.hat)
    Eta.par.reg[5] = p0.hat
    Eta.par.reg[3] = p1.hat
    Eta.par.reg[2] = length(non_missing_sites) *length(pheno)
    Eta.par.reg[4] = length(non_missing_sites) *sum(pheno==1)
    Eta.par.reg[6] = length(non_missing_sites) *sum(pheno==0)
    
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
  
  
  dataset = data.frame(variants=variants_per_individual , pheno=pheno  )
  if(permuteSamples > 0) dataset$pheno = sample(dataset$pheno)

  #sample.order = order(dataset$pheno)
  #sample.order
  #dataset = dataset[sample.order ,]
  #snps = snps[, sample.order ]
  
  par = NA
  

  if( class(hyper) == "function") {
    
    par = hyper(variants_per_individual, non_missing_sites, pheno, method)
    
  } else if(is.na(hyper)) {
    
    par = compute_hyper_parameters(variants_per_individual, non_missing_sites, pheno, method )
    
  } else {
    
    par = hyper  
    
  }
  
  
  #last_par  <<- par
  #cat("par=", par)
  
  #print(Eta.par.reg)
  #print(t - Eta.par.reg)
  #sfdsd();
  
  
  
  if(verbose) cat("#### Dataset:\n")
  if(verbose) print(dataset)
  if(verbose) cat("KK=",KK, "\n");
  if(verbose) cat("params = ", par , "\n" )
  
  
  BayesFactor = NA
  
  
  if(method == "reg_eta_miss") {
    dataset$sites.num = non_missing_sites
    Eta.par.reg = par
    BayesFactor = try( BF_reg_eta_miss(dataset, Eta.par.reg, KK) )
  }
  
  if(method == "reg_eta") {
    Eta.par = par
    BayesFactor = try( BF_reg_eta(dataset, Eta.par, KK,  non_missing_sites[1]  ) )
  }
  
  if(method == "mix_eta") {
    Eta.par = par
    BayesFactor = try( BF_mix_eta(dataset, Eta.par,  non_missing_sites[1]  ) )
  }
  
  if(method == "mix_both") {
    Both.par = par
    BayesFactor = try( BF_mix_eta(dataset, Both.par, non_missing_sites[1]  ) )
  }
 
  if(method == "mix_w0") {
    W.par = par
    BayesFactor = try( BF_mix_w0(dataset, W.par,  non_missing_sites[1]  ) )
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





#' Adaptive truncated permuatation methods
#' 
#' @export
adaptivePermutation = function(f, maxPerm = 3e6, parapply = lapply)  {
  
  n = 0
  b = 0
  
  b = b + sum ( unlist( parapply(1:10, f ) ) )
  n = n + 10
  cat("b=",b,n,"\n")
  
  if(b > 3) return ( c(b/n,b,n) )
  
  permToRun = 40
  
  #while( b < 10 && n < 500000) {
  while( b < 20 && n < maxPerm) {
    b = b + sum ( unlist( parapply(1:permToRun, f ) ) )
    n = n + permToRun
    permToRun = round(1.3*permToRun)
    if(permToRun > maxPerm) permToRun = maxPerm
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


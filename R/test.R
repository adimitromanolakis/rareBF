

test_null = function()  {
  
source("R/BF.R")
source("R/BF-package.r")
source("R/R-bf-permute.R")
source("R/BF_reg_eta_miss_15Sep2016_mod.R")


set.seed(101)

Nsamples = 40
Nsites = 50

pheno = ( runif(Nsamples) > 0.5 ) ^ 1

v = round ( rexp(Nsamples * Nsites, rate=0.1) /50 ) 
variants = matrix(v, ncol=Nsamples, nrow=Nsites)
BF(variants,pheno,verbose=F)

# Test for resampling pheno

ns = sample(1:Nsamples)

pheno = pheno[ns]
variants = variants[,ns]
ret = BF(variants,pheno,verbose=T)

# expected  1.002853

cat(ret," ", 1.002853, "\n");


}



test_null2 = function()  {
  
  source("R/BF.R")
  source("R/BF-package.r")
  source("R/R-bf-permute.R")
  source("R/BF_reg_eta_miss_15Sep2016_mod.R")
  
  
  set.seed(101)
  
  Nsamples = 20
  Nsites = 50
  
  pheno = ( runif(Nsamples) > 0.5 ) ^ 1
  
  v = round ( rexp(Nsamples * Nsites, rate=0.2) /50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  BF(variants,pheno,verbose=F)
  
  # Test for resampling pheno
  
  ns = sample(1:Nsamples)
  
  pheno = pheno[ns]
  variants = variants[,ns]
  BF(variants,pheno,verbose=T)
  
  # expected  1.427134
  
  cat(r," ", 1.427134, "\n");
  
  
}


test3 = function()  {
  
  source("R/BF.R")
  source("R/BF-package.r")
  source("R/R-bf-permute.R")
  source("R/BF_reg_eta_miss_15Sep2016_mod.R")
  
  
  set.seed(10)
  
  Nsamples = 20
  Nsites = 50
  
  pheno = ( runif(Nsamples) > 0.5 ) ^ 1
  
  
  v = round ( rexp(Nsamples * Nsites, rate=0.1) / 50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  s = which(pheno == 1)
  for(i in s) variants[,i] = round ( rexp( Nsites, rate=0.15) / 50 )   
  
  
  cat(apply(variants[,pheno==0],1,sum))
  cat(apply(variants[,pheno==1],1,sum))
  
  
  t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "reg_eta_miss") )
  t1  
  r
}



test4 = function()  {
  
  source("R/BF.R")
  source("R/BF-package.r")
  source("R/R-bf-permute.R")
  source("R/BF_reg_eta_miss_15Sep2016_mod.R")
  
  
  set.seed(10)
  
  Nsamples = 150
  Nsites = 1500
  
  pheno = ( runif(Nsamples) > 0.5 ) ^ 1
  
  
  v = round ( rexp(Nsamples * Nsites, rate=0.1) / 50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  
  
  s = which(pheno == 1)
  for(i in s) variants[,i] = round ( rexp( Nsites, rate=0.105) / 50 )   
  
  t1=(apply(variants[,pheno==0],1,sum))
  t2=(apply(variants[,pheno==1],1,sum))
  cat(mean(t1),mean(t2),"\n")
  
  
  t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "reg_eta_miss") )
  t1  
  r
  cat(r," ", 3036207, "\n");
  
  #expected 3036207
}



run_tests = function() {
  test_null()
  test_null2()
  test3()
  test4()
  
}

#run_tests()

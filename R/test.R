library(LearnBayes)


#' Example 1 for Bayes factor methods
#'
#' @return BF
#' @export
example1 = function() {
  
  
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

cat("Result: " , ret," ", 1.002853, "\n");


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
  r = BF(variants,pheno,verbose=T)
  
  # expected  1.427134
  
  cat("Result: " , r," ", 1.427134, "\n");
  
  
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
  cat("Result: " , r," ", 3036207, "\n");
  
  
  vars = apply(variants,2,function(x) sum(x>0))
  sites = apply(variants,2,length)
  sites
  
  t1 = system.time( r <- BFvector(vars,sites,pheno) )
  t1
  r
  
  
  #expected 3036207
}


test5 = function()  {
  
  source("R/BF.R")
  source("R/BF-package.r")
  source("R/R-bf-permute.R")
  source("R/BF_reg_eta_miss_15Sep2016_mod.R")
  
  
  set.seed(10)
  
  Nsamples = 550
  Nsites = 5500
  
  pheno = ( runif(Nsamples) > 0.5 ) ^ 1
  
  
  v = round ( rexp(Nsamples * Nsites, rate=0.1) / 50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  
  
  s = which(pheno == 1)
  for(i in s) variants[,i] = round ( rexp( Nsites, rate=0.103) / 50 )   
  
  t1=(apply(variants[,pheno==0],1,sum))
  t2=(apply(variants[,pheno==1],1,sum))
  cat(mean(t1),mean(t2),"\n")
  cat("Missing rate:"  , sum(is.na(variants)) / (nrow(variants)*ncol(variants)) , "\n")
  
  
  t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "reg_eta_miss") )
  t1  
  r
  cat("Result: " , r," ", 162654.3, "\n");
  
  #expected 162654.3
}



test6 = function()  {
  
  source("R/BF.R")
  source("R/BF-package.r")
  source("R/R-bf-permute.R")
  source("R/BF_reg_eta_miss_15Sep2016_mod.R")
  
  
  
  set.seed(10)
  
  Nsamples = 550
  Nsites = 5500
  
  pheno = ( runif(Nsamples) > 0.5 ) + 0
  
  
  v = round ( rexp(Nsamples * Nsites, rate=0.1) / 50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  
  
  s = which(pheno == 1)
  for(i in s) variants[,i] = round ( rexp( Nsites, rate=0.103) / 50 )   
  
  
  for(i in 1:Nsamples) {
    variants[ which( runif(Nsites) > 0.9 )    ,i] = NA
  }
  
  
  t1=(apply(variants[,pheno==0],1,sum))
  t2=(apply(variants[,pheno==1],1,sum))
  cat(mean(t1,na.rm=T),mean(t2,na.rm=T),"\n")
  
  print(variants[1:10,1:10])
  cat("Missing rate:"  , sum(is.na(variants)) / (nrow(variants)*ncol(variants)) , "\n")
  
  
  t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "reg_eta_miss") )
  t1  
  r
  cat("Result: " , r," ", 106242, "\n");
  
  #expected 162654.3
}



test7 = function()  {
  
  source("R/BF.R")
  source("R/BF-package.r")
  source("R/R-bf-permute.R")
  source("R/BF_reg_eta_miss_15Sep2016_mod.R")
  
  
  
  set.seed(10)
  
  Nsamples = 2550
  Nsites = 9500
  
  pheno = ( runif(Nsamples) > 0.5 ) ^ 1
  
  
  v = round ( rexp(Nsamples * Nsites, rate=0.1) / 50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  
  
  s = which(pheno == 1)
  for(i in s) variants[,i] = round ( rexp( Nsites, rate=0.1) / 50 )   
  
  
  for(i in 1:Nsamples) {
    variants[ which( runif(Nsites) > 0.9 )    ,i] = NA
  }
  
  
  t1=(apply(variants[,pheno==0],1,sum))
  t2=(apply(variants[,pheno==1],1,sum))
  cat(mean(t1,na.rm=T),mean(t2,na.rm=T),"\n")
  
  
  t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "reg_eta_miss") )
  cat("Time=", t1, "\n");  
  r
  cat(r," ", 106242, "\n");
  
  #expected 162654.3
}


test_mix_eta_1 = function()  {
  
  source("R/BF.R")
  source("R/BF-package.r")
  source("R/R-bf-permute.R")
  source("R/BF_reg_eta_miss_15Sep2016_mod.R")
  source("R/BF_mix_eta_25Apr2016.R")
  
  
  set.seed(10)
  
  Nsamples = 40
  Nsites = 50
  
  pheno = ( runif(Nsamples) > 0.5 ) ^ 1
  
  
  v = round ( rexp(Nsamples * Nsites, rate=0.1) / 50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  
  
  s = which(pheno == 1)
  for(i in s) variants[,i] = round ( rexp( Nsites, rate=0.2) / 50 )   
  
  
  if(0) for(i in 1:Nsamples) {
    variants[ which( runif(Nsites) > 0.9 )    ,i] = NA
  }
  
  
  t1=(apply(variants[,pheno==0],1,sum))
  t2=(apply(variants[,pheno==1],1,sum))
  cat(mean(t1,na.rm=T),mean(t2,na.rm=T),"\n")
  
  
  t1 = system.time( r <- BF(variants,pheno,verbose=T, method = "mix_eta") )
  t1  
  r
  cat(r," ", 106242, "\n");
  
  #expected 162654.3
}




test_custom = function(Nsamples,Nsites, missrate, r1,r2, nrepeat=1 )  {
  
  source("R/BF.R")
  source("R/BF-package.r")
  source("R/R-bf-permute.R")
  source("R/BF_reg_eta_miss_15Sep2016_mod.R")
  
  
  
  set.seed(10)
 
  pheno = ( runif(Nsamples) > 0.5 ) + 0
  
  
  v = round ( rexp(Nsamples * Nsites, rate=r1) / 50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  
  
  s = which(pheno == 1)
  for(i in s) variants[,i] = round ( rexp( Nsites, rate=r2) / 50 )   
  
  
  for(i in 1:Nsamples) {
    variants[ which( runif(Nsites) > (1-missrate) )    ,i] = NA
  }
  
  
  t1=(apply(variants[,pheno==0],1,sum))
  t2=(apply(variants[,pheno==1],1,sum))
  
  print(variants[1:10, which( pheno==0)[1:10]])
  print(variants[1:10, which( pheno==1)[1:10]])
  cat("Mean variants: " , mean(t1,na.rm=T),mean(t2,na.rm=T),"\n")
  cat("Total variants: " , sum(t1,na.rm=T),sum(t2,na.rm=T),"\n")
  
  cat("Missing rate:"  , sum(is.na(variants)) / (nrow(variants)*ncol(variants)) , "\n")
  
  
  t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "reg_eta_miss") )
  cat("Time=", t1, "\n");  
  
  r
  cat("Result: " , r," ", 106242, "\n");
  
  
  if(nrepeat > 1) {
    f = function() {
      r=0
      for(i in 1:nrepeat) r = r + BF(variants,pheno,verbose=F, method = "reg_eta_miss") 
      return(r);
    }
    t1 = system.time( r <- f() );
    cat("Time x",nrepeat, " = ", t1, "\n");  
    
  }
  
  #expected 162654.3
}



run_tests = function() {
  test_null()
  test_null2()
  test3()
  test4()
  
}

# run_tests()

# rareBF


R package for Bayesian Models for Rare Variant Association Analysis.

#### Author: Laurent Briollais, Jingxiong Xu 
#### Maintainer: Apostolos Dimitromanolakis

## Installation

Users can install the most recent version by running (in R):

```R
install.packages("https://github.com/adimitromanolakis/rareBF/releases/download/v1.03/BF_1.03.tar.gz",repos=NULL, verbose=T)
```

## Basic Usage

The following example demonstrates simulating data from the null model:


```R
  library(BF)
  
  set.seed(101)
  Nsamples = 40
  Nsites = 50
  
  pheno = ( runif(Nsamples) > 0.5 ) ^ 1
  
  v = round ( rexp(Nsamples * Nsites, rate=0.1) /50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
  
  
  bf = BF(variants, pheno, method="reg_eta_miss", verbose=F)
 
  # expected  1.002853
  
  cat(bf," ", 1.002853, "\n");

```


## Testing of mixed eta prior


```R

  library(BF)
  set.seed(101)
  
  Nsamples = 40
  Nsites = 50
  
  pheno = ( runif(Nsamples) > 0.5 ) ^ 1
  
  
  v = round ( rexp(Nsamples * Nsites, rate=0.1) / 50 ) 
  variants = matrix(v, ncol=Nsamples, nrow=Nsites)
    
  s = which(pheno == 1)
  for(i in s) variants[,i] = round ( rexp( Nsites, rate=0.2) / 50 )   
  
  t1=(apply(variants[,pheno==0],1,sum))
  t2=(apply(variants[,pheno==1],1,sum))
  cat(mean(t1,na.rm=T),mean(t2,na.rm=T),"\n")
  
  
  t1 = system.time( r <- BF(variants,pheno,verbose=T, method = "mix_eta") )
  t1  
  r
  cat(r," ", 54637, "\n");

```



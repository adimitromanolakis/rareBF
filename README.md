# rareBF


R package for Bayesian Models for Rare Variant Association Analysis.

#### Author: Laurent Briollais, Jingxiong Xu 
#### Maintainer: Apostolos Dimitromanolakis

## Installation

Users can install the most recent version by running (in R):

```R
install.packages("https://github.com/adimitromanolakis/rareBF/releases/download/v1.04/rareBF_1.04.tar.gz",repos=NULL, verbose=T)
```

## Basic Usage

The following example demonstrates simulating data from the null model:


```R
  library(rareBF)
  
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

  library(rareBF)
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


## Use with modified hyper parameters

The following example illustrates how you can modify the default hyper parameters when running the BF function:

```R
library(rareBF)

set.seed(101)
Nsamples = 40
Nsites = 50

pheno = ( runif(Nsamples) > 0.5 ) ^ 1

v = round ( rexp(Nsamples * Nsites, rate=0.1) /50 ) 
variants = matrix(v, ncol=Nsamples, nrow=Nsites)

hyper_parameters = c(0.08210526 ,44968.96, 0.08190476 ,37916.96, 0.08210526, 44968.96 )

bf = BF(variants, pheno, method="reg_eta_miss", verbose=F, hyper=hyper_parameters)

# expected  1.002853

cat(bf," ", 1.002853, "\n");
```

## Comparison of methods


```R
library(rareBF)


set.seed(10)

Nsamples = 150
Nsites = 100

pheno = ( runif(Nsamples) > 0.5 ) ^ 1


v = round ( rexp(Nsamples * Nsites, rate=0.01) / 50 ) 
variants = matrix(v, ncol=Nsamples, nrow=Nsites)



s = which(pheno == 1)
for(i in s) variants[,i] = round ( rexp( Nsites, rate=0.012) / 50 )   

t1=(apply(variants[,pheno==0],1,sum))
t2=(apply(variants[,pheno==1],1,sum))
cat(mean(t1),mean(t2),"\n")


t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "reg_eta_miss") )

cat("------ reg_eta_miss  BF=", r , "   time spent=", t1[1])

t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "mix_w0") )

cat("------ mix_w0  BF=", r , "   time spent=", t1[1])

t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "mix_eta") )

cat("------ mix_eta  BF=", r , "   time spent=", t1[1])

t1 = system.time( r <- BF(variants,pheno,verbose=F, method = "mix_both") )

cat("------ mix_both  BF=", r , "   time spent=", t1[1])

```



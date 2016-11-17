# rareBF


R package for Bayesian Models for Rare Variant Association Analysis.

#### Author: Laurent Briollais, Jingxiong Xu 
#### Maintainer: Apostolos Dimitromanolakis

## Installation

Users can install the most recent version by running (in R):

```R
 install.packages("https://github.com/adimitromanolakis/rareBF/releases/download/1.01/BF_1.01.tar.gz")

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
  bf = BF(variants,pheno,verbose=F)
 
  # expected  1.002853
  
  cat(bf," ", 1.002853, "\n");

```

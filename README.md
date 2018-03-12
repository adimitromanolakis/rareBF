# rareBF


R package for Bayesian Models for Rare Variant Association Analysis.

#### Author: Laurent Briollais, Jingxiong Xu 
#### Maintainer: Apostolos Dimitromanolakis

## Installation

Users can install the most recent version by running (in R):

```R
install.packages("LearnBayes")
devtools::install_github("adimitromanolakis/rareBF")
```

If that fails, you can try downloading the package source from: [https://github.com/adimitromanolakis/rareBF/files/1803537/rareBF_2.0.tar.gz](https://github.com/adimitromanolakis/rareBF/files/1803537/rareBF_2.0.tar.gz) and installing manually:

```R
install.packages("rareBF_2.0.tar.gz")
```



On MACOS or Windows, you can download the [zip file](https://github.com/adimitromanolakis/rareBF/archive/v2.0.zip) and install it.



## Basic Usage



The following example demonstrates simulating data from the null and an alternative model:


```R


library(rareBF)

set.seed(101)

nsites = 500

ncase = 200
ncontrol = 200

maf = 0.005
rexp_param = 2


model1 = function(i) {
    
    p1 = rbinom( ncontrol  , nsites , maf   )
    p2 = rbinom( ncase     , nsites , maf )  + H1 * round(rexp(ncase,rexp_param))
    
    disease_status = c (  rep(0,ncontrol),rep(1,ncase)  )
    data = data.frame(c(p1,p2), disease_status)
    
    bf = BFmixture(data, nsites)[3]
    bf
}

nrep = 20

H1 = 0
bf1 = ( replicate(nrep,model1() ) )

H1 = 1
bf2 = ( replicate(nrep,model1() ) )

col = c( rep("blue",nrep) , rep("red", nrep))


plot(log(c(bf1,bf2)),col = col ,pch=20,cex=1.5)



```


A more involved example can be found at the inst/example directory or at github:

[https://github.com/adimitromanolakis/rareBF/tree/master/inst/examples]



library(rareBF)

set.seed(101)

nsites = 500

ncase = 200
ncontrol = 200

maf = 0.005

model1 = function(i) {
    
    p1 = rbinom( ncontrol  , nsites , maf   )
    p2 = rbinom( ncase     , nsites , maf )  + round(rexp(ncase,rexp_param))
    
    disease_status = c (  rep(0,ncontrol),rep(1,ncase)  )
    data = data.frame(c(p1,p2), disease_status)
    
    bf = BFmixture(data, nsites)[3]
    bf
}

rexp_param = 100000

nrep = 20

bf1 = ( replicate(nrep,model1() ) )
rexp_param = 2

bf2 = ( replicate(nrep,model1() ) )

col = c( rep("blue",nrep) , rep("red", nrep))


plot(log(c(bf1,bf2)),col = col ,pch=20,cex=1.5)



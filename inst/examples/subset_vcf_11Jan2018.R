subset_vcf = function(seed_number)  {
  # 
  # # plot(vcf$maf)
  # quantile(vcf$maf)
  
  
  subsetVCF = function(vcf, var_index = NA) {
    
    
    vcf2 = as.list(vcf)
    
    s = var_index
    
    vcf2$maf = vcf2$maf [s]
    vcf2$varid = vcf2$varid [s]
    
    vcf2$gt1 = vcf2$gt1 [s,]
    vcf2$gt2 = vcf2$gt2 [s,]
    
    vcf2$vcf = vcf2$vcf [s,]
    
    as.environment(vcf2)
   }
  # str(as.list(vcf))
  # # find number of variants in the vcf file
  num_variants = nrow(vcf$gt1)
  # num_individuals = ncol(vcf$gt1)
  # 
  # 
  # select first 300 variants
  set.seed(seed_number)
  var_select = sample(1:num_variants,nvariants)
  vcf2 = subsetVCF(vcf, var_select)
  
  return(vcf2)
}



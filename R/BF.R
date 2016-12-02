


#' Run Bayes factor analysis
#'
#' @param nvariants matrix Data frame of sites ( rows ) x individuals ( columns )
#' @param pheno vector Phenotypes (0/1)
#' @param method One of the Bayes Factor methods to use (default reg_eta_miss)
#' @param param Fine tuning of parameter set
#' @param KK Default is 500
#' @param hyper Specify hyper parameters or function returning hyper parameters
#' @param verbose Print additional debugging information
#'
#' @return Bayes Factor (numeric)
#' 
#' @details
#'  nvariants contains a matrix of variants per site. Missing data are coded as NA. 
#'  The following methods are available:
##' \itemize{
##'  \item{ \code{reg_eta_miss} : }{Regular prior, able to handle missing data}
##'  \item{ \code{mix_eta} : }{Mixed prior}
##'  \item{ \code{mix_both} : }{Mixed prior both}
##'  \item{ \code{mix_w0} : }{Mixed prior w0}

##' }
#'  More information is available at the following link: \url{https://adimitromanolakis.github.io/rareBF/}
#' 
#' 
#' 
#' @examples
#' #################################################################################
#' ## Example simulating data from null model
#' #################################################################################
#'
#' Nsamples = 20
#' Nsites = 50
#' pheno = ( runif(Nsamples) > 0.2 ) ^ 1
#'
#' v = runif(Nsamples * Nsites) > 0.3
#' variants = matrix(v, ncol=Nsamples, nrow=Nsites)
#' BF(variants,pheno,verbose=TRUE)
#'
#' @seealso \code{\link{BF}}
#' 
#' @export
BF = function(variants, pheno, method = "reg_eta_miss",  KK = 500, hyper = NA, verbose = FALSE) {
  
  if(! ( class(pheno) == "numeric" || class(pheno) == "integer" ) ) {
    stop("pheno must be a numeric vector")
  }
  
  if(! identical( sort(unique(pheno)) , c(0,1) ) ) {
    stop("phenotype vector must contain only 0 or 1")
  }
  
  if(! ( class(variants) == "data.frame" || class(variants) == "matrix" ) ) {
    stop("variants must be a data.frame or matrix")
  } 
  
  
  if(1) if(ncol(variants) != length(pheno)) {
    stop("pheno and variants specify different number of individuals")  
  }

  available_methods = c("reg_eta_miss", "mix_eta", "mix_both", "mix_w0")
  
  if(!(method %in% available_methods)) {
    stop("method argument not recognized")  
  }
  
  
  
  # snps : genotypes coded as NA/0/1/2
  # pheno : phenotypes
  
  
  
  # Recode as 0: no alleles 
  # snps = ( snps > 0 ) ^ 1     
  # causal.pool = 1:nrow(snps)
  
  variants_per_individual = apply(variants,2, function(x) sum(x>0,na.rm=T))
  
  # Number of non missing sites per individual
  # For all methods except reg_eta_miss , use just number of sites
  
  if(method == "reg_eta_miss")
  { non_missing_sites = apply(variants,2,function(x) sum(!is.na(x))) }
  else
  { sites = nrow(variants); non_missing_sites = rep(sites, length(variants_per_individual) ); }
  
  #non_missing_sites = apply(snps,2,function(x) sum(!is.na(x)))
  #cat("Site = ", non_missing_sites,"\n")
  #{ sites = nrow(snps); non_missing_sites = rep(sites, length(variants_per_individual) ); }
  #cat("Site = ", non_missing_sites,"\n")
  
  
  returnvalue = NA
  
  returnvalue = run_BF(variants_per_individual, non_missing_sites, pheno, method, F, KK, hyper, verbose)
  
   
  return(returnvalue)
}







#' BF method for vector data
#'
#' @param variants vector (required) Number of sites per individual
#' @param nsites vector  (required) Number of (non-missing) sites per individual
#' @param pheno vector  (required) Phenotypes (0/1)
#' @param method One of the Bayes Factor methods to use (default reg_eta_miss)
#' @param param Fine tuning of parameter set
#' @param KK Default is 500
#' @param hyper Specify hyper parameters or function returning hyper parameters
#' @param verbose Print additional debugging information
#'
#' @return Bayes Factor (numeric)
#' 
#' @details
#'  nvariants contains a matrix of variants per site. Missing data are coded as NA. 
#'  The following methods are available:
##' \itemize{
##'  \item{ \code{reg_eta_miss} : }{Regular prior, able to handle missing data}
##'  \item{ \code{mix_eta} : }{Mixed prior}
##'  \item{ \code{mix_both} : }{Mixed prior both}
##'  \item{ \code{mix_w0} : }{Mixed prior w0}

##' }
#'  More information is available at the following link: \url{https://adimitromanolakis.github.io/rareBF/}
#' 
#' 
#' 
#' @examples
#' #################################################################################
#' ## Example simulating data from null model
#' #################################################################################
#'
#' Nsamples = 20
#' Nsites = 50
#' pheno = ( runif(Nsamples) > 0.2 ) ^ 1
#'
#' v = runif(Nsamples * Nsites) > 0.3
#' variants = matrix(v, ncol=Nsamples, nrow=Nsites)
#' BF(variants,pheno,verbose=TRUE)
#'
#' @seealso \code{\link{BF}}
#' 
#' @export
BFvector = function(variants, nsites, pheno, method = "reg_eta_miss",  KK = 500, hyper = NA, permutations = NA, verbose = FALSE) {
  
  
  if(! ( class(pheno) == "numeric" || class(pheno) == "integer" ) ) {
    stop("pheno must be a numeric vector")
  }
  
  if(! identical( sort(unique(pheno)) , c(0,1) ) ) {
    stop("phenotype vector must contain only 0 or 1")
  }
  
  if(1) if(length(variants) != length(pheno) || length(nsites) != length(pheno)  )  {
    stop("pheno and variants specify different number of individuals")  
  }
  
  available_methods = c("reg_eta_miss", "mix_eta", "mix_both", "mix_w0")
  
  if(!(method %in% available_methods)) {
    stop("method argument not recognized")  
  }
  
  
  
  returnvalue = NA
  
  if(!is.na(permutations)) {
    
    returnvalue = run_BF(variants, nsites, pheno, method, T, KK, hyper, verbose)
    return (returnvalue)
    
  }
  
  returnvalue = run_BF(variants, nsites, pheno, method, F, KK, hyper, verbose)
  
  

  
  
  return(returnvalue)
}









#' Run permutations for BF method for vector data
#'
#' @param variants vector (required) Number of sites per individual
#' @param nsites vector  (required) Number of (non-missing) sites per individual
#' @param pheno vector  (required) Phenotypes (0/1)
#' @param method One of the Bayes Factor methods to use (default reg_eta_miss)
#' @param param Fine tuning of parameter set
#' @param KK Default is 500
#' @param hyper Specify hyper parameters or function returning hyper parameters
#' @param verbose Print additional debugging information
#'
#' @return Bayes Factor (numeric)
#' 
#' @details
#'  More information is available at the following link: \url{https://adimitromanolakis.github.io/rareBF/}
#' @seealso \code{\link{BFvector}}
#' 
#' @export
BFvectorPermutations = function(variants, nsites, pheno, method = "reg_eta_miss",  KK = 500, hyper = NA, permutations = NA, verbose = FALSE) {
  
  
  if(! ( class(pheno) == "numeric" || class(pheno) == "integer" ) ) {
    stop("pheno must be a numeric vector")
  }
  
  if(! identical( sort(unique(pheno)) , c(0,1) ) ) {
    stop("phenotype vector must contain only 0 or 1")
  }
  
  if(1) if(length(variants) != length(pheno) || length(nsites) != length(pheno)  )  {
    stop("pheno and variants specify different number of individuals")  
  }
  
  available_methods = c("reg_eta_miss", "mix_eta", "mix_both", "mix_w0")
  
  if(!(method %in% available_methods)) {
    stop("method argument not recognized")  
  }
  
  if(is.na(permutations)) {
    stop("permutations number must be specified")
  }
  

  
  originalBF = run_BF(variants, nsites, pheno, method, F, KK, hyper, verbose)
  
  cat("originalBF=", originalBF, "\n");
  
  f = function(x) { 
    cat(originalBF); 
    r = run_BF(variants, nsites, pheno, method, T, KK, hyper, F) > originalBF 
    r
  }
  
  permResult = adaptivePermutation(f, permutations, lapply)

  ret = list(originalBF=originalBF, p.value = permResult[1], n.success = permResult[2], n.total = permResult[3])
  
  return(ret)
}




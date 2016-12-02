


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
BF = function(variants, pheno, method = "reg_eta_miss", param, KK = 500, hyper = NA, verbose = FALSE) {
  
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
  
  returnvalue = NA
  
  #variants = variants^1
  returnvalue = run_BF(variants, pheno, method, F, KK, hyper, verbose)
  
  
  return(returnvalue)
}







#' BF method for vector data
#'
#' @param variants vector Number of sites per individual
#' @param nsites vector Number of (non-missing) sites per individual
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
BFvector = function(variants, nsites, pheno, method = "reg_eta_miss", param, KK = 500, hyper = NA, verbose = FALSE) {
  
}






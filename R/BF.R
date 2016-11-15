


#' Run Bayes factor analysis
#'
#' @param nvariants matrix Data frame of sites (columns) x individuals ( rows )
#' @param pheno vector Phenotypes (0/1)
#' @param method One of the Bayes Factor methods to use (default reg_eta_miss)
#' @param param Fine tuning of parameter set
#' @param KK Default is 500
#'
#' @return Bayes Factor (numeric)
#' 
#' @details
#'  Documentation for function BF.
#'  The following methods are available:
##' \itemize{
##'  \item{ \code{reg_eta_miss} : }{--}
##'  \item{ \code{reg_eta} : }{--}
##' }
#'  \url{http://www.r-project.org}
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
BF = function(variants, pheno, method = "reg_eta_miss", param, KK = 500, verbose = FALSE) {
  
  if(! ( class(pheno) == "numeric" || class(pheno) == "integer" ) ) {
    stop("pheno must be a numeric vector")
  } 
  
  if(! ( class(variants) == "data.frame" || class(variants) == "matrix" ) ) {
    stop("variants must be a data.frame or matrix")
  } 
  
  
  if(1) if(ncol(variants) != length(pheno)) {
    stop("pheno and variants specify different number of individuals")  
  }

  returnvalue = NA
  returnvalue = run_BF(variants, pheno, method, F, KK, verbose)
  
  available_methods = c("reg_eta_miss")
  
  if(!(method %in% available_methods)) {
    stop("method argument not recognized")  
  }  
  
  
  return(returnvalue)
}







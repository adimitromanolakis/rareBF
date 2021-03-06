% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fun_reg_MLE_K_03Jan2018_marg_lik.R
\name{BFbeta}
\alias{BFbeta}
\title{calculate BF with beta prior}
\usage{
BFbeta(obs.data)
}
\arguments{
\item{obs.data}{N*2 data matrix, N is the sample size, 1st column is sum of rare variants for each individual, 2nd column is the disease status}
}
\value{
Vector of 2 elements. 1st: precision parameter of beta distribution; 2nd: BF
}
\description{
calculate BF with beta prior
}

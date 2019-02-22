#' A spatial simulation function
#'
#' This function allows you to simulate a spatial process by specifying the known parameters.
#' @name sprocess$fitmodel
#' @usage sprocess$fitmodel
#' @param x1 a vector of longitudes
#' @param x2 a vector of latitudes
#' @param sigmasq.true true sigma squared value
#' @param tausq.true true tau squared value
#' @param phi.true true phi value
#' @param phi.upper upper bound of a discrete uniform prior on phi
#' @param tausq.upper upper bound of a discrete uniform prior on tau squared
#' @keywords spatial
#' @export
#' @examples fitmodel(x1,x2,10,0.15,0.5,0.7,1)
#' sprocess_function()




library(readr)
library(devtools)
library(roxygen2)
library(spBayes)
library(geoR)
library(dplyr)
library(ggplot2)
library(mnormt)
library(rjags)

fitmodel<-function(x1,x2, sigmasq.true, tausq.true, phi.true, phi.upper,tausq.upper){

  x1<-x1
  x2<-x2
  ex.data <- grf(length(x1), cov.pars=c(sigmasq.true, phi.true), cov.model="exponential", nugget = tausq.true)
  # defining the grid of prediction locations:
  ex.grid <- as.matrix(expand.grid(seq(0,1,l=15), seq(0,1,l=15)))
  #
  # computing posterior and predictive distributions
  # (warning: the next command can be time demanding)
  ex.bayes <- krige.bayes(ex.data, loc=ex.grid,
                          model = model.control(cov.m="exponential"),
                          prior = prior.control(beta.prior = 'flat',
                                                sigmasq.prior = 'reciprocal',
                                                phi.discrete=seq(0, phi.upper, l=25),
                                                phi.prior="uniform",
                                                tausq.rel.discrete = seq(0, tausq.upper, l=25),
                                                tausq.rel.prior = 'uniform'))
  par(mfrow=c(2,2))
  hist(ex.bayes)
  par(mfrow=c(1,1))


}

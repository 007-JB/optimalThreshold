#################################################################################################
###                                  ClassUndefined.R                                         ###
#################################################################################################
#' An S4 class to represent an 'undefined' distribution.
#'
#' This class allows to fit an undefined distribution on the marker values \code{x}.
#' @slot x a vector of marker values.
#' @slot n Length of x vector (including NA values). Numeric argument.
#' @slot mcmc This slot allows the main function to k,now whether an MCMC algorithm must be performed to sample the distribution parameters from their posterior distribution.
#' @details This class is automatically created when the user applies the \code{fit} function with the argument \code{distr="undefined"}. You never have to create manually this class, it is created internally.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit an undefined distribution.
#' @name undefined-class
#' @aliases undefined
#' @exportClass undefined
setClass("undefined", slots = c(x = "numeric", n = "numeric", mcmc = "logical"))
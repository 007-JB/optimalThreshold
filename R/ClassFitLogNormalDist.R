#################################################################################################
###                                  ClassFitLogNormalDist.R                                  ###
#################################################################################################
#' An S4 class to fit a log-normal distribution on a vector of marker values.
#'
#' This class allows to fit a log-normal distribution on the marker values \code{x}.
#' @slot x This slot takes in argument the marker values. Numeric argument.
#' @slot n Length of x vector (including NA values). Numeric argument.
#' @slot mcmc This slot allows the main function to k,now whether an MCMC algorithm must be performed to sample the distribution parameters from their posterior distribution.
#' @details This class is automatically created when the user applies the \code{fit} function with the argument \code{distr="lnorm"}. You never have to create manually this class, it is created internally.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a log-normal distribution.
#' @name fitLogNormalDist-class
#' @aliases fitLogNormalDist
#' @exportClass fitLogNormalDist
setClass("fitLogNormalDist", slots = c(x = "numeric", n = "numeric", mcmc = "logical"))
setValidity("fitLogNormalDist",
	function(object) {
		if (length(object@x) == 0) return("'x' slot of 'fitLogNormalDist' object cannot be ignored.")
		if (length(object@n) == 0) return("'n' slot of 'fitLogNormalDist' object cannot be ignored.")
		if (length(object@mcmc) == 0) return("'mcmc' slot of 'fitLogNormalDist' object cannot be ignored.")
		if (any(is.na(object@x))) {
			return("NA values in 'x' slot of 'fitLogNormalDist' object.")
		}
		if (any(object@x <= 0)) return("Negative marker values are not supported by the log-normal distribution.")
		else return(TRUE)
	}
)
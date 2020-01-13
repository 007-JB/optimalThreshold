#################################################################################################
###                                  ClassFitLogisticDist.R                                   ###
#################################################################################################
#' An S4 class to fit a logistic distribution on a vector of marker values.
#'
#' This class allows to fit a logistic distribution on the marker values \code{x}.
#' @slot x This slot takes in argument the marker values. Numeric argument.
#' @slot n Length of x vector (including NA values). Numeric argument.
#' @slot ini This slot is a list of initial values passed to the MCMC algorithm. List argument.
#' @slot thin This slot is a strictly positive integer value that specifies the 'thin' applied to the MCMC algorithm.
#' @slot burnin This slot is a positive integer value that specifies the length of the burnin period in the MCMC algorithm.
#' @slot model This slot is a character string that specifies the model passed to the JAGS software to perform the MCMC algorithm.
#' @slot mcmc This slot allows the main function to k,now whether an MCMC algorithm must be performed to sample the distribution parameters from their posterior distribution.
#' @details This class is automatically created when the user applies the \code{fit} function with the argument \code{distr="logis"}. You never have to create manually this class, it is created internally.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a logistic distribution.
#' @name fitLogisticDist-class
#' @aliases fitLogisticDist
#' @exportClass fitLogisticDist
setClass("fitLogisticDist", slots = c(x = "numeric", n = "numeric", ini = "list", thin = "integer", burnin = "integer", model = "character", mcmc = "logical"))
setValidity("fitLogisticDist",
	function(object) {
		if (length(object@x) == 0) return("'x' slot of 'fitLogisticDist' object cannot be ignored.")
		if (length(object@n) == 0) return("'n' slot of 'fitLogisticDist' object cannot be ignored.")
		if (length(object@ini) == 0) return("'ini' slot of 'fitLogisticDist' object cannot be ignored.")
		if (length(object@thin) == 0) return("'thin' slot of 'fitLogisticDist' object cannot be ignored.")
		if (length(object@burnin) == 0) return("'burnin' slot of 'fitLogisticDist' object cannot be ignored.")
		if (length(object@model) == 0) return("'model' slot of 'fitLogisticDist' object cannot be ignored.")
		if (length(object@mcmc) == 0) return("'mcmc' slot of 'fitLogisticDist' object cannot be ignored.")
		if (any(is.na(object@x))) return("NA values in 'x' slot of 'fitLogisticDist' object.")
		if (object@thin <= 0) return("'thin' slot must be positive.")
		if (object@burnin < 0) return("'burnin' slot must be positive, or equal to zero.")
		if (object@model == "") return("'model' slot must be a model that can be used in JAGS.")
		testParaNames <- function(ini) {
			namesParaIni <- names(ini)
			if (length(namesParaIni) != 2) return("Two parameters must be specified: 'location' and 'scale'.")
			if (any(!(namesParaIni %in% c("location", "scale")))) return("Parameters names must be: 'location' and 'scale'.")
			else return(TRUE)
		}
		testParaNamesResults <- sapply(object@ini, testParaNames)
		if (all(is.logical(testParaNamesResults))) {
			if(any(sapply(object@ini, function(x) x$scale <= 0))) return("'scale' parameter must be positive.")
			return(TRUE)
		}
		else return(testParaNamesResults[testParaNamesResults != "TRUE"])
	}
)
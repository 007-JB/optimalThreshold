#################################################################################################
###                                  ClassFitGammaDist.R                                      ###
################################################################################################# 
#' An S4 class to fit a gamma distribution on a vector of marker values.
#'
#' This class allows to fit a gamma distribution on the marker values \code{x}.
#' @slot x This slot takes in argument the marker values. Numeric argument.
#' @slot n Length of x vector (including NA values). Numeric argument.
#' @slot ini This slot is a list of initial values passed to the MCMC algorithm. List argument.
#' @slot thin This slot is a strictly positive integer value that specifies the 'thin' applied to the MCMC algorithm.
#' @slot burnin This slot is a positive integer value that specifies the length of the burnin period in the MCMC algorithm.
#' @slot mcmc This slot allows the main function to k,now whether an MCMC algorithm must be performed to sample the distribution parameters from their posterior distribution.
#' @details This class is automatically created when the user applies the \code{fit} function with the argument \code{distr="gamma"}. You never have to create manually this class, it is created internally.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a gamma distribution.
#' @name fitGammaDist-class
#' @aliases fitGammaDist
#' @exportClass fitGammaDist
setClass("fitGammaDist", slots = c(x = "numeric", n = "numeric", ini = "list", thin = "integer", burnin = "integer", mcmc = "logical"))
setValidity("fitGammaDist",
	function(object) {
		if (length(object@x) == 0) return("'x' slot of 'fitGammaDist' object cannot be ignored.")
		if (length(object@n) == 0) return("'n' slot of 'fitGammaDist' object cannot be ignored.")
		if (length(object@ini) == 0) return("'ini' slot of 'fitGammaDist' object cannot be ignored.")
		if (length(object@thin) == 0) return("'thin' slot of 'fitGammaDist' object cannot be ignored.")
		if (length(object@burnin) == 0) return("'burnin' slot of 'fitGammaDist' object cannot be ignored.")
		if (length(object@mcmc) == 0) return("'mcmc' slot of 'fitGammaDist' object cannot be ignored.")
		if (any(is.na(object@x))) return("NA values in 'x' slot of 'fitGammaDist' object.")
		submittedParaNamesIni <- unique(sapply(object@ini, names))
		if (length(submittedParaNamesIni) == 0) return("No parameter name submitted. Names must be 'shape' for all elements in 'ini'.")
		if (length(submittedParaNamesIni) > 1) return("Multiple parameter names submitted. Names must be 'shape' for all elements in 'ini'.")
		if (submittedParaNamesIni != "shape") return("Parameter names must be 'shape' for all elements in 'ini'.")
		if (any(sapply(object@ini, function(x) x$shape <= 0))) return("'shape' parameters in 'ini' must be positive.")
		if (object@thin <= 0) return("'thin' slot must be positive.")
		if (object@burnin < 0) return("'burnin' slot must be positive, or equal to zero.")
		if (any(object@x < 0)) return("Negative marker values are not supported by the gamma distribution.")
		return(TRUE)
	}
)
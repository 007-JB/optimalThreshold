#################################################################################################
###                                  ClassFitStudentDist.R                                    ###
#################################################################################################
#' An S4 class to fit a Student distribution on a vector of marker values.
#'
#' This class allows to fit a t distribution on the marker values \code{x}.
#' @slot x This slot takes in argument the marker values. Numeric argument.
#' @slot n Length of x vector (including NA values). Numeric argument.
#' @slot ini This slot is a list of initial values passed to the MCMC algorithm. List argument.
#' @slot thin This slot is a strictly positive integer value that specifies the 'thin' applied to the MCMC algorithm.
#' @slot burnin This slot is a positive integer value that specifies the length of the burnin period in the MCMC algorithm.
#' @slot model This slot is a character string that specifies the model passed to the JAGS software to perform the MCMC algorithm.
#' @slot mcmc This slot allows the main function to know whether an MCMC algorithm must be performed to sample the distribution parameters from their posterior distribution.
#' @details This class is automatically created when the user applies the \code{fit} function with the argument \code{distr="t"}. You never have to create manually this class, it is created internally.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a t distribution.
#' @name fitStudentDist-class
#' @aliases fitStudentDist
#' @exportClass fitStudentDist
setClass("fitStudentDist", slots = c(x = "numeric", n = "numeric", ini = "list", thin = "integer", burnin = "integer", model = "character", mcmc = "logical"))
setValidity("fitStudentDist",
	function(object) {
		if (length(object@x) == 0) return("'x' slot of 'fitStudentDist' object cannot be ignored.")
		if (length(object@n) == 0) return("'n' slot of 'fitStudentDist' object cannot be ignored.")
		if (length(object@ini) == 0) return("'ini' slot of 'fitStudentDist' object cannot be ignored.")
		if (length(object@thin) == 0) return("'thin' slot of 'fitStudentDist' object cannot be ignored.")
		if (length(object@burnin) == 0) return("'burnin' slot of 'fitStudentDist' object cannot be ignored.")
		if (length(object@model) == 0) return("'model' slot of 'fitStudentDist' object cannot be ignored.")
		if (length(object@mcmc) == 0) return("'mcmc' slot of 'fitStudentDist' object cannot be ignored.")
		if (any(is.na(object@x))) return("NA values in 'x' slot of 'fitStudentDist' object.")
		if (object@thin <= 0) return("'thin' slot must be positive.")
		if (object@burnin < 0) return("'burnin' slot must be positive, or equal to zero.")
		if (object@model == "") return("'model' slot must be a model that can be used in JAGS.")
		testParaNames <- function(ini) {
			namesParaIni <- names(ini)
			if (length(namesParaIni) != 3) return("Three parameters must be specified: 'mu', 'sd', and 'df'.")
			if (any(!(namesParaIni %in% c("mu", "sd", "df")))) return("Parameters names must be: 'mu', 'sd', and 'df'.")
			else return(TRUE)
		}
		testParaNamesResults <- sapply(object@ini, testParaNames)
		if (all(is.logical(testParaNamesResults))) {
			if (any(sapply(object@ini, function(x) x$sd <= 0))) return("'sd' parameter must be positive.")
			if (any(sapply(object@ini, function(x) x$df <= 0))) return("'df' parameter must be positive.")
			return(TRUE)
		}
		else return(testParaNamesResults[testParaNamesResults != "TRUE"])
	}
)
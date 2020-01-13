#################################################################################################
###                                  ClassFitUserDefinedDist.R                                ###
#################################################################################################
#' An S4 class to fit a user-defined distribution on a vector of marker values.
#'
#' This class allows to fit a user-defined distribution on the marker values \code{x}.
#' @slot x This slot takes in argument the marker values. Numeric argument.
#' @slot n Length of x vector (including NA values). Numeric argument.
#' @slot paraNames This slot is a character vector of distribution parameter names.
#' @slot mcmcList This slot is an mcmc.list object summing up all the sampled parameters of the user-defined distribution.
#' @slot cdf This slot is a function that describes the cumulative distribution function of the user-defined distribution.
#' @slot gradient This slot is a function that describes the probability density function of the user-defined distribution.
#' @slot hessian This slot is a function that describes the fisrt derivative of the probability density function of the user-defined distribution.
#' @slot mcmc This slot allows the main function to k,now whether an MCMC algorithm must be performed to sample the distribution parameters from their posterior distribution.
#' @details This class is automatically created when the user applies the \code{fit} function with the argument \code{distr="user"}. You never have to create manually this class, it is created internally.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a user-defined distribution.
#' @name fitUserDefinedDist-class
#' @aliases fitUserDefinedDist
#' @exportClass fitUserDefinedDist
setClass("fitUserDefinedDist", slots = c(x = "numeric", n = "numeric", paraNames = "character", mcmcList = "mcmc.list", cdf = "function", gradient = "function", hessian = "function", mcmc = "logical"))
setValidity("fitUserDefinedDist",
	function(object) {
		if (length(object@x) == 0) return("You cannot ignore 'x' argument")
		if (length(object@n) == 0) return("You cannot ignore 'n' argument")
		if (length(object@paraNames) == 0) return("You cannot ignore 'paraNames' argument")
		if (length(object@mcmcList) == 0) return("You cannot ignore 'mcmcList' argument")
		if (length(object@cdf) == 0) return("You cannot ignore 'cdf' argument")
		if (length(object@gradient) == 0) return("You cannot ignore 'gradient' argument")
		if (length(object@hessian) == 0) return("You cannot ignore 'hessian' argument")
		if (length(object@mcmc) == 0) return("You cannot ignore 'mcmc' argument")
		if (!(methods::is(object@cdf, "MethodDefinition"))) if (!(all(object@paraNames %in% names(formals(object@cdf))))) return("Arguments of 'cdf' function must include parameter names specified in 'paraNames'.")
		if (!(methods::is(object@gradient, "MethodDefinition"))) if (!(all(object@paraNames %in% names(formals(object@gradient))))) return("Arguments of 'gradient' function must include parameter names specified in 'paraNames'.")
		if (!(methods::is(object@hessian, "MethodDefinition"))) if (!(all(object@paraNames %in% names(formals(object@hessian))))) return("Arguments of 'hessian' function must include parameter names specified in 'paraNames'.")
		if (!(all(sapply(object@mcmcList, function(i) identical(sort(coda::varnames(i)), sort(object@paraNames)))))) return("Error: names in 'paraNames' argument do not match with names of 'mcmcList' argument.")
		lMethods <- c(NA,NA,NA)
		if (methods::is(object@cdf, "MethodDefinition")) {
			lMethods[1] <- attributes(object@cdf)$target[[1]]
			if (length(object@paraNames) != length(methods::slotNames(attributes(object@cdf)$target[[1]]))) return(paste0("Error: length of 'paraNames' is different from length of the parameters defined in an object ", attributes(object@cdf)$target[[1]], ". Maybe you did not choose the right cdf method?"))
			if (!(all(sapply(object@mcmcList, function(i) coda::varnames(i) %in% methods::slotNames(attributes(object@cdf)$target[[1]]))))) return(paste("Error: names in 'paraNames' must be", paste(methods::slotNames(attributes(object@cdf)$target[[1]]), collapse = ", and "), "to match with the definition of an object", attributes(object@cdf)$target[[1]]))
		}
		if (methods::is(object@gradient, "MethodDefinition")) {
			lMethods[2] <- attributes(object@gradient)$target[[1]]
			if (length(object@paraNames) != length(methods::slotNames(attributes(object@gradient)$target[[1]]))) return(paste0("Error: length of 'paraNames' is different from length of the parameters defined in an object ", attributes(object@gradient)$target[[1]], ". Maybe you did not choose the right gradient method?"))
			if (!(all(sapply(object@mcmcList, function(i) coda::varnames(i) %in% methods::slotNames(attributes(object@gradient)$target[[1]]))))) return(paste("Error: names in 'paraNames' must be", paste(methods::slotNames(attributes(object@gradient)$target[[1]]), collapse = ", and "), "to match with the definition of an object", attributes(object@gradient)$target[[1]]))
		}
		if (methods::is(object@hessian, "MethodDefinition")) {
			lMethods[3] <- attributes(object@hessian)$target[[1]]
			if (length(object@paraNames) != length(methods::slotNames(attributes(object@hessian)$target[[1]]))) return(paste0("Error: length of 'paraNames' is different from length of the parameters defined in an object ", attributes(object@hessian)$target[[1]], ". Maybe you did not choose the right hessian method?"))
			if (!(all(sapply(object@mcmcList, function(i) coda::varnames(i) %in% methods::slotNames(attributes(object@hessian)$target[[1]]))))) return(paste("Error: names in 'paraNames' must be", paste(methods::slotNames(attributes(object@hessian)$target[[1]]), collapse = ", and "), "to match with the definition of an object", attributes(object@hessian)$target[[1]]))
		}
		if (!(length(unique(lMethods[!is.na(lMethods)])) %in% c(0, 1))) return("You cannot use methods of different distribution types.")
		return(TRUE)
	}
)
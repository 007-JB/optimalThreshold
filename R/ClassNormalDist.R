#################################################################################################
###                                  ClassNormalDist.R                                        ###
#################################################################################################
#' An S4 class to represent a normal distribution.
#'
#' This S4 class describes the normal distribution that is fitted to the marker values. The normal distribution is characterized by the \code{mu} and the \code{sd} parameters.
#' @slot mu mu parameter.
#' @slot sd standard deviation parameter. Must be strictly positive.
#' @details You never have to create this class manually. This class is created internally when a normal distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a normal distribution.
#' @name normalDist-class
#' @aliases normalDist
#' @exportClass normalDist
setClass("normalDist", slots = c(mu = "numeric", sd = "numeric"))
setValidity("normalDist",
	function(object) {
		if (object@sd <= 0) return("Standard deviation lower or equal to zero.")
		else return(TRUE)
	}
)
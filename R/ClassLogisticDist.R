#################################################################################################
###                                  ClassLogisticDist.R                                      ###
#################################################################################################
#' An S4 class to represent a logistic distribution.
#'
#' This S4 class describes the logistic distribution that is fitted to the marker values. The logistic distribution is characterized by the \code{location} and the \code{scale} parameters.
#' @slot location location parameter.
#' @slot scale scale parameter. Must be strictly positive.
#' @details You never have to create this class manually. This class is created internally when a logistic distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a logistic distribution.
#' @name logisticDist-class
#' @aliases logisticDist
#' @exportClass logisticDist
setClass("logisticDist", slots = c(location = "numeric", scale = "numeric"))
setValidity("logisticDist",
	function(object) {
		if (object@scale <= 0) return("'scale' parameter must be positive.")
		return(TRUE)
	}
)
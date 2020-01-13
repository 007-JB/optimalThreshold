#################################################################################################
###                                  ClassGammaDist.R                                         ###
#################################################################################################
#' An S4 class to represent a gamma distribution.
#'
#' This S4 class describes the gamma distribution that is fitted to the marker values. The gamma distribution is characterized by the \code{shape} and the \code{scale} parameters.
#' @slot shape shape parameter. Must be positive.
#' @slot scale scale parameter. Must be strictly positive.
#' @details You never have to create this class manually. This class is created internally when a gamma distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a gamma distribution.
#' @name gammaDist-class
#' @aliases gammaDist
#' @exportClass gammaDist
setClass("gammaDist", slots = c(shape = "numeric", scale = "numeric"))
setValidity("gammaDist",
	function(object) {
		if (object@shape <= 0 | object@scale <= 0) return("'shape' and 'scale' parameters must be positive.")
		else return(TRUE)
	}
)
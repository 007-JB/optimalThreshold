#################################################################################################
###                                  ClassStudentDist.R                                       ###
#################################################################################################
#' An S4 class to represent a scaled Student distribution.
#'
#' This S4 class describes the scaled t distribution that is fitted to the marker values. The scaled t distribution is characterized by the \code{df} (degrees of freedom), the \code{mu}, and the \code{sd} parameters.
#' @slot df degrees of freedom (> 0, maybe non-integer). 
#' @slot mu mu parameter.
#' @slot sd standard deviation parameter. Must be strictly positive.
#' @details You never have to create this class manually. This class is created internally when a scaled t distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a t distribution.
#' @name studentDist-class
#' @aliases studentDist
#' @exportClass studentDist
setClass("studentDist", slots = c(df = "numeric", mu = "numeric", sd = "numeric"))
setValidity("studentDist",
	function(object) {
		if (object@sd <= 0) return("sd parameter must be positive.")
		if (object@df <= 0) return("df parameter must be positive.")
		return(TRUE)
	}
)
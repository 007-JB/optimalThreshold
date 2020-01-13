#################################################################################################
###                                  ClassLogNormalDist.R                                     ###
#################################################################################################
#' An S4 class to represent a log-normal distribution.
#'
#' This S4 class describes the log-normal distribution that is fitted to the marker values. The log-normal distribution is characterized by the \code{mu} and the \code{sd} parameters.
#' @slot mu mu parameter. 
#' @slot sd standard deviation parameter. Must be strictly positive.
#' @details You never have to create this class manually. This class is created internally when a log-normal distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit a log-normal distribution.
#' @name logNormalDist-class
#' @aliases logNormalDist
#' @exportClass logNormalDist
setClass("logNormalDist", slots = c(mu = "numeric", sd = "numeric"))
setValidity("logNormalDist",
	function(object) {
		if (object@sd <= 0) return("Standard deviation lower than zero.")
		else return(TRUE)
	}
)
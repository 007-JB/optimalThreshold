#################################################################################################
###                                  gradient.R                                               ###
#################################################################################################
#' Probability density function of a specified distribution
#'
#' The \code{gradient} function returns the probability density function relative to the S4 object passed in its argument. See details to know on what kind of S4 objects this function could be applied.
#' @name gradient
#' @param object Any S4 object for which a \code{gradient} method is defined. Should match with the definition of an S4 distribution object as defined in the \code{optimalThreshold} package.
#' @details This method can be applied to the S4 distribution objects that are supported in the \code{optimalThreshold} package: \code{normalDist}, \code{logNormalDist}, \code{gammaDist}, \code{studentDist}, \code{logisticDist}, and \code{userDefinedDist}. These methods are applied internally, and you have no need to use it outside of the main functions \code{trtSelThresh} and \code{diagThresh}. 
#' \itemize{
#' \item Normal distribution: the \code{gradient} method applied to a \code{normalDist} object is simply the \code{dnorm} function (see help on this function to have more details).
#' \item Log-normal distribution: the \code{gradient} method applied to a \code{logNormalDist} object is simply the \code{dlnorm} function (see help on this function to have more details).
#' \item Gamma distribution: the \code{gradient} method applied to a \code{gammaDist} object is simply the \code{dgamma} function (see help on this function to have more details).
#' \item Scaled t distribution: the scaled t distribution with \code{df} = n, \code{mu} = \eqn{\mu}, and \code{sd} = \eqn{\sigma} has density:
#' \deqn{f(x)=(\Gamma((n+1)/2)/(\sqrt{n\pi}\Gamma(n/2))(1+((x-\mu)/\sigma)^2/n)^-((n+1)/2))/\sigma}
#' \item Logistic distribution: the \code{gradient} method applied to a \code{logisticDist} object is simply the \code{dlogis} function (see help on this function to have more details).
#' \item User-defined distribution: the \code{gradient} method applied to a \code{userDefinedDist} object is simply the gradient function provided by the user when fitting a user-defined distribution with the \code{fit} function.
#' }
#' The S4 objects \code{compoundEvtRefDist}, \code{compoundNoEvtRefDist}, \code{compoundEvtInnovDist}, and \code{compoundNoEvtInnovDist} are created internally. The \code{gradient} function applied to these objects is defined dynamically depending on what types of distribution are fitted. The definition of the \code{gradient} function relies on the expression of the randomization constraint of a clinical trial that enforces the distribution of the marker in each treatment arm to be identical (see References for more details).
#' @return Returns the probability density function of the specified distribution.
#' @section References:
#' Blangero, Y, Rabilloud, M, Ecochard, R, and Subtil, F. A Bayesian method to estimate the optimal threshold of a marker used to select patients' treatment. \emph{Statistical Methods in Medical Research}. 2019.
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}}, \code{\link[stats]{dnorm}}, \code{\link[stats]{dlnorm}}, \code{\link[stats]{dgamma}}, \code{\link[stats]{dlogis}}, \code{\link[optimalThreshold]{fit}}
#' @exportMethod gradient
setGeneric(name = "gradient",
	def = function(object) {standardGeneric("gradient")}
)

#' @rdname gradient
#' @aliases gradient.normalDist
setMethod("gradient", "normalDist",
	function(object) {
		function(x) stats::dnorm(x, object@mu, object@sd)
	}
)

#' @rdname gradient
#' @aliases gradient.logNormalDist
setMethod("gradient", "logNormalDist",
	function(object) {
		function(x) stats::dlnorm(x, object@mu, object@sd)
	}
)

#' @rdname gradient
#' @aliases gradient.gammaDist
setMethod("gradient", "gammaDist",
	function(object) {
		function(x) stats::dgamma(x, object@shape, 1 / object@scale)
	}
)

#' @rdname gradient
#' @aliases gradient.studentDist
setMethod("gradient", "studentDist",
	function(object) {
		function(x) dt.scaled(x, object@df, object@mu, object@sd)
	}
)

#' @rdname gradient
#' @aliases gradient.logisticDist
setMethod("gradient", "logisticDist",
	function(object) {
		function(x) stats::dlogis(x, object@location, object@scale)
	}
)

#' @rdname gradient
#' @aliases gradient.compoundEvtRefDist
setMethod("gradient", "compoundEvtRefDist",
	function(object) {
		function(x) (gradient(object@EvtInnovDist)(x) * object@r1 + gradient(object@NoEvtInnovDist)(x) * (1 - object@r1) - gradient(object@NoEvtRefDist)(x) * (1 - object@r0)) / object@r0
	}
)

#' @rdname gradient
#' @aliases gradient.compoundNoEvtRefDist
setMethod("gradient", "compoundNoEvtRefDist",
	function(object) {
		function(x) (gradient(object@EvtInnovDist)(x) * object@r1 + gradient(object@NoEvtInnovDist)(x) * (1 - object@r1) - gradient(object@EvtRefDist)(x) * object@r0) / (1 - object@r0)
	}
)

#' @rdname gradient
#' @aliases gradient.compoundEvtInnovDist
setMethod("gradient", "compoundEvtInnovDist",
	function(object) {
		function(x) (gradient(object@EvtRefDist)(x) * object@r0 + gradient(object@NoEvtRefDist)(x) * (1 - object@r0) - gradient(object@NoEvtInnovDist)(x) * (1 - object@r1)) / object@r1
	}
)

#' @rdname gradient
#' @aliases gradient.compoundNoEvtInnovDist
setMethod("gradient", "compoundNoEvtInnovDist",
	function(object) {
		function(x) (gradient(object@EvtRefDist)(x) * object@r0 + gradient(object@NoEvtRefDist)(x) * (1 - object@r0) - gradient(object@EvtInnovDist)(x) * object@r1) / (1 - object@r1)
	}
)







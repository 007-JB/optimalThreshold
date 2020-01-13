#################################################################################################
###                                  cdf.R                                                    ###
#################################################################################################
#' Cumulative distribution function of a specified distribution
#'
#' The \code{cdf} function returns the cumulative distribution function relative to the S4 object passed in its argument. See details to know on what kind of S4 objects this function could be applied.
#' @name cdf
#' @param object Any S4 object for which a \code{cdf} method is defined. Should match with the definition of an S4 distribution object as defined in the \code{optimalThreshold} package.
#' @details This method can be applied to the S4 distribution objects that are supported in the \code{optimalThreshold} package: \code{normalDist}, \code{logNormalDist}, \code{gammaDist}, \code{studentDist}, \code{logisticDist}, and \code{userDefinedDist}. These methods are applied internally, and you have no need to use it outside of the main functions \code{trtSelThresh} and \code{diagThresh}. 
#' \itemize{
#' \item Normal distribution: the \code{cdf} method applied to a \code{normalDist} object is simply the \code{pnorm} function (see help on this function to have more details).
#' \item Log-normal distribution: the \code{cdf} method applied to a \code{logNormalDist} object is simply the \code{plnorm} function (see help on this function to have more details).
#' \item Gamma distribution: the \code{cdf} method applied to a \code{gammaDist} object is simply the \code{pgamma} function (see help on this function to have more details).
#' \item Scaled t distribution: the scaled t distribution with \code{df} = n, \code{mu} = \eqn{\mu}, and \code{sd} = \eqn{\sigma} has density:
#' \deqn{f(x)=(\Gamma((n+1)/2)/(\sqrt{n\pi}\Gamma(n/2))(1+((x-\mu)/\sigma)^2/n)^-((n+1)/2))/\sigma}
#' \item Logistic distribution: the \code{cdf} method applied to a \code{logisticDist} object is simply the \code{plogis} function (see help on this function to have more details).
#' \item User-defined distribution: the \code{cdf} method applied to a \code{userDefinedDist} object is simply the cumulative distribution function provided by the user when fitting a user-defined distribution with the \code{fit} function.
#' }
#' The S4 objects \code{compoundEvtRefDist}, \code{compoundNoEvtRefDist}, \code{compoundEvtInnovDist}, and \code{compoundNoEvtInnovDist} are created internally. The \code{cdf} function applied to these objects is defined dynamically depending on what types of distribution are fitted. The definition of the \code{cdf} function relies on the expression of the randomization constraint of a clinical trial that enforces the distribution of the marker in each treatment arm to be identical (see References for more details).
#' @return Returns the cumulative distribution function of the specified distribution.
#' @section References:
#' Blangero, Y, Rabilloud, M, Ecochard, R, and Subtil, F. A Bayesian method to estimate the optimal threshold of a marker used to select patients' treatment. \emph{Statistical Methods in Medical Research}. 2019.
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}}, \code{\link[optimalThreshold]{diagThresh}}, \code{\link[stats]{pnorm}}, \code{\link[stats]{plnorm}}, \code{\link[stats]{pgamma}}, \code{\link[stats]{plogis}}, \code{\link[optimalThreshold]{fit}}.
#' @exportMethod cdf
setGeneric(name = "cdf",
	def = function(object) {standardGeneric("cdf")}
)

#' @rdname cdf
#' @aliases cdf.normalDist
setMethod("cdf", "normalDist",
	function(object) {
		function(x) stats::pnorm(x, object@mu, object@sd)
	}
)

#' @rdname cdf
#' @aliases cdf.logNormalDist
setMethod("cdf", "logNormalDist",
	function(object) {
		function(x) stats::plnorm(x, object@mu, object@sd)
	}
)

#' @rdname cdf
#' @aliases cdf.gammaDist
setMethod("cdf", "gammaDist",
	function(object) {
		function(x) stats::pgamma(x, object@shape, 1 / object@scale)
	}
)

#' @rdname cdf
#' @aliases cdf.studentDist
setMethod("cdf", "studentDist",
	function(object) {
		function(x) pt.scaled(x, object@df, object@mu, object@sd)
	}
)

#' @rdname cdf
#' @aliases cdf.logisticDist
setMethod("cdf", "logisticDist",
	function(object) {
		function(x) stats::plogis(x, object@location, object@scale)
	}
)

#' @rdname cdf
#' @aliases cdf.compoundEvtRefDist
setMethod("cdf", "compoundEvtRefDist",
	function(object) {
		function(x) (cdf(object@EvtInnovDist)(x) * object@r1 + cdf(object@NoEvtInnovDist)(x) * (1 - object@r1) - cdf(object@NoEvtRefDist)(x) * (1 - object@r0)) / object@r0
	}
)

#' @rdname cdf
#' @aliases cdf.compoundNoEvtRefDist
setMethod("cdf", "compoundNoEvtRefDist",
	function(object) {
		function(x) (cdf(object@EvtInnovDist)(x) * object@r1 + cdf(object@NoEvtInnovDist)(x) * (1 - object@r1) - cdf(object@EvtRefDist)(x) * object@r0) / (1 - object@r0)
	}
)

#' @rdname cdf
#' @aliases cdf.compoundEvtInnovDist
setMethod("cdf", "compoundEvtInnovDist",
	function(object) {
		function(x) (cdf(object@EvtRefDist)(x) * object@r0 + cdf(object@NoEvtRefDist)(x) * (1 - object@r0) - cdf(object@NoEvtInnovDist)(x) * (1 - object@r1)) / object@r1
	}
)

#' @rdname cdf
#' @aliases cdf.compoundNoEvtInnovDist
setMethod("cdf", "compoundNoEvtInnovDist",
	function(object) {
		function(x) (cdf(object@EvtRefDist)(x) * object@r0 + cdf(object@NoEvtRefDist)(x) * (1 - object@r0) - cdf(object@EvtInnovDist)(x) * object@r1) / (1 - object@r1)
	}
)
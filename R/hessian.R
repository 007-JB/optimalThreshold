#################################################################################################
###                                  hessian.R                                                ###
#################################################################################################
#' Second derivative of the cumulative distribution function of a specified distribution
#'
#' The \code{hessian} function returns the second derivative of the cumulative distribution function relative to the S4 object passed in its argument. See details to know on what kind of S4 objects this function could be applied.
#' @name hessian
#' @param object A distribution object.
#' @details This method can be applied to the S4 distribution objects that are supported in the \code{optimalThreshold} package: \code{normalDist}, \code{logNormalDist}, \code{gammaDist}, \code{studentDist}, \code{logisticDist}, and \code{userDefinedDist}. These methods are applied internally, and you have no need to use it outside of the main functions \code{trtSelThresh} and \code{diagThresh}. 
#' \itemize{
#' \item Normal distribution: the \code{hessian} method applied to a \code{normalDist} object is simply the second derivative of the cumulative distribution function of a normal distribution, with \code{mu}=\eqn{\mu} and \code{sd}=\eqn{\sigma}, and expressed as: 
#' \deqn{f'(x)=((\mu-x)/\sigma^2)*f(x)}
#' \item Log-normal distribution: the \code{hessian} method applied to a \code{logNormalDist} object is simply the second derivative of the cumulative distribution function of a log-normal distribution, with \code{mu}=\eqn{\mu} and \code{sd}=\eqn{\sigma}, and expressed as:
#' \deqn{f'(x)=(((\mu-\log(x))/(x*\sigma^2))-1/x)*f(x)}
#' \item Gamma distribution: the \code{hessian} method applied to a \code{gammaDist} object is simply the second derivative of the cumulative distribution function of a gamma distribution, with \code{shape}=\eqn{\alpha} and \code{scale}=\eqn{\beta}, and expressed as:
#' \deqn{f'(x)=((\alpha-1)/x-1/\beta)*f(x)}
#' \item Scaled t distribution: the \code{hessian} method applied to a \code{studentDist} object is simply the second derivative of the cumulative distribution function of a t scaled distribution, with \code{df}=n, \code{mu}=\eqn{\mu} and \code{sd}=\eqn{\sigma}, and expressed as:
#' \deqn{f'(x)=(-(n+1))*((x-\mu)/(\sigma^2*(n+((x-\mu)/\sigma)^2)))*f(x)}
#' \item Logistic distribution: the \code{hessian} method applied to a \code{logisticDist} object is simply the second derivative of the cumulative distribution function of a logistic distribution, with \code{location}=\eqn{\mu}, and \code{scale}=\eqn{\sigma}, and expressed as:
#' \deqn{f'(x)=((\exp(-(x-\mu)/\sigma)^2-1)/(\sigma*(1+\exp(-(x-\mu)/\sigma))^2))*f(x)}
#' \item User-defined distribution: the \code{hessin} method applied to a \code{userDefinedDist} object is simply the hessian function provided by the user when fitting a user-defined distribution with the \code{fit} function.
#' }
#' The S4 objects \code{compoundEvtRefDist}, \code{compoundNoEvtRefDist}, \code{compoundEvtInnovDist}, and \code{compoundNoEvtInnovDist} are created internally. The \code{hessian} function applied to these objects is defined dynamically depending on what types of distribution are fitted. The definition of the \code{hessian} function relies on the expression of the randomization constraint of a clinical trial that enforces the distribution of the marker in each treatment arm to be identical (see References for more details).
#' @return Returns the second derivative of the cumulative distribution function of the specified distribution.
#' @section References:
#' Blangero, Y, Rabilloud, M, Ecochard, R, and Subtil, F. A Bayesian method to estimate the optimal threshold of a marker used to select patients' treatment. \emph{Statistical Methods in Medical Research}. 2019.
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}}, \code{\link[optimalThreshold]{fit}}
#' @exportMethod hessian
setGeneric(name = "hessian",
	def = function(object) {standardGeneric("hessian")}
) 

#' @rdname hessian
#' @aliases hessian.normalDist
setMethod("hessian", "normalDist",
	function(object) {
		function(x) ((object@mu - x) / object@sd^2) * stats::dnorm(x, object@mu, object@sd)
	}
)

#' @rdname hessian
#' @aliases hessian.logNormalDist
setMethod("hessian", "logNormalDist",
	function(object) {
		function(x) stats::dlnorm(x, object@mu, object@sd) * (((object@mu - log(x)) / (x * object@sd^2)) - (1 / x))
	}
)

#' @rdname hessian
#' @aliases hessian.gammaDist
setMethod("hessian", "gammaDist",
	function(object) {
		function(x) stats::dgamma(x, object@shape, 1 / object@scale) * ((object@shape - 1) / x - 1 / object@scale)
	}
)

#' @rdname hessian
#' @aliases hessian.studentDist
setMethod("hessian", "studentDist",
	function(object) {
		function(x) dt.scaled(x, object@df, object@mu, object@sd) * (- (object@df + 1) * (x - object@mu)) / (object@sd^2 * (object@df + ((x - object@mu) / object@sd)^2))
	}
)

#' @rdname hessian
#' @aliases hessian.logisticDist
setMethod("hessian", "logisticDist",
	function(object) {
		function(x) stats::dlogis(x, object@location, object@scale) * (exp(- (x - object@location) / object@scale)^2 - 1) / (object@scale * (1 + exp(- (x - object@location) / object@scale))^2)
	}
)

#' @rdname hessian
#' @aliases hessian.compoundEvtRefDist
setMethod("hessian", "compoundEvtRefDist",
	function(object) {
		function(x) (hessian(object@EvtInnovDist)(x) * object@r1 + hessian(object@NoEvtInnovDist)(x) * (1 - object@r1) - hessian(object@NoEvtRefDist)(x) * (1 - object@r0)) / object@r0
	}
)

#' @rdname hessian
#' @aliases hessian.compoundNoEvtRefDist
setMethod("hessian", "compoundNoEvtRefDist",
	function(object) {
		function(x) (hessian(object@EvtInnovDist)(x) * object@r1 + hessian(object@NoEvtInnovDist)(x) * (1 - object@r1) - hessian(object@EvtRefDist)(x) * object@r0) / (1 - object@r0)
	}
)

#' @rdname hessian
#' @aliases hessian.compoundEvtInnovDist
setMethod("hessian", "compoundEvtInnovDist",
	function(object) {
		function(x) (hessian(object@EvtRefDist)(x) * object@r0 + hessian(object@NoEvtRefDist)(x) * (1 - object@r0) - hessian(object@NoEvtInnovDist)(x) * (1 - object@r1)) / object@r1
	}
)

#' @rdname hessian
#' @aliases hessian.compoundNoEvtInnovDist
setMethod("hessian", "compoundNoEvtInnovDist",
	function(object) {
		function(x) (hessian(object@EvtRefDist)(x) * object@r0 + hessian(object@NoEvtRefDist)(x) * (1 - object@r0) - hessian(object@EvtInnovDist)(x) * object@r1) / (1 - object@r1)
	}
)
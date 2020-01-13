#################################################################################################
###                                  ClassDiagOptThresh.R                                     ###
#################################################################################################

#' An S4 class to describe the optimal threshold of a diagnostic marker.
#'
#' @slot optThresh This slot is an object that takes in argument the sampled optimal threshold values. Numeric argument.
#' @slot p Prevalence specified by the user. Numeric argument.
#' @slot r Risk threshold preference. Numeric argument.
#' @slot N Sample size.
#' @slot xEvt This slot is an object that takes in argument the marker values in the subgroup of patients that developed the event. Numeric argument.
#' @slot xNoEvt This slot is an object that takes in argument the marker values in the subgroup of patients that did not develop the event. Numeric argument.
#' @slot lowEvt This slot is a logical argument that specifies whether the low values of the marker are associated with the presence of the disease or not.
#' @slot mcmcChainEvt This slot is an object that takes in argument the sampled distribution objects in the subgroup of patients that developed the event. list argument.
#' @slot mcmcChainNoEvt This slot is an object that takes in argument the sampled distribution objects in the subgroup of patients that did not develop the event. list argument.
#' @slot tabMCMCChain This slot is an object that takes in argument all the distribution parameters that were sampled using the MCMC algorithm. mcmc.listOrNull argument.
#' @slot paraNamesUserDefined This slot is an object that takes in argument the list of the distribution parameter names defined by the user in a 'fitUserDefinedDist' object. list argument.
#' @slot cdfUserDefined This slot is an object that takes in argument the list of cumulative distribution functions defined by the user in 'fitUserDefinedDist' objects. list argument.
#' @slot gradientUserDefined This slot is an object that takes in argument the list of gradient functions defined by the user in 'fitUserDefinedDist' objects. list argument.
#' @slot hessianUserDefined This slot is an object that takes in argument the list of hessian functions defined by the user in 'fitUserDefinedDist' objects. list argument.
#' @slot percentNA This slot is a numeric object that indicates the percentage of NA values contained in the 'optThresh' slot.
#' @details You never have to create this class manually. This class is created internally when the \code{diagThresh} function is used.
#' @seealso \code{\link[optimalThreshold]{diagThresh}} for more details on how to estimate the optimal threshold of a diagnostic marker.
#' @name diagOptThresh-class
#' @aliases diagOptThresh
#' @exportClass diagOptThresh
setClass("diagOptThresh", slots = c(optThresh = "numeric", p = "numeric", r = "numeric", N = "numeric", xEvt = "numeric", xNoEvt = "numeric", lowEvt = "logical", mcmcChainEvt = "list", mcmcChainNoEvt = "list", tabMCMCChain = "mcmc.listOrNull", paraNamesUserDefined = "list", cdfUserDefined = "list", gradientUserDefined = "list", hessianUserDefined = "list", percentNA = "numeric"))

#' @name summary-methods
#' @aliases summary.diagOptThresh summary,diagOptThresh-method
#' @return This function returns an object of class 'summaryDiagOptThresh'.
#' @details For a \code{diagOpthThresh} object, it prints:
#' \itemize{
#' 	\item The decision rule: is the reference treatment recommended for low values of the marker?
#'	\item The median (default), mean, or mode risk of event occurrence in each treatment arm, and their credible interval.
#'	\item Some summary statistics of the marker under study (min, max, quartiles and mean)
#'	\item The optimal threshold estimate and its credible interval (percentile and highest posterior density).
#'	\item The median (default), mean, or mode risk in each arm under the marker-based strategy.
#'	\item The median (default), mean, or mode benefit estimate under each treatment arm.
#'	\item The percentage of NA values returned during the optimal threshold estimation process.
#'}
#' @section References :
#' Subtil, F, and Rabilloud. A Bayesian method to estimate the optimal threshold of a longitudinal biomarker. \emph{Biometrical Journal}. 2010.
#' @seealso \code{\link[optimalThreshold]{diagThresh}} for more details on how to estimate the optimal threshold of a diagnostic marker.
#' @exportMethod summary
setMethod("summary", "diagOptThresh",
	function(object, alpha = 0.05, method = "median") {
		if (!(method %in% c("median", "mean", "mode"))) stop("'method' must be mean, median or mode.")
		if (object@lowEvt) cat("Decision rule: Low values of the marker are associated with the presence of the event.\n")
		else cat("Decision rule: High values of the marker are associated with the presence of the event.\n")
		cat(paste0("\nThe prevalence of the event is ", object@p, ".\n"))
		if (method == "median") thresholdEst <- stats::median(object@optThresh, na.rm = TRUE)
		else {
			if (method == "mean") thresholdEst <- mean(object@optThresh, na.rm = TRUE)
			else thresholdEst <- Mode(object@optThresh, na.rm = TRUE)
		}
		markerRange <- summary(c(object@xEvt, object@xNoEvt))
		thresholdEstIC <- stats::quantile(object@optThresh, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		thresholdEstICHPD <- HDInterval::hdi(object@optThresh, credMass = 1 - alpha)
		attributes(thresholdEstICHPD) <- NULL
		cat("\nOptimal threshold estimate:\n")
		print(thresholdEst)
		cat("\nCredible intervals of the optimal threshold\n")
		cat(paste("- Percentile method (", (1 - alpha) * 100, "%)", " \n", sep = ""))
		print(unname(thresholdEstIC))
		cat(paste("- Highest Posterior Density (", (1 - alpha) * 100, "%)", " \n", sep = ""))
		print(thresholdEstICHPD)
		cat(paste("\nPercentage of NA values returned during the estimation process: ", object@percentNA * 100, "%\n", sep = ""))
	}
)

#' @rdname decisionCurve
#' @aliases decisionCurve,diagOptThresh-method decisionCurve.diagOptThresh
setMethod("decisionCurve", "diagOptThresh",
	function(object, r) {
		if (any(r < 0 | r >= 1)) stop("r must be in [0;1[.")
		if (!(is.null(object@cdfUserDefined$Evt))) {
			.optimalThreshold_envEvtDist <- new.env()
			paraTypes <- rep("numeric", length(object@paraNamesUserDefined$Evt))
			names(paraTypes) <- object@paraNamesUserDefined$Evt
			methods::setClass("userDefinedDist0E", slots = paraTypes, where = .optimalThreshold_envEvtDist)
			methods::setMethod("cdf", "userDefinedDist0E",
				object@cdfUserDefined$Evt,
				where = .optimalThreshold_envEvtDist
			)
			methods::setMethod("gradient", "userDefinedDist0E",
				object@gradientUserDefined$Evt,
				where = .optimalThreshold_envEvtDist
			)
			methods::setMethod("hessian", "userDefinedDist0E",
				object@hessianUserDefined$Evt,
				where = .optimalThreshold_envEvtDist
			)
		}
		if (!(is.null(object@cdfUserDefined$NoEvt))) {
			.optimalThreshold_envNoEvtDist <- new.env()
			paraTypes <- rep("numeric", length(object@paraNamesUserDefined$NoEvt))
			names(paraTypes) <- object@paraNamesUserDefined$NoEvt
			methods::setClass("userDefinedDist0Eb", slots = paraTypes, where = .optimalThreshold_envNoEvtDist)
			methods::setMethod("cdf", "userDefinedDist0Eb",
				object@cdfUserDefined$NoEvt,
				where = .optimalThreshold_envNoEvtDist
			)
			methods::setMethod("gradient", "userDefinedDist0Eb",
				object@gradientUserDefined$NoEvt,
				where = .optimalThreshold_envNoEvtDist
			)
			methods::setMethod("hessian", "userDefinedDist0Eb",
				object@hessianUserDefined$NoEvt,
				where = .optimalThreshold_envNoEvtDist
			)
		}
		if (object@lowEvt) coefU <- (- 1)
		else coefU <- 1
		mcmcChainEvt <- unlist(object@mcmcChainEvt)
		mcmcChainNoEvt <- unlist(object@mcmcChainNoEvt)
		if (coefU == 1) {
			diagEB <- function(x, r, prev, EvtDist, NoEvtDist, coefU) {
				U <- (- ((1 - cdf(EvtDist)(x)) * prev - (1 - cdf(NoEvtDist)(x)) * (1 - prev) * (r / (1 - r))))
				attr(U, 'gradient') <- (- (- gradient(EvtDist)(x) * prev + gradient(NoEvtDist)(x) * (1 - prev) * (r / (1 - r))))
				attr(U, 'hessian') <- (- (- hessian(EvtDist)(x) * prev + hessian(NoEvtDist)(x) * (1 - prev) * (r / (1 - r))))
				U
			}
		}
		else {
			diagEB <- function(x, r, prev, EvtDist, NoEvtDist, coefU) {
				U <- (- (cdf(EvtDist)(x) * prev - cdf(NoEvtDist)(x) * (1 - prev) * (r / (1 - r))))
				attr(U, 'gradient') <- (- (gradient(EvtDist)(x) * prev - gradient(NoEvtDist)(x) * (1 - prev) * (r / (1 - r))))
				attr(U, 'hessian') <- (- (hessian(EvtDist)(x) * prev - hessian(NoEvtDist)(x) * (1 - prev) * (r / (1 - r))))
				U
			}
		}
		defineStartVal2 <- function(gridStartVal, utility, ...) {
			Us <- utility(x = gridStartVal, ...)
			startVal <- gridStartVal[which(Us == min(Us, na.rm = TRUE))][1]
			return(startVal)
		}
		maxDiagUtilityValue <- function(EvtDist, NoEvtDist, prev, r, minMarker, maxMarker, utility, coefU) {
			gridStartVal <- seq(minMarker, maxMarker, length.out = 1000)
			startVal <- defineStartVal2(gridStartVal, utility, r, prev, EvtDist, NoEvtDist, coefU)
			thresholdEst <- stats::nlm(utility, startVal, r = r, prev = prev, EvtDist = EvtDist, NoEvtDist = NoEvtDist, coefU = coefU)
			return(- thresholdEst$minimum)
		}
		minMarker <- min(c(object@xEvt, object@xNoEvt))
		maxMarker <- max(c(object@xEvt, object@xNoEvt))
		Uc <- sapply(r, function(i) mapply(maxDiagUtilityValue, EvtDist = mcmcChainEvt, NoEvtDist = mcmcChainNoEvt, MoreArgs = list(prev = object@p, r = i, minMarker = minMarker, maxMarker = maxMarker, utility = diagEB, coefU = coefU)))
		if (object@lowEvt) {
			UNoTreat <- sapply(r, function(r) object@p - (1 - object@p) * (r / (1 - r)))
			UTreatAll <- rep(0, length(r))
		}
		else {
			UNoTreat <- rep(0, length(r))
			UTreatAll <- sapply(r, function(r) object@p - (1 - object@p) * (r / (1 - r)))
		}
		res <- methods::new("diagRelUtility", r = r, U = Uc, UNoTreat = as.matrix(UNoTreat), UTreatAll = as.matrix(UTreatAll))
		plot(res)
		if (!(is.null(object@cdfUserDefined$Evt))) {
			methods::removeMethod(cdf, "userDefinedDist0E", where = .optimalThreshold_envEvtDist)
			methods::removeMethod(gradient, "userDefinedDist0E", where = .optimalThreshold_envEvtDist)
			methods::removeMethod(hessian, "userDefinedDist0E", where = .optimalThreshold_envEvtDist)
			methods::removeClass("userDefinedDist0E", where = .optimalThreshold_envEvtDist)			
		}
		if (!(is.null(object@cdfUserDefined$NoEvt))) {
			methods::removeMethod(cdf, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
			methods::removeMethod(gradient, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
			methods::removeMethod(hessian, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
			methods::removeClass("userDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)			
		}
		invisible(res)
	}
)

#' estimates method
#'
#' @rdname estimates
#' @aliases estimates,diagOptThresh-method estimates.diagOptThresh
setMethod("estimates", "diagOptThresh",
	function(object, method = "median") {
		if (!(method %in% c("median", "mean", "mode"))) stop("'method' must be mean, median or mode.")
		if (method == "median") thresholdEst <- stats::median(object@optThresh, na.rm = TRUE)
		if (method == "mean") thresholdEst <- mean(object@optThresh, na.rm = TRUE)
		if(method == "mode") thresholdEst <- Mode(object@optThresh, na.rm = TRUE)
		return(c(optThresh = thresholdEst))
	}
)

#' @rdname credibleIntervals
#' @aliases credibleIntervals,diagOptThresh-method credibleIntervals.diagOptThresh
setMethod("credibleIntervals", "diagOptThresh",
	function(object, alpha = 0.05, hpd = FALSE) {
		if (alpha < 0 | alpha > 1) stop("'alpha' parameter must be in [0;1].")
		if (!(is.logical(hpd))) stop("'hpd' parameter must be logical.")
		if (hpd) thresholdEstIC <- HDInterval::hdi(object@optThresh, credMass = 1 - alpha)
		else thresholdEstIC <- stats::quantile(object@optThresh, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		return(thresholdEstIC)
	}
)

#' @rdname plot-methods
#' @aliases plot,diagOptThresh-method plot.diagOptThresh
#' @exportMethod plot
setMethod("plot", "diagOptThresh",
	function(x, y, main = "MCMC sample distribution of optimal threshold", col = "gray85", border.col = "darkgrey", xlab = "Optimal threshold estimate", yaxs = "i", freq = FALSE, breaks = seq(min(x@optThresh, na.rm = TRUE), max(x@optThresh, na.rm = TRUE), length.out = 20), ...) {
		opar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(opar))
		graphics::par(mfrow = c(1, 1), lwd = 2)
		graphics::hist(x@optThresh, yaxs = yaxs, main = main, col = col, border = border.col, xlab = xlab, freq = freq, breaks = breaks, ...)
		colAlpha <- t(grDevices::col2rgb("gray85", alpha = TRUE))
		colAlpha[4] <- 160
		graphics::polygon(c(stats::density(x@optThresh, na.rm = TRUE)$x, rev(stats::density(x@optThresh, na.rm = TRUE)$x)), c(stats::density(x@optThresh, na.rm = TRUE)$y, rep(0, length(stats::density(x@optThresh, na.rm = TRUE)$x))), col = grDevices::rgb(colAlpha[1], colAlpha[2], colAlpha[3], colAlpha[4], maxColorValue = 255), border = NA)
		graphics::lines(stats::density(x@optThresh, na.rm = TRUE), col = border.col, lwd = 2)
		graphics::box(bty = "L", lwd = 1)
	}
)

#' @name show-methods
#' @aliases show,diagOptThresh-method show.diagOptThresh
#' @exportMethod show
setMethod("show", "diagOptThresh",
	function(object) {
		cat("object@optThresh (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@optThresh))
		print(object@optThresh[1:lengthShow])
		if (lengthShow < length(object@optThresh)) cat("...\n")
		cat("\nobject@p\n")
		print(object@p)
		cat("\nobject@r\n")
		print(object@r)
		cat("\nobject@xEvt (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@xEvt))
		print(object@xEvt[1:lengthShow])
		if (lengthShow<length(object@xEvt)) cat("...\n")
		cat("\nobject@xNoEvt (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@xNoEvt))
		print(object@xNoEvt[1:lengthShow])
		if (lengthShow<length(object@xNoEvt)) cat("...\n")
		cat("\nobject@mcmcChainEvt (limited to the first value of each MCMC chain)\n")
		print(lapply(object@mcmcChainEvt, function(x) x[1]))
		cat("\nobject@mcmcChainNoEvt (limited to the first value of each MCMC chain)\n")
		print(lapply(object@mcmcChainNoEvt, function(x) x[1]))
		cat("\nobject@paraNamesUserDefined\n")
		print(object@paraNamesUserDefined)
		cat("\nobject@cdfUserDefined\n")
		print(object@cdfUserDefined)
		cat("\nobject@gradientUserDefined\n")
		print(object@gradientUserDefined)
		cat("\nobject@hessianUserDefined\n")
		print(object@hessianUserDefined)
	}
)
#################################################################################################
###                                  ClassTrtSelOptThresh.R                                   ###
#################################################################################################

#' An S4 class to describe the optimal threshold of a treatment selection marker.
#'
#' @slot optThresh This slot is an object that takes in argument the sampled optimal threshold values. Numeric argument.
#' @slot r0 This slot is an object that takes in argument the sampled mean risks of event occurrence in the reference arm. Numeric argument.
#' @slot r1 This slot is an object that takes in argument the sampled mean risks of event occurrence in the innovative arm. Numeric argument.
#' @slot xEvtRef This slot is an object that takes in argument the marker values in the subgroup of patients that developed the event in the reference arm. Numeric argument.
#' @slot xNoEvtRef This slot is an object that takes in argument the marker values in the subgroup of patients that did not develop the event in the reference arm. Numeric argument.
#' @slot xEvtInnov This slot is an object that takes in argument the marker values in the subgroup of patients that developed the event in the innovative arm. Numeric argument.
#' @slot xNoEvtInnov This slot is an object that takes in argument the marker values in the subgroup of patients that did not develop the event in the innovative arm. Numeric argument.
#' @slot lowRef This slot is a logical argument that specifies whether the reference treatment is recommended for low values of the marker.
#' @slot toxRef This slot is a logical argument that specifies whether the reference treatment is the most toxic treatment option at equal efficacy with the innovative treatment.
#' @slot markerBasedRiskRef This slot is an object that takes in argument the sampled mean risks of event occurrence in the reference treatment under the marker-based allocation rule. Numeric argument.
#' @slot markerBasedRiskInnov This slot is an object that takes in argument the sampled mean risks of event occurrence in the innovative treatment under the marker-based allocation rule. Numeric argument.
#' @slot mcmcChainEvtRef This slot is an object that takes in argument the sampled distribution objects in the subgroup of patients that developed the event in the reference arm. list argument.
#' @slot mcmcChainNoEvtRef This slot is an object that takes in argument the sampled distribution objects in the subgroup of patients that did not develop the event in the reference arm. list argument.
#' @slot mcmcChainEvtInnov This slot is an object that takes in argument the sampled distribution objects in the subgroup of patients that developed the event in the innovative arm. list argument.
#' @slot mcmcChainNoEvtInnov This slot is an object that takes in argument the sampled distribution objects in the subgroup of patients that did not develop the event in the innovative arm. list argument.
#' @slot tabMCMCChain This slot is an object that takes in argument all the distribution parameters that were sampled using the MCMC algorithm. mcmc.listOrNull argument.
#' @slot paraNamesUserDefined This slot is an object that takes in argument the list of the distribution parameter names defined by the user in a 'fitUserDefinedDist' object. list argument.
#' @slot cdfUserDefined This slot is an object that takes in argument the list of cumulative distribution functions defined by the user in 'fitUserDefinedDist' objects. list argument.
#' @slot gradientUserDefined This slot is an object that takes in argument the list of gradient functions defined by the user in 'fitUserDefinedDist' objects. list argument.
#' @slot hessianUserDefined This slot is an object that takes in argument the list of hessian functions defined by the user in 'fitUserDefinedDist' objects. list argument.
#' @slot percentNA This slot is a numeric object that indicates the percentage of NA values contained in the 'optThresh' slot.
#' @details You never have to create this class manually. This class is created internally when the \code{trtSelThresh} function is used.
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}} for more details on how to estimate the optimal threshold of a treatment selection marker.
#' @name trtSelOptThresh-class
#' @aliases trtSelOptThresh
#' @exportClass trtSelOptThresh
setClass("trtSelOptThresh", slots = c(optThresh = "numeric", r0 = "numeric", r1 = "numeric", xEvtRef = "numeric", xNoEvtRef = "numeric", xEvtInnov = "numeric", xNoEvtInnov = "numeric", lowRef = "logical", toxRef = "logical", markerBasedRiskRef = "numeric", markerBasedRiskInnov = "numeric", mcmcChainEvtRef = "list", mcmcChainNoEvtRef = "list", mcmcChainEvtInnov = "list", mcmcChainNoEvtInnov = "list", tabMCMCChain = "mcmc.listOrNull", paraNamesUserDefined = "list", cdfUserDefined = "list", gradientUserDefined = "list", hessianUserDefined = "list", percentNA = "numeric"))

#' An S4 method that summarizes the results of a \code{trtSelOpthThresh} or a \code{diagOpthThresh} object.
#'
#' @name summary-methods
#' @aliases summary.trtSelOptThresh summary,trtSelOptThresh-method
#' @param object a \code{trtSelOptThresh} S4 class object for which a summary is desired.
#' @param alpha alpha parameter for the confidence level required.
#' @param method which method to use: median, mean or mode (median is the default).
#' @return This function returns an object of class 'summaryTrtSelOptThresh'.
#' @details This function presents the results stocked in a \code{trtSelOpthThresh} object, or in a \code{diagOpthThresh} object. For a \code{trtSelOpthThresh} object it prints:
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
#' Blangero, Y, Rabilloud, M, Ecochard, R, and Subtil, F. A Bayesian method to estimate the optimal threshold of a marker used to select patients' treatment. \emph{Statistical Methods in Medical Research}. 2019.
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}} for more details on how to estimate the optimal threshold of a treatment selection marker.
#' @exportMethod summary
setMethod("summary", "trtSelOptThresh",
	function(object, alpha = 0.05, method = "median") {
		if (!(method %in% c("median", "mean", "mode"))) stop("'method' must be mean, median or mode.")
		if (object@lowRef) cat("Decision rule: The reference treatment is preferred for low values of the marker.\n")
		else cat("Decision rule: The innovative treatment is preferred for low values of the marker.\n")
		if (method == "median") {
			thresholdEst <- stats::median(object@optThresh, na.rm = TRUE)
			r0 <- stats::median(object@r0, na.rm = TRUE)
			r1 <- stats::median(object@r1, na.rm = TRUE)
			r0x <- stats::median(object@markerBasedRiskRef, na.rm = TRUE)
			r1x <- stats::median(object@markerBasedRiskInnov, na.rm = TRUE)
			benefitRef <- stats::median(object@r0 - object@markerBasedRiskRef, na.rm = TRUE)
			benefitInnov <- stats::median(object@r1 - object@markerBasedRiskInnov, na.rm = TRUE)
		}
		else {
			if (method == "mean") {
				thresholdEst <- mean(object@optThresh, na.rm = TRUE)
				r0 <- mean(object@r0, na.rm = TRUE)
				r1 <- mean(object@r1, na.rm = TRUE)
				r0x < -mean(object@markerBasedRiskRef, na.rm = TRUE)
				r1x <- mean(object@markerBasedRiskInnov, na.rm = TRUE)
				benefitRef <- mean(object@r0 - object@markerBasedRiskRef, na.rm = TRUE)
				benefitInnov <- mean(object@r1 - object@markerBasedRiskInnov, na.rm = TRUE)
			}
			else {
				thresholdEst <- Mode(object@optThresh, na.rm = TRUE)
				r0 <- Mode(object@r0, na.rm = TRUE)
				r1 <- Mode(object@r1, na.rm = TRUE)
				r0x <- Mode(object@markerBasedRiskRef, na.rm = TRUE)
				r1x <- Mode(object@markerBasedRiskInnov, na.rm = TRUE)
				benefitRef <- Mode(object@r0 - object@markerBasedRiskRef, na.rm = TRUE)
				benefitInnov <- Mode(object@r1 - object@markerBasedRiskInnov, na.rm = TRUE)
			}
		}
		markerRange <- summary(c(object@xEvtRef, object@xNoEvtRef, object@xEvtInnov, object@xNoEvtInnov))
		thresholdEstIC <- stats::quantile(object@optThresh, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		thresholdEstICHPD <- HDInterval::hdi(object@optThresh, credMass = 1 - alpha)
		attributes(thresholdEstICHPD) <- NULL
		r0IC <- stats::quantile(object@r0, probs = c(alpha / 2, 1 - alpha / 2))
		r1IC <- stats::quantile(object@r1, probs = c(alpha / 2, 1 - alpha / 2))
		r0xIC <- stats::quantile(object@markerBasedRiskRef, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		r1xIC <- stats::quantile(object@markerBasedRiskInnov, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		benefitRefIC <- stats::quantile(object@r0 - object@markerBasedRiskRef, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		benefitInnovIC <- stats::quantile(object@r1 - object@markerBasedRiskInnov, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		if (method == "median") {
			cat(paste("\nMedian risk in the reference group: ", round(r0, 2), " [", round(r0IC[[1]], 2), "; ", round(r0IC[[2]], 2), "]\n", sep = ""))
			cat(paste("Median risk in the innovative group: ", round(r1, 2), " [", round(r1IC[[1]], 2), "; ", round(r1IC[[2]], 2), "]\n", sep = ""))
		}
		if (method == "mean") {
			cat(paste("\nMean risk in the reference group: ", round(r0, 2), " [", round(r0IC[[1]], 2), "; ", round(r0IC[[2]], 2), "]\n", sep = ""))
			cat(paste("Mean risk in the innovative group: ", round(r1, 2), " [", round(r1IC[[1]], 2), "; ", round(r1IC[[2]], 2), "]\n", sep = ""))
		}
		if (method == "mode") {
			cat(paste("\nMode risk in the reference group: ", round(r0, 2), " [", round(r0IC[[1]], 2), "; ", round(r0IC[[2]], 2), "]\n", sep = ""))
			cat(paste("Mode risk in the innovative group: ", round(r1, 2), " [", round(r1IC[[1]], 2), "; ", round(r1IC[[2]], 2), "]\n", sep = ""))
		}
		cat("\nSummary statistics of the marker:\n")
		print(markerRange)
		cat("\nOptimal threshold estimate:\n")
		print(thresholdEst)
		cat("\nCredible intervals of the optimal threshold\n")
		cat(paste("- Percentile method (", (1 - alpha) * 100, "%)", " \n", sep = ""))
		print(unname(thresholdEstIC))
		cat(paste("- Highest Posterior Density (", (1 - alpha) * 100, "%)", " \n", sep = ""))
		print(thresholdEstICHPD)
		if (method == "median") {
			cat(paste("\nMedian risk in the reference group under marker-based strategy: ", round(r0x, 2), " [", round(r0xIC[[1]], 2), "; ", round(r0xIC[[2]], 2), "]\n", sep = ""))
			cat(paste("Median risk in the innovative group under marker-based strategy: ", round(r1x, 2), " [", round(r1xIC[[1]], 2), "; ", round(r1xIC[[2]], 2), "]\n", sep = ""))
			cat("\nMedian benefit of the marker-based strategy \n")
			cat(paste("- In the reference group: ", round(benefitRef, 2), " [", round(benefitRefIC[[1]], 2), "; ", round(benefitRefIC[[2]], 2), "]\n", sep = ""))
			cat(paste("- In the innovative group: ", round(benefitInnov, 2), " [", round(benefitInnovIC[[1]], 2), "; ", round(benefitInnovIC[[2]], 2), "]\n", sep = ""))
		}		
		if (method == "mean") {
			cat(paste("\nMean risk in the reference group under marker-based strategy: ", round(r0x, 2), " [", round(r0xIC[[1]], 2), "; ", round(r0xIC[[2]], 2), "]\n", sep = ""))
			cat(paste("Mean risk in the innovative group under marker-based strategy: ", round(r1x, 2), " [", round(r1xIC[[1]], 2), "; ", round(r1xIC[[2]], 2), "]\n", sep = ""))
			cat("\nMean benefit of the marker-based strategy \n")
			cat(paste("- In the reference group: ", round(benefitRef, 2), " [", round(benefitRefIC[[1]], 2), "; ", round(benefitRefIC[[2]], 2), "]\n", sep = ""))
			cat(paste("- In the innovative group: ", round(benefitInnov, 2), " [", round(benefitInnovIC[[1]], 2), "; ", round(benefitInnovIC[[2]], 2), "]\n", sep = ""))
		}
		if (method == "mode") {
			cat(paste("\nMode risk in the reference group under marker-based strategy: ", round(r0x, 2), " [", round(r0xIC[[1]], 2), "; ", round(r0xIC[[2]], 2), "]\n", sep = ""))
			cat(paste("Mode risk in the innovative group under marker-based strategy: ", round(r1x, 2), " [", round(r1xIC[[1]], 2), "; ", round(r1xIC[[2]], 2), "]\n", sep = ""))
			cat("\nMode benefit of the marker-based strategy \n")
			cat(paste("- In the reference group: ", round(benefitRef, 2), " [", round(benefitRefIC[[1]], 2), "; ", round(benefitRefIC[[2]], 2), "]\n", sep = ""))
			cat(paste("- In the innovative group: ", round(benefitInnov, 2), " [", round(benefitInnovIC[[1]], 2), "; ", round(benefitInnovIC[[2]], 2), "]\n", sep = ""))
		}
		cat(paste("\nPercentage of NA values returned during the estimation process: ", object@percentNA * 100, "%\n", sep = ""))
	}
)

#' Decision curve plot
#'
#' This S4 method allows to plot the decision curve associated with a treatment selection marker for several treatment/event cost ratios, or with a diagnostic marker for several risk thresholds.
#' @name decisionCurve
#' @param object a \code{trtSelOptThresh} or a \code{diagOptThresh} S4 class object from which the decision curve must be plotted.
#' @param r Ratio of treatment/event costs (for treatment selection markers), or risk preference (for diagnostic marker). Numeric argument.
#' @param ... other arguments passed to methods.
#' @return Returns an object of class \code{trtSelRelUtility} if applied to a \code{trtSelOptThresh} object, and an object of class \code{diagRelUtility} if applied to a \code{diagOptThresh} object.
#' @details This function plots the decision curves according to the type of marker under study (treatment selection or diagnostic) and the \code{r} argument. For treatment selection markers, it plots two graphs: the first one is a classical decision curves graph comparing the utilities of: the marker under study, the perfect marker, the strategy 'Treat everyone with the reference treatment', and the strategy 'Treat everyone with the innovative treatment'; the second one is the relative decision curve that plots the relative utility of the marker under study (0 meaning that using the marker to guide treatment decisions is not better than treating everyone with the overall best treatment, and 1 meaning that the marker under study has the same utility as the perfect marker). The decision curves are calculated using the mean risk of event in each treatment arm provided in the \code{trtSelOptThresh} object.
#' For diagnostic markers it calculates the expected benefit of the marker and compares it with the strategies "Treat everyone" and "Do not treat anyone". The decision curve is calculated for a population with a disease prevalence provided in the \code{diagOptThresh} object.
#' @section References:
#' Blangero, Y, Rabilloud, M, Ecochard, R, and Subtil, F. A Bayesian method to estimate the optimal threshold of a marker used to select patients' treatment. \emph{Statistical Methods in Medical Research}. 2019.
#' Huang, Y, Laber, EB, and Janes H. Characterizing expected benefits of biomarkers in treatment selection. \emph{Biostatistics}. 2015; 16(2): 383-399.
#' Vickers, AJ, and Elkin, EB. Decision curve analysis: a novel method for evaluating prediction models. \emph{Medical Decision Making}. 2006; 26(6): 565-574.
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}} and \code{\link[optimalThreshold]{diagThresh}} for some examples on how to use this function.
#' @exportMethod decisionCurve
setGeneric(name = "decisionCurve",
	def = function(object, r, ...) {standardGeneric("decisionCurve")}
)

#' @rdname decisionCurve
#' @aliases decisionCurve,trtSelOptThresh-method decisionCurve.trtSelOptThresh
#' @param alpha alpha parameter for the confidence level required.
setMethod("decisionCurve","trtSelOptThresh",
	function(object, r, alpha = 0.05) {
		if (any(r < 0 | r >= 1)) stop("r must be in [0;1[.")
		if (!(is.null(object@cdfUserDefined$EvtRef))) {
			.optimalThreshold_envEvtRefDist <- new.env()
			paraTypes <- rep("numeric", length(object@paraNamesUserDefined$EvtRef))
			names(paraTypes) <- object@paraNamesUserDefined$EvtRef
			methods::setClass("userDefinedDist0E", slots = paraTypes, where = .optimalThreshold_envEvtRefDist)
			methods::setMethod("cdf", "userDefinedDist0E",
				object@cdfUserDefined$EvtRef,
				where = .optimalThreshold_envEvtRefDist
			)
			methods::setMethod("gradient", "userDefinedDist0E",
				object@gradientUserDefined$EvtRef,
				where = .optimalThreshold_envEvtRefDist
			)
			methods::setMethod("hessian", "userDefinedDist0E",
				object@hessianUserDefined$EvtRef,
				where = .optimalThreshold_envEvtRefDist
			)
		}
		if (!(is.null(object@cdfUserDefined$NoEvtRef))) {
			.optimalThreshold_envNoEvtRefDist <- new.env()
			paraTypes <- rep("numeric", length(object@paraNamesUserDefined$NoEvtRef))
			names(paraTypes) <- object@paraNamesUserDefined$NoEvtRef
			methods::setClass("userDefinedDist0Eb", slots = paraTypes, where = .optimalThreshold_envNoEvtRefDist)
			methods::setMethod("cdf", "userDefinedDist0Eb",
				object@cdfUserDefined$NoEvtRef,
				where = .optimalThreshold_envNoEvtRefDist
			)
			methods::setMethod("gradient", "userDefinedDist0Eb",
				object@gradientUserDefined$NoEvtRef,
				where = .optimalThreshold_envNoEvtRefDist
			)
			methods::setMethod("hessian", "userDefinedDist0Eb",
				object@hessianUserDefined$NoEvtRef,
				where = .optimalThreshold_envNoEvtRefDist
			)
		}
		if (!(is.null(object@cdfUserDefined$EvtInnov))) {
			.optimalThreshold_envEvtInnovDist <- new.env()
			paraTypes <- rep("numeric", length(object@paraNamesUserDefined$EvtInnov))
			names(paraTypes) <- object@paraNamesUserDefined$EvtInnov
			methods::setClass("userDefinedDist1E", slots = paraTypes, where = .optimalThreshold_envEvtInnovDist)
			methods::setMethod("cdf", "userDefinedDist1E",
				object@cdfUserDefined$EvtInnov,
				where = .optimalThreshold_envEvtInnovDist
			)
			methods::setMethod("gradient", "userDefinedDist1E",
				object@gradientUserDefined$EvtInnov,
				where = .optimalThreshold_envEvtInnovDist
			)
			methods::setMethod("hessian", "userDefinedDist1E",
				object@hessianUserDefined$EvtInnov,
				where = .optimalThreshold_envEvtInnovDist
			)
		}
		if (!(is.null(object@cdfUserDefined$NoEvtInnov))) {
			.optimalThreshold_envNoEvtInnovDist <- new.env()
			paraTypes <- rep("numeric", length(object@paraNamesUserDefined$NoEvtInnov))
			names(paraTypes) <- object@paraNamesUserDefined$NoEvtInnov
			methods::setClass("userDefinedDist1Eb", slots = paraTypes, where = .optimalThreshold_envNoEvtInnovDist)
			methods::setMethod("cdf", "userDefinedDist1Eb",
				object@cdfUserDefined$NoEvtInnov,
				where = .optimalThreshold_envNoEvtInnovDist
			)
			methods::setMethod("gradient", "userDefinedDist1Eb",
				object@gradientUserDefined$NoEvtInnov,
				where = .optimalThreshold_envNoEvtInnovDist
			)
			methods::setMethod("hessian", "userDefinedDist1Eb",
				object@hessianUserDefined$NoEvtInnov,
				where = .optimalThreshold_envNoEvtInnovDist
			)
		}
		mcmcChainEvtRef <- unlist(object@mcmcChainEvtRef)
		mcmcChainNoEvtRef <- unlist(object@mcmcChainNoEvtRef)
		mcmcChainEvtInnov <- unlist(object@mcmcChainEvtInnov)
		mcmcChainNoEvtInnov <- unlist(object@mcmcChainNoEvtInnov)
		defineStartVal2 <- function(gridStartVal, utility, ...) {
			Us <- utility(x = gridStartVal, ...)
			startVal <- gridStartVal[which(Us == min(Us, na.rm = TRUE))][1]
			return(startVal)
		}
		maxUtilityValue <- function(EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist, r0, r1, pT0, pT1, r, coefU, minMarker, maxMarker, utility) {
			gridStartVal <- seq(minMarker, maxMarker, length.out = 1000)
			startVal <- defineStartVal2(gridStartVal, utility, r0, r1, r, pT0, pT1, coefU, EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist)
			thresholdEst <- stats::nlm(utility, startVal, r0 = r0, r1 = r1, r = r, pT0 = pT0, pT1 = pT1, coefU = coefU, EvtRefDist = EvtRefDist, NoEvtRefDist = NoEvtRefDist, EvtInnovDist = EvtInnovDist, NoEvtInnovDist = NoEvtInnovDist)
			return(- thresholdEst$minimum)
		}
		if (object@lowRef == TRUE & object@toxRef == TRUE) {
			Up <- sapply(r, function(r) (object@r1 - pmax(0, object@r0 + object@r1 - 1)) * (1 - r))
			UT1 <- sapply(r, function(r) object@r1 - object@r0 - r)
			UT0 <- 0
			r <- (-r)
			coefU <- 1
		}
		else {
			if (object@lowRef == TRUE & object@toxRef == FALSE) {
				Up <- sapply(r, function(r) (object@r1 - pmax(0, object@r0 + object@r1 - 1)) + (1 - object@r0 + pmax(0, object@r0 + object@r1 - 1)) * r)
				UT1 <- sapply(r, function(r) object@r1 - object@r0 + r)
				UT0 <- 0
				coefU <- 1
			}
			else {
				if (object@lowRef == FALSE & object@toxRef == TRUE) {
					Up <- sapply(r, function(r) (object@r0 - pmax(0, object@r0 + object@r1 - 1)) + (1 - object@r1 + pmax(0, object@r0 + object@r1 - 1)) * r)
					UT1 <- sapply(r, function(r) object@r0 - object@r1 + r)
					UT0 <- 0
					r <- (-r)
					coefU <- (-1)
				}
				else {
					Up <- sapply(r, function(r) (object@r0 - pmax(0, object@r0 + object@r1 - 1)) * (1 - r))
					UT1 <- sapply(r, function(r) object@r0 - object@r1 - r)
					UT0 <- 0
					coefU <- (-1)
				}
			}
		}
		pT0 <- length(c(object@xEvtRef, object@xNoEvtRef)) / length(c(object@xEvtRef, object@xNoEvtRef, object@xEvtInnov, object@xNoEvtInnov))
		pT1 <- length(c(object@xEvtInnov, object@xNoEvtInnov)) / length(c(object@xEvtRef, object@xNoEvtRef, object@xEvtInnov, object@xNoEvtInnov))
		minMarker <- min(c(object@xEvtRef, object@xNoEvtRef, object@xEvtInnov, object@xNoEvtInnov))
		maxMarker <- max(c(object@xEvtRef, object@xNoEvtRef, object@xEvtInnov, object@xNoEvtInnov))
		Uc <- sapply(r, function(i) mapply(maxUtilityValue, mcmcChainEvtRef, mcmcChainNoEvtRef, mcmcChainEvtInnov, mcmcChainNoEvtInnov, object@r0, object@r1, MoreArgs = list(r = i, pT0 = pT0, pT1 = pT1, coefU = coefU, minMarker = minMarker, maxMarker = maxMarker, utility = utility)))
		RU <- (Uc - pmax(UT0, UT1)) / (Up - pmax(UT0, UT1))
		if (any(apply(RU, 2, function(x) mean(x > 1 | x < 0)) > 0.05)) {
			warning("At least 5% of iterations lead to numerical instability and relative utilities greater than 1 or lower than 0 for some given r ratios.")
		} 
		r <- abs(r)
		res <- methods::new("trtSelRelUtility", RU = RU, r = r, U = Uc, UT0 = as.matrix(UT0), UT1 = as.matrix(UT1), Up = Up)
		plot(res, alpha = alpha)
		if (!(is.null(object@cdfUserDefined$EvtRef))) {
			methods::removeMethod(cdf, "userDefinedDist0E", where = .optimalThreshold_envEvtRefDist)
			methods::removeMethod(gradient, "userDefinedDist0E", where = .optimalThreshold_envEvtRefDist)
			methods::removeMethod(hessian, "userDefinedDist0E", where = .optimalThreshold_envEvtRefDist)
			methods::removeClass("userDefinedDist0E", where = .optimalThreshold_envEvtRefDist)			
		}
		if (!(is.null(object@cdfUserDefined$NoEvtRef))) {
			methods::removeMethod(cdf, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
			methods::removeMethod(gradient, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
			methods::removeMethod(hessian, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
			methods::removeClass("userDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)			
		}
		if (!(is.null(object@cdfUserDefined$EvtInnov))) {
			methods::removeMethod(cdf, "userDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
			methods::removeMethod(gradient, "userDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
			methods::removeMethod(hessian, "userDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
			methods::removeClass("userDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)			
		}
		if (!(is.null(object@cdfUserDefined$NoEvtInnov))) {
			methods::removeMethod(cdf, "userDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
			methods::removeMethod(gradient, "userDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
			methods::removeMethod(hessian, "userDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
			methods::removeClass("userDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)			
		}
		invisible(res)
	}
)

#' Indicator estimates
#'
#' This function calculates the punctual estimates of several indicators.
#' @name estimates
#' @param object a \code{trtSelOptThresh} S4 class object from which the decision curve must be plotted.
#' @param ... other arguments to pass to estimates method.
#' @return Returns a list of several indicator estimates.
#' @details This function calculates the punctual estimates of several indicators (median, mean, or mode) depending on the type of marker under study (treatment selection or diagnostic marker).
#' @return Returns the estimates of several indicators.
#' @exportMethod estimates
setGeneric(name = "estimates",
	def = function(object, ...) {standardGeneric("estimates")}
)

#' @rdname estimates
#' @param method which method to use: median, mean or mode (median is the default).
#' @aliases estimates,trtSelOptThresh-method estimates.trtSelOptThresh
setMethod("estimates", "trtSelOptThresh",
	function(object, method = "median") {
		if (!(method %in% c("median", "mean", "mode"))) stop("'method' must be mean, median or mode.")
		if (method == "median") {
			thresholdEst <- stats::median(object@optThresh, na.rm = TRUE)
			r0 <- stats::median(object@r0, na.rm = TRUE)
			r1 <- stats::median(object@r1, na.rm = TRUE)
			r0x <- stats::median(object@markerBasedRiskRef, na.rm = TRUE)
			r1x <- stats::median(object@markerBasedRiskInnov, na.rm = TRUE)
			benefitRef <- stats::median(object@r0 - object@markerBasedRiskRef, na.rm = TRUE)
			benefitInnov <- stats::median(object@r1 - object@markerBasedRiskInnov, na.rm = TRUE)
		}
		if (method == "mean") {
			thresholdEst <- mean(object@optThresh, na.rm = TRUE)
			r0 <- mean(object@r0, na.rm = TRUE)
			r1 <- mean(object@r1, na.rm = TRUE)
			r0x <- mean(object@markerBasedRiskRef, na.rm = TRUE)
			r1x <- mean(object@markerBasedRiskInnov, na.rm = TRUE)
			benefitRef <- mean(object@r0 - object@markerBasedRiskRef, na.rm = TRUE)
			benefitInnov <- mean(object@r1 - object@markerBasedRiskInnov, na.rm = TRUE)
		}
		if (method == "mode") {
			thresholdEst <- Mode(object@optThresh, na.rm = TRUE)
			r0 <- Mode(object@r0, na.rm = TRUE)
			r1 <- Mode(object@r1, na.rm = TRUE)
			r0x <- Mode(object@markerBasedRiskRef, na.rm = TRUE)
			r1x <- Mode(object@markerBasedRiskInnov, na.rm = TRUE)
			benefitRef <- Mode(object@r0 - object@markerBasedRiskRef, na.rm = TRUE)
			benefitInnov <- Mode(object@r1 - object@markerBasedRiskInnov, na.rm = TRUE)
		}
		return(list(optThresh = thresholdEst, risks = c(Reference = r0, Innovative = r1), markerBasedRisks = c(Reference = r0x, Innovative = r1x), benefits = c(Reference = benefitRef, Innovative = benefitInnov)))
	}
)

#' Credible intervals estimation
#'
#' @name credibleIntervals
#' @param object a \code{trtSelOptThresh} or a \code{diagOptThresh} S4 class object from which the credible intervals of several indicators (including the optimal threshold) must be calculated.
#' @param ... other arguments passed to the method.
#' @return Returns a matrix with the credible intervals of several indicators.
#' @details This function calculates the credible intervals of several indicators depending on the type of marker under study (treatment selection or diagnostic marker). The user may specify the alpha risk for the confidence level (default is 5%), and the method to calculate the credible intervals: percentile (default with \code{hpd = FALSE}) or using the highest posterior density method (with \code{hpd = TRUE}). When \code{hpd = TRUE}, the function \code{hdi} from the \code{HDInterval} package is used.
#' @section References:
#' Gelman, A, et al. 2014. \emph{Bayesian Data Analysis}. 3rd edition, CRC Press, Boca Raton, section 2.3.
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}}, \code{\link[optimalThreshold]{diagThresh}}, and \code{\link[HDInterval]{hdi}}.
#' @return Returns the credible intervals of several indicators.
#' @exportMethod credibleIntervals
setGeneric(name = "credibleIntervals",
	def = function(object, ...) {standardGeneric("credibleIntervals")}
) 

#' @rdname credibleIntervals
#' @aliases credibleIntervals,trtSelOptThresh-method credibleIntervals.trtSelOptThresh
#' @param alpha alpha parameter for the confidence level required.
#' @param hpd logical value to specify whether the function has to return Highest Posterior Density interval or not.
setMethod("credibleIntervals", "trtSelOptThresh",
	function(object, alpha = 0.05, hpd = FALSE) {
		if (alpha < 0 | alpha > 1) stop("'alpha' parameter must be in [0;1].")
		if (!(is.logical(hpd))) stop("'hpd' parameter must be logical.")
		if (hpd) {
			thresholdEstIC <- HDInterval::hdi(object@optThresh, credMass = 1 - alpha)
			r0IC <- HDInterval::hdi(object@r0, credMass = 1 - alpha)
			r1IC <- HDInterval::hdi(object@r1, credMass = 1 - alpha)
			r0xIC <- HDInterval::hdi(object@markerBasedRiskRef, credMass = 1 - alpha)
			r1xIC <- HDInterval::hdi(object@markerBasedRiskInnov, credMass = 1 - alpha)
			benefitRefIC <- HDInterval::hdi(object@r0 - object@markerBasedRiskRef, credMass = 1 - alpha)
			benefitInnovIC <- HDInterval::hdi(object@r1 - object@markerBasedRiskInnov, credMass = 1 - alpha)
		}
		else {
			thresholdEstIC <- stats::quantile(object@optThresh, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
			r0IC <- stats::quantile(object@r0, probs = c(alpha / 2, 1 - alpha / 2))
			r1IC <- stats::quantile(object@r1, probs = c(alpha / 2, 1 - alpha / 2))
			r0xIC <- stats::quantile(object@markerBasedRiskRef, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
			r1xIC <- stats::quantile(object@markerBasedRiskInnov, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
			benefitRefIC <- stats::quantile(object@r0 - object@markerBasedRiskRef, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
			benefitInnovIC <- stats::quantile(object@r1 - object@markerBasedRiskInnov, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		}
		ans <- rbind(thresholdEstIC, r0IC, r1IC, r0xIC, r1xIC, benefitRefIC, benefitInnovIC)
		rownames(ans) <- c("OptThresh CI", "Reference risk CI", "Innovative risk CI", "Reference marker-based risk CI", "Innovative marker-based risk CI", "Reference benefit CI", "Innovative benefit CI")
		return(ans)
	}
)

#' Plot method
#'
#' @name plot-methods
#' @aliases plot,trtSelOptThresh-method plot.trtSelOptThresh
#' @param x a \code{trtSelOptThresh} or a \code{diagOptThresh} object.
#' @param y unused parameter.
#' @param main an overall title for the plot.
#' @param col the color of the histogram.
#' @param border.col the color of the histogram border.
#' @param xlab a label for the x axis of the plot.
#' @param yaxs The style of axis interval calculation to be used for the y-axis.
#' @param freq logical; if TRUE, the histogram graphic is a representation of frequencies; if FALSE, probability densities are plotted (so that the histogram has a total area of one). 
#' @param breaks one of: 
#' \itemize{
#' \item a vector giving the breakpoints between histogram cells, 
#' \item a function to compute the vector of breakpoints, 
#' \item a single number giving the number of cells for the histogram, 
#' \item a character string naming an algorithm to compute the number of cells, 
#' \item a function to compute the number of cells. 
#' }
#' In the last three cases the number is a suggestion only; as the breakpoints will be set to pretty values, the number is limited to 1e6 (with a warning if it was larger). If breaks is a function, the x vector is supplied to it as the only argument (and the number of breaks is only limited by the amount of available memory). 
#' @param ... other graphical parameters.
#' @return None
#' @exportMethod plot
setMethod("plot", "trtSelOptThresh",
	function(x, y, main = "MCMC sample distribution of optimal threshold", col = "gray85", border.col = "darkgrey", xlab = "Optimal threshold estimate", yaxs = "i", freq = FALSE, breaks = seq(min(x@optThresh, na.rm = TRUE), max(x@optThresh, na.rm = TRUE), length.out = 20), ...) {
		opar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(opar))
		graphics::par(mfrow = c(1, 1), lwd = 2)
		graphics::hist(x@optThresh,yaxs = yaxs, main = main, col = col, border = border.col, xlab = xlab, freq = freq, breaks = breaks, ...)
		colAlpha <- t(grDevices::col2rgb("gray85", alpha = TRUE))
		colAlpha[4] <- 160
		graphics::polygon(c(stats::density(x@optThresh, na.rm = TRUE)$x, rev(stats::density(x@optThresh, na.rm = TRUE)$x)), c(stats::density(x@optThresh, na.rm = TRUE)$y, rep(0, length(stats::density(x@optThresh, na.rm = TRUE)$x))), col = grDevices::rgb(colAlpha[1], colAlpha[2], colAlpha[3], colAlpha[4], maxColorValue = 255), border = NA)
		graphics::lines(stats::density(x@optThresh, na.rm = TRUE), col = border.col, lwd = 2)
		graphics::box(bty = "L", lwd = 1)
	}
)

#' Show method
#'
#' Show some of the slots of a \code{trtSelOptThresh} or a \code{diagOptThresh} objects.
#' @name show-methods
#' @aliases show,trtSelOptThresh-method show.trtSelOptThresh
#' @param object a \code{trtSelOptThresh} or a \code{diagOptThresh} S4 object.
#' @return None
#' @exportMethod show
setMethod("show", "trtSelOptThresh",
	function(object) {
		cat("object@optThresh (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@optThresh))
		print(object@optThresh[1:lengthShow])
		if (lengthShow < length(object@optThresh)) cat("...\n")
		cat("\nobject@r0 (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@r0))
		print(object@r0[1:lengthShow])
		if (lengthShow < length(object@r0)) cat("...\n")
		cat("\nobject@r1 (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@r1))
		print(object@r1[1:lengthShow])
		if (lengthShow < length(object@r1)) cat("...\n")
		cat("\nobject@xEvtRef (limited to first 10 values)\n")
		lengthShow <- min(10,length(object@xEvtRef))
		print(object@xEvtRef[1:lengthShow])
		if (lengthShow < length(object@xEvtRef)) cat("...\n")
		cat("\nobject@xNoEvtRef (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@xNoEvtRef))
		print(object@xNoEvtRef[1:lengthShow])
		if (lengthShow < length(object@xNoEvtRef)) cat("...\n")
		cat("\nobject@xEvtInnov (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@xEvtInnov))
		print(object@xEvtInnov[1:lengthShow])
		if (lengthShow < length(object@xEvtInnov)) cat("...\n")
		cat("\nobject@xNoEvtInnov (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@xNoEvtInnov))
		print(object@xNoEvtInnov[1:lengthShow])
		if (lengthShow < length(object@xNoEvtInnov)) cat("...\n")
		cat("\nobject@lowRef\n")
		print(object@lowRef)
		cat("\nobject@toxRef\n")
		print(object@toxRef)
		cat("\nobject@markerBasedRiskRef (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@markerBasedRiskRef))
		print(object@markerBasedRiskRef[1:lengthShow])
		if (lengthShow < length(object@markerBasedRiskRef)) cat("...\n")
		cat("\nobject@markerBasedRiskInnov (limited to first 10 values)\n")
		lengthShow <- min(10, length(object@markerBasedRiskInnov))
		print(object@markerBasedRiskInnov[1:lengthShow])
		if (lengthShow < length(object@markerBasedRiskInnov)) cat("...\n")
		cat("\nobject@mcmcChainEvtRef (limited to the first value of each MCMC chain)\n")
		print(lapply(object@mcmcChainEvtRef, function(x) x[1]))
		cat("\nobject@mcmcChainNoEvtRef (limited to the first value of each MCMC chain)\n")
		print(lapply(object@mcmcChainNoEvtRef, function(x) x[1]))
		cat("\nobject@mcmcChainEvtInnov (limited to the first value of each MCMC chain)\n")
		print(lapply(object@mcmcChainEvtInnov, function(x) x[1]))
		cat("\nobject@mcmcChainNoEvtInnov (limited to the first value of each MCMC chain)\n")
		print(lapply(object@mcmcChainNoEvtInnov, function(x) x[1]))
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
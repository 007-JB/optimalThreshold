#################################################################################################
###                                  samplePosteriorDist.R                                    ###
#################################################################################################
#' Sample in the posterior distribution of the parameters of a given theoretical distribution.
#'
#' The \code{samplePosteriorDist} function samples the parameters of a given theoretical distribution using explicit posterior distribution (if it exists), or a Markov Chain Monte Carlo (MCMC) algorithm when the posterior distribution is unknown. See details to know on what kind of S4 objects this function could be applied.
#' @name samplePosteriorDist
#' @param object A distribution object.
#' @param K A numerical value indicating the length of the sample.
#' @param ... other parameters passed to methods.
#' @details This method can be applied to the S4 distribution objects that are supported in the \code{optimalThreshold} package: \code{fitNormalDist}, \code{fitLogNormalDist}, \code{fitGammaDist}, \code{fitStudentDist}, and \code{fitLogisticDist}. These methods are applied internally, and you have no need to use it outside of the main function \code{optThresEst}. See below to have details on the expression of the \code{samplePosteriorDist} function according to the type of distribution.
#' \itemize{
#' \item Normal distribution: a noninformative prior is used for the parameters of the normal distribution (\code{mu}=\eqn{\mu}, and \code{sd}=\eqn{\sigma}). The \eqn{\sigma^2} parameter is sampled from an inverse Chi-squared distribution, and the \eqn{\mu} parameter is sampled from a normal distribution with known variance. So, sampling in the posterior distribution of \eqn{\mu} and \eqn{\sigma} does not involve an MCMC algorithm (see References for more details and justification).
#' \item Log-normal distribution: a noninformative prior is used for the parameters of the log-normal distribution (\code{mu}=\eqn{\mu}, and \code{sd}=\eqn{\sigma}). The \eqn{\sigma^2} parameter is sampled from an inverse Chi-squared distribution, and the \eqn{\mu} parameter is sampled from a normal distribution with known variance. So, sampling in the posterior distribution of \eqn{\mu} and \eqn{\sigma} does not involve an MCMC algorithm (see References for more details and justification).
#' \item Gamma distribution: a noninformative prior is used for the parameters of the gamma distribution (\code{shape}=\eqn{\alpha}, and \code{scale}=\eqn{\beta}). The parameters are sampled using an adaptive rejection sampling (ARS) algorithm. The \eqn{\beta} parameter is sampled at the first iteration from an inverse gamma distribution using the initial value of the \eqn{\alpha} parameter provided by the user. Then the ARS algorithm is performed to sample \eqn{\alpha} from its posterior distribution (see References for more details and justification).
#' \item Scaled t distribution: a vague prior is used for the parameters of the scaled t distribution as a default. However, the user can write its own JAGS model to use different priors (see the \code{fit} function for more details). Sampling from the posterior distribution of the parameters of a scaled t distribution requires JAGS to be installed.
#' \item Logistic distribution: a vague prior is used for the parameters of the logistic distribution as a default. However, the user can write its own JAGS model to use different priors (see the \code{fit} function for more details). Sampling from the posterior distribution of the parameters of a logistic distribution requires JAGS to be installed.
#' }
#' @return Returns an object of class list.
#' @section References:
#' Gelman, A, et al. 2014. \emph{Bayesian Data Analysis}. 3rd edition, CRC Press, Boca Raton, section 2.8.
#' Sook, Y, and Oh, M. Bayesian estimation of the two-parameter Gamma distribution. \emph{Communications in Statistics - Simulation and Computation}. 2006; 35: 285-293.
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}}
#' @exportMethod samplePosteriorDist
setGeneric(name = "samplePosteriorDist",
	def = function(object, K, ...) {standardGeneric("samplePosteriorDist")}
)

#' @rdname samplePosteriorDist
#' @aliases samplePosteriorDist.fitNormalDist
#' @param n number of MCMC chains.
setMethod("samplePosteriorDist", "fitNormalDist",
	function(object, K, n) {
		mc <- function(object, K) {
			N <- length(object@x)
			standDev <- sqrt((N - 1) * stats::var(object@x) / stats::rchisq(K, N - 1))
			mu <- stats::rnorm(K, mean(object@x), standDev / sqrt(N))
			listPara <- mapply(function(mu, sd) methods::new("normalDist", mu = mu, sd = sd), mu = mu, sd = standDev)
			listPara
		}
		replicate(n, mc(object, K), simplify = FALSE)
	}
)

#' @rdname samplePosteriorDist
#' @aliases samplePosteriorDist.fitLogNormalDist
setMethod("samplePosteriorDist", "fitLogNormalDist",
	function(object, K, n) {
		mc <- function(object, K) {
			N <- length(object@x)
			standDev <- sqrt((N - 1) * stats::var(log(object@x)) / stats::rchisq(K, N - 1))
			mu <- stats::rnorm(K, mean(log(object@x)), standDev / sqrt(N))
			listPara <- mapply(function(mu, sd) methods::new("logNormalDist", mu = mu, sd = sd), mu = mu, sd = standDev)
			listPara
		}
		replicate(n, mc(object, K), simplify = FALSE)
	}
)

#' @rdname samplePosteriorDist
#' @aliases samplePosteriorDist.fitGammaDist
#' @param do.pb Indicates whther progressing bar or not
#' @param seed seed for the random number generator. Integer.
setMethod("samplePosteriorDist", "fitGammaDist",
	function(object, K, do.pb, seed) {
		nchains <- length(object@ini)
		scaleMaxL <- stats::var(object@x) / mean(object@x)
		mcmc <- function(i, ini, object, K, scaleMaxL, pb, nchains) {
			N <- length(object@x)
			if (!is.null(seed)) set.seed(seed + (nchains - (nchains + 1 - i)))
			shapeVec <- vector(length = K * object@thin + object@burnin + 1)
			scaleVec <- vector(length = K * object@thin + object@burnin)
			shapeVec[1] <- ini$shape
			gridArs <- seq(min(object@x) / scaleMaxL, max(object@x) / scaleMaxL, length.out = 10)
			if (!(is.null(pb))) {
				for (j in 1:(K * object@thin + object@burnin)) {
					utils::setTxtProgressBar(pb, (i - 1) * (K * object@thin + object@burnin) + j)
					scaleVec[j] <- (1 / stats::rgamma(1, N * shapeVec[j], scale = 1 / sum(object@x)))
					shapeVec[j+1] <- ars::ars(n = 1, f = fInArs, fprima = fPrimaInArs, xVec = object@x, scale = scaleVec[j], x = gridArs, m = 10, lb = TRUE, xlb = 0)
				}
			}
			else {
				for (j in 1:(K * object@thin + object@burnin)) {
					scaleVec[j] <- (1 / stats::rgamma(1, N * shapeVec[j], scale = 1 / sum(object@x)))
					shapeVec[j+1] <- ars::ars(n = 1, f = fInArs, fprima = fPrimaInArs, xVec = object@x, scale = scaleVec[j], x = gridArs, m = 10, lb = TRUE, xlb = 0)
				}
			}
			shapeVec <- shapeVec[(object@burnin + 2):(K * object@thin + object@burnin + 1)]
			shapeVec <- shapeVec[seq(1, length(shapeVec), by = object@thin)]
			scaleVec <- scaleVec[(object@burnin + 1):(K * object@thin + object@burnin)]
			scaleVec <- scaleVec[seq(1, length(scaleVec), by = object@thin)]
			listPara <- mapply(function(shape, scale) methods::new("gammaDist", shape = shape, scale = scale), shape = shapeVec, scale = scaleVec)
			listPara
		}
		if (do.pb) {
			cat("\nSampling\n\n")
			pb <- utils::txtProgressBar(1, nchains * (K * object@thin + object@burnin), initial = 1, style = 3, width = 50, char = "*")
		}
		else pb <- NULL
		res <- mapply(mcmc, 1:nchains, object@ini, MoreArgs = list(object = object, K = K, scaleMaxL = scaleMaxL, pb = pb, nchains = nchains), SIMPLIFY = FALSE)
		if (do.pb) close(pb)
		res
	}
)

#' @rdname samplePosteriorDist
#' @aliases samplePosteriorDist.fitStudentDist
setMethod("samplePosteriorDist", "fitStudentDist",
	function(object, K, do.pb, seed) {
		nchains <- length(object@ini)
		convertIni <- function(l) {
			l$sd <- log(l$sd)
			names(l)[names(l) == "sd"] <- "log_sd"
			l$df <- 1 / l$df
			names(l)[names(l) == "df"] <- "inv_df"
			l
		}
		inits <- lapply(object@ini, convertIni)
		if (!(is.null(seed))) {
			addSeed <- function(i, inits) {
				newInits <- inits[[i]]
				newInits[".RNG.name"] <- "base::Wichmann-Hill"
				newInits[".RNG.seed"] <- seed + (nchains - (nchains - i)) * (nchains * K)
				newInits
			}
			inits <- lapply(1:nchains, addSeed, inits = inits)
		}
		if (do.pb) {
			model <- rjags::jags.model(file = textConnection(object@model), data = list(x = object@x, N = length(object@x)), inits = inits, n.chains = nchains, quiet = FALSE)
			cat("\nUpdating\n \n")
			if (object@burnin > 0) update(model, object@burnin, progress.bar = "text")
			cat("\nSampling\n \n")
			mcmcpara <- rjags::coda.samples(model, c("mu", "log_sd", "inv_df"), n.iter = K * object@thin, thin = object@thin)
		}
		else {
			model <- rjags::jags.model(file = textConnection(object@model), data = list(x = object@x, N = length(object@x)), inits = inits, n.chains = nchains, quiet = TRUE)
			if (object@burnin > 0) update(model, object@burnin, progress.bar = "none")
			mcmcpara <- rjags::coda.samples(model, c("mu", "log_sd", "inv_df"), n.iter = K * object@thin, thin = object@thin, progress.bar = "none")
		}
		listPara <- lapply(mcmcpara, function(l) mapply(function(mu, sd, df) methods::new("studentDist", mu = mu, sd = sd, df = df), mu = l[, "mu"], sd = exp(l[, "log_sd"]), df = 1 / l[, "inv_df"]))
		listPara
	}
)

#' @rdname samplePosteriorDist
#' @aliases samplePosteriorDist.fitLogisticDist
setMethod("samplePosteriorDist", "fitLogisticDist",
	function(object, K, do.pb, seed) {
		nchains <- length(object@ini)
		convertIni <- function(l) {
			l$scale <- log(l$scale)
			names(l)[names(l) == "scale"] <- "log_scale"
			l
		}
		inits <- lapply(object@ini, convertIni)
		if (!(is.null(seed))) {
			addSeed <- function(i, inits) {
				newInits <- inits[[i]]
				newInits[".RNG.name"] <- "base::Wichmann-Hill"
				newInits[".RNG.seed"] <- seed + (nchains - (nchains - i)) * (nchains * K)
				newInits
			}
			inits <- lapply(1:nchains, addSeed, inits = inits)
		}
		if (do.pb) {
			model <- rjags::jags.model(file = textConnection(object@model), data = list(x = object@x, N = length(object@x)), inits = inits, n.chains = nchains, quiet = FALSE)
			cat("\nUpdating\n \n")
			if (object@burnin > 0) update(model, object@burnin, progress.bar = "text")
			cat("\nSampling\n \n")
			mcmcpara <- rjags::coda.samples(model, c("location", "log_scale"), n.iter = K * object@thin, thin = object@thin)
		}
		else {
			model <- rjags::jags.model(file = textConnection(object@model), data = list(x = object@x, N = length(object@x)), inits = inits, n.chains = nchains, quiet = TRUE)
			if (object@burnin > 0) update(model, object@burnin, progress.bar = "none")
			mcmcpara <- rjags::coda.samples(model, c("location", "log_scale"), n.iter = K * object@thin, thin = object@thin, progress.bar = "none")
		}
		listPara <- lapply(mcmcpara, function(l) mapply(function(location, scale) methods::new("logisticDist", location = location, scale = scale), location = l[, "location"], scale = exp(l[, "log_scale"])))
		listPara
	}
)
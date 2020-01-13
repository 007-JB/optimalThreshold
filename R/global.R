#################################################################################################
###                                  global.R                                                 ###
#################################################################################################
#' Specify which distribution to fit on the marker values
#'
#' This function is a wrapper to create an S4 object to specify a distribution to fit the marker values.
#'
#' @param x a vector of marker values (NA values allowed, see Details).
#' @param distr a character that specifies the distribution to fit (normal, log-normal, scaled t, gamma, logistic, user-defined or undefined, see Details).
#' @param ini specification of initial values for the parameters of the marker distribution in the form of a list. Each list must be named. A list should be provided for each MCMC chain. NULL for "norm" and "lnorm".
#' @param thin the thinning interval between consecutive observations. NULL for "norm" and "lnorm".
#' @param burnin a positive integer that defines the length of the burn-in iterations when performing the MCMC algorithm. NULL for "norm" and "lnorm".
#' @param model a character string used to define the model. Must match with the definition of a model compatible with JAGS. Necessary only for the t and logistic distributions (see Details).
#' @param paraNames a string vector containing the names of the parameters of the submitted distribution. Should be provided only for "user" defined distribution.
#' @param mcmcList an object of class mcmc.list where each list contains an MCMC chain. To be provided only for "user" defined distribution.
#' @param cdf a function that characterizes the cumulative distribution. To be provided only for "user" defined distribution (see Details).
#' @param gradient a function that characterizes the density distribution. To be provided only for "user" defined distribution (see Details).
#' @param hessian a function that characterizes the first derivative of the probability density function. To be provided only for "user" defined distribution (see Details).
#' @return Returns an object to be passed to the \code{trtSelThresh} and \code{diagThresh} functions.
#' @details This function allows the user to specify which distribution should be fitted to the marker values. If NA values are present in the \code{x} argument passed to the function, a warning is produced. However, the user should not discard the NA values from the original data because the length of the \code{x} argument is calculated internally to to estimate the mean risk of event occurrence in each treatment arm. So NA values are managed internally by the function.
#'Five theoretical distributions are implemented by the package: normal, log-normal, gamma, scaled t, and logistic. This is here that the user must specify which of the four distributions must be of type 'undefined' (or in other words which distribution must be expressed as a function of the three other distributions and mean risks of event). The user may also define its own theoretical distribution. The details for each theoretical distribution are provided hereafter:
#' \itemize{
#' 	\item Fit a normal distribution: when specifying \code{distr="norm"} you fit a normal distribution to the marker values passed to the \code{x} argument of the function. Non-informative priors are used (\eqn{p(\mu,\sigma^2) \propto (\sigma^2)^(-1)}). Posterior values of the normal distribution parameters are sampled directly from the exact posterior distributions. If you don't want to use non-informative priors, see the explanation on how to fit a user-defined distribution.
#'	\item Fit a log-normal distribution: when specifying \code{distr="lnorm"} you fit a log-normal distribution to the marker values passed to the \code{x} argument of the function. Non-informative priors are used (\eqn{p(\mu,\sigma^2) \propto (\sigma^2)^(-1)}). Posterior values of the log-normal distribution parameters are sampled directly from the exact posterior distributions. If you don't want to use non-informative priors, see the explanation on how to fit a user-defined distribution.
#'	\item Fit a gamma distribution: when specifying \code{distr="gamma"} you fit a gamma distribution to the marker values passed to the \code{x} argument of the function. Non-informative priors are used (\eqn{p(shape,scale) \propto 1/scale}). Posterior values of the gamma distribution parameters are sampled using the ARS method. This method requires that the user specifies a list of initial values passed to the \code{ini} argument of the function. Each element of this list must be a list with one element named 'shape'. It also requires the \code{thin} of the MCMC chain, and the length of the burnin phase passed to the \code{burnin} argument. If you don't want to use non-informative priors, see the explanation on how to fit a user-defined distribution.
#'	\item Fit a scaled t distribution: when specifying \code{distr="t"} you fit a scaled t distribution to the marker values passed to the \code{x} argument of the function. Posterior values of the scaled t distribution parameters are sampled using an MCMC algorithm through the JAGS software, so the function requires the user to provide the JAGS model as a character string through the \code{model} argument of the function. If \code{NULL}, a model with vague priors is provided to the function automatically:
#'	\deqn{mu ~ U(min(x),max(x))}
#'	\deqn{log(sd) ~ U(-10,10)}
#'	\deqn{1/df ~ U(0,1)}
#' This method requires that the user specifies a list of initial values passed to the \code{ini} argument of the function. Each element of this list must be a list with three elements named 'mu', 'sd', and 'df'. It also requires the \code{thin} of the MCMC chain, and the length of the burnin phase passed to the \code{burnin} argument. 
#'	\item Fit a logistic distribution: when specifying \code{distr="logis"} you fit a logistic distribution to the marker values passed to the \code{x} argument of the function. Posterior values of the logistic distribution parameters are sampled using a MCMC algorithm through the JAGS software, so the function requires the user to provide the JAGS model as a character string through the \code{model} argument of the function. If \code{NULL}, a model with vague priors is provided to the function automatically:
#'	\deqn{location ~ U(min(x),max(x))}
#'	\deqn{log(scale) ~ U(-10,10)}
#' This method requires that the user specifies a list of initial values passed to the \code{ini} argument of the function. Each element of this list must be a list with two elements named 'location', and 'scale'. It also requires the \code{thin} of the MCMC chain, and the length of the burnin phase passed to the \code{burnin} argument. 
#'	\item Fit a user-defined distribution: when specifying \code{distr="user"} you fit a user-defined distribution to the marker values passed to the \code{x} argument of the function. First of all, the user must give the parameters name in the argument \code{paraNames} of the function using a character vector. Then, the user provides a posterior sample of the parameters of the distribution obtained using JAGS or another software through an object of class \code{mcmc.list} to the argument \code{mcmcList} of the function (this implies that the user performed the Bayesian inference himself). Note that the names passed to the \code{mcmc.list} object must match with the names given in the \code{paraNames} argument. Then, the user must specify the \code{cdf}, \code{gradient}, and \code{hessian} functions associated with the fitted distribution. The \code{cdf} function is the cumulative distribution function that is fitted to the marker values, the \code{gradient} function is its first derivative which corresponds to the probability density function fitted to the marker values, and the \code{hessian} function is the second derivative of \code{cdf}. When the fitted distribution is a supported distribution (e.g. a normal distribution with informative priors), the user may use the \code{getMethod(cdf,"normalDist")} function to use the standard method for normal distribution used in the package. When the fitted distribution is not supported, the user must specify directly the \code{cdf} as \code{function(x,mu,sd) pnorm(x,mu,sd)} (if we keep the example of the normal distribution). The same idea may be used for the \code{gradient} and \code{hessian} functions (see the examples to have more details).
#'	\item Specify which marker distribution is expressed as a function of the three others and the mean risks of event using \code{distr="undefined"}.
#' }
#' @seealso \code{\link[optimalThreshold]{trtSelThresh}} and \code{\link[optimalThreshold]{diagThresh}}.
#' @examples
#' #Fit a normal distribution
#' x <- rnorm(250)
#' fitX <- fit(x, "norm")
#'
#' #Fit a log-normal distribution
#' x <- rlnorm(250)
#' fitX <- fit(x, "lnorm")
#'
#' #Fit a gamma distribution
#' x <- rgamma(250, shape = 2, scale = 1.2)
#' fitX <- fit(x, "gamma", 
#'             ini = list(list(shape = 1), 
#'                        list(shape = 2), 
#'                        list(shape = 3)),
#'             thin = 1, burnin = 1000)
#'
#' #Fit a scaled t distribution
#' x <- optimalThreshold:::rt.scaled(250, df = 4, mean = 2.5, sd = 2)
#' fitX <- fit(x, "t",
#'             ini = list(list(mu = 1, sd = 1, df = 2), 
#'                        list(mu = 2, sd = 2, df = 4), 
#'                        list(mu = 3, sd = 3, df = 6)),
#'             thin = 1, burnin = 1000, model = NULL)
#' 
#' #Fit a logistic distribution
#' x <- rlogis(250)
#' fitX <- fit(x, "logis", 
#'             ini = list(list(location = 0.3, scale = 0.5), 
#'                        list(location = 1, scale = 1), 
#'                        list(location = 2, scale = 2)), 
#'             thin = 1, burnin = 1000, model = NULL)
#'
#' #Specify which distribution is 'undefined'
#' x <- rnorm(250)
#' fitX <- fit(x, "undefined")
#'
#' #Fit a user-defined normal distribution with informative priors
#' library(rjags)
#' x <- rnorm(250, mean = 2, sd = 1)
#' model <- "model
#'			{
#'				mu ~ dunif(0, 4)
#'				log_sd ~ dunif(-1, 1)
#'				sd <- exp(log_sd)
#'				tau <- 1 / (sd^2)
#'				for (i in 1:N)
#'				{
#'					x[i] ~ dnorm(mu, tau)
#'				}
#'			}
#'			"
#' modelJAGS <- jags.model(file = textConnection(model), data = list(x = x, N = length(x)), 
#'                         inits = list(list(mu = 1, log_sd = -0.5),list(mu = 3.5, log_sd = 0.5)),
#'                         n.chains = 2, quiet = TRUE)
#' update(modelJAGS, 1000, progress.bar = "text")
#' mcmcpara <- coda.samples(modelJAGS, c("mu", "log_sd"), n.iter = 2000, thin = 1)
#' varnames(mcmcpara) <- c("mu", "sd")
#' mcmcpara[[1]][, "sd"] <- exp(mcmcpara[[1]][, "sd"])
#' mcmcpara[[2]][, "sd"] <- exp(mcmcpara[[2]][, "sd"])
#' fitX <- fit(x, "user", paraNames = varnames(mcmcpara), mcmcList = mcmcpara, 
#'             cdf = function(x, mu, sd) pnorm(x, mu, sd), 
#'             gradient = getMethod(gradient, "normalDist"), 
#'             hessian = function(x, mu, sd) ((mu - x) / sd^2) * dnorm(x, mu, sd))
#'
#' @export fit
fit <- function(x, distr, ini = NULL, thin = NULL, burnin = NULL, model = NULL, paraNames = NULL, mcmcList = NULL, cdf = NULL, gradient = NULL, hessian = NULL) {
	if (length(distr) != 1) stop("You must supply one 'distr' argument.")
	if (!(distr %in% c("norm", "lnorm", "gamma", "t", "logis", "user", "undefined"))) stop("You must specify a supported distribution type.")
	n <- length(x)
	if (any(is.na(x))) {
		x <- x[!(is.na(x))]
		warning("NA values were discarded from 'x' argument.")
	}
	if (distr == "norm") { 
		if (!(is.null(ini))) warning("'ini' argument is ignored as a 'norm' distribution does not need initial values.")
		if (!(is.null(thin))) warning("'thin' argument is ignored as a 'norm' distribution does not need 'thin' to be specified.")
		if (!(is.null(burnin))) warning("'burnin' argument is ignored as a 'norm' distribution does not need 'burnin' to be specified.")
		if (!(is.null(model))) warning("'model' argument is ignored as a 'norm' distribution does not need 'model' to be specified.")
		if (!(is.null(paraNames))) warning("'paraNames' argument is ignored as a 'norm' distribution does not need 'paraNames' to be specified.")
		if (!(is.null(cdf))) warning("'cdf' argument is ignored as a 'norm' distribution does not need 'cdf' to be specified.")
		if (!(is.null(gradient))) warning("'gradient' argument is ignored as a 'norm' distribution does not need 'gradient' to be specified.")
		if (!(is.null(hessian))) warning("'hessian' argument is ignored as a 'norm' distribution does not need 'hessian' to be specified.")
		return(methods::new("fitNormalDist", x = x, n = n, mcmc = FALSE))
	}
	if (distr == "lnorm") {	
		if (!(is.null(ini))) warning("'ini' argument is ignored as a 'lnorm' distribution does not need initial values.")
		if (!(is.null(thin))) warning("'thin' argument is ignored as a 'lnorm' distribution does not need 'thin' to be specified.")
		if (!(is.null(burnin))) warning("'burnin' argument is ignored as a 'lnorm' distribution does not need 'burnin' to be specified.")
		if (!(is.null(model))) warning("'model' argument is ignored as a 'lnorm' distribution does not need 'model' to be specified.")
		if (!(is.null(paraNames))) warning("'paraNames' argument is ignored as a 'lnorm' distribution does not need 'paraNames' to be specified.")
		if (!(is.null(cdf))) warning("'cdf' argument is ignored as a 'lnorm' distribution does not need 'cdf' to be specified.")
		if (!(is.null(gradient))) warning("'gradient' argument is ignored as a 'lnorm' distribution does not need 'gradient' to be specified.")
		if (!(is.null(hessian))) warning("'hessian' argument is ignored as a 'lnorm' distribution does not need 'hessian' to be specified.")
		return(methods::new("fitLogNormalDist", x = x, n = n, mcmc = FALSE))
	}
	if (distr == "gamma") {
		if (!(is.null(model))) warning("'model' argument is ignored as a 'gamma' distribution does not need 'model' to be specified.")
		if (!(is.null(paraNames))) warning("'paraNames' argument is ignored as a 'gamma' distribution does not need 'paraNames' to be specified.")
		if (!(is.null(cdf))) warning("'cdf' argument is ignored as a 'gamma' distribution does not need 'cdf' to be specified.")
		if (!(is.null(gradient))) warning("'gradient' argument is ignored as a 'gamma' distribution does not need 'gradient' to be specified.")
		if (!(is.null(hessian))) warning("'hessian' argument is ignored as a 'gamma' distribution does not need 'hessian' to be specified.")
		return(methods::new("fitGammaDist", x = x, n = n, ini = ini, thin = as.integer(thin), burnin = as.integer(burnin), mcmc = TRUE))
	}
	if (distr == "t") {
		if (!(is.null(paraNames))) warning("'paraNames' argument is ignored as a 't' distribution does not need 'paraNames' to be specified.")
		if (!(is.null(cdf))) warning("'cdf' argument is ignored as a 't' distribution does not need 'cdf' to be specified.")
		if (!(is.null(gradient))) warning("'gradient' argument is ignored as a 't' distribution does not need 'gradient' to be specified.")
		if (!(is.null(hessian))) warning("'hessian' argument is ignored as a 't' distribution does not need 'hessian' to be specified.")
		if (is.null(model)) {
		rangex <- range(x)
		model = paste("model
			{
				mu ~ dunif(",rangex[1],", ",rangex[2],")
				log_sd ~ dunif(-10, 10)
				sd <- exp(log_sd)
				tau <- 1 / (sd ^ 2)
				inv_df ~ dunif(0, 1)
				df <- 1 / inv_df
		
				for (i in 1:N)
				{
					x[i] ~ dt(mu, tau, df)
				}
			}
			",sep="")
		}
		return(methods::new("fitStudentDist", x = x, n = n, ini = ini, thin = as.integer(thin), burnin = as.integer(burnin), model = model, mcmc = TRUE))
	}
	if (distr == "logis") {
		if (!(is.null(paraNames))) warning("'paraNames' argument is ignored as a 'logis' distribution does not need 'paraNames' to be specified.")
		if (!(is.null(cdf))) warning("'cdf' argument is ignored as a 'logis' distribution does not need 'cdf' to be specified.")
		if (!(is.null(gradient))) warning("'gradient' argument is ignored as a 'logis' distribution does not need 'gradient' to be specified.")
		if (!(is.null(hessian))) warning("'hessian' argument is ignored as a 'logis' distribution does not need 'hessian' to be specified.")
		if (is.null(model)) {
		rangex <- range(x)
		model = paste("model
			{
				location ~ dunif(",rangex[1],", ",rangex[2],")
				log_scale ~ dunif(-10, 10)
				scale <- exp(log_scale)
				tau <- 1 / scale
				
				for (i in 1:N)
				{
					x[i] ~ dlogis(location, tau)
				}
			}
			",sep="")
		}
		return(methods::new("fitLogisticDist", x = x, n = n, ini = ini, thin = as.integer(thin), burnin = as.integer(burnin), model = model, mcmc = TRUE))
	}
	if (distr == "user") {
		if (!(is.null(ini))) warning("'ini' argument is ignored as a 'user' distribution does not need initial values.")
		if (!(is.null(thin))) warning("'thin' argument is ignored as a 'user' distribution does not need 'thin' to be specified.")
		if (!(is.null(burnin))) warning("'burnin' argument is ignored as a 'user' distribution does not need 'burnin' to be specified.")
		if (!(is.null(model))) warning("'model' argument is ignored as a 'user' distribution does not need 'model' to be specified.")
		return(methods::new("fitUserDefinedDist", x = x, n = n, paraNames = paraNames, mcmcList = mcmcList, cdf = cdf, gradient = gradient, hessian = hessian, mcmc = FALSE))
	}
	if (distr == "undefined") {
		return(methods::new("undefined", x = x, n = n, mcmc = FALSE))
	}
}

compoundEvtRefDist <- function(NoEvtRefDist, EvtInnovDist, NoEvtInnovDist, r0, r1) {
	methods::new("compoundEvtRefDist", NoEvtRefDist = NoEvtRefDist, EvtInnovDist = EvtInnovDist, NoEvtInnovDist = NoEvtInnovDist, r0 = r0, r1 = r1)
}
compoundNoEvtRefDist <- function(EvtRefDist, EvtInnovDist, NoEvtInnovDist, r0, r1) {
	methods::new("compoundNoEvtRefDist", EvtRefDist = EvtRefDist, EvtInnovDist = EvtInnovDist, NoEvtInnovDist = NoEvtInnovDist, r0 = r0, r1 = r1)
}
compoundEvtInnovDist <- function(EvtRefDist, NoEvtRefDist, NoEvtInnovDist, r0, r1) {
	methods::new("compoundEvtInnovDist", EvtRefDist = EvtRefDist, NoEvtRefDist = NoEvtRefDist, NoEvtInnovDist = NoEvtInnovDist, r0 = r0, r1 = r1)
}
compoundNoEvtInnovDist <- function(EvtRefDist, NoEvtRefDist, EvtInnovDist, r0, r1) {
	methods::new("compoundNoEvtInnovDist", EvtRefDist = EvtRefDist, NoEvtRefDist = NoEvtRefDist, EvtInnovDist = EvtInnovDist, r0 = r0, r1 = r1)
}
utility <- function(x, r0, r1, r, pT0, pT1, coefU, EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist) {
		U <- (- coefU * (cdf(EvtInnovDist)(x) * r1 - cdf(EvtRefDist)(x) * r0 + r * (pT0 * (cdf(EvtRefDist)(x) * r0 + cdf(NoEvtRefDist)(x) * (1 - r0)) + pT1 * (cdf(EvtInnovDist)(x) * r1 + cdf(NoEvtInnovDist)(x) * (1 - r1)))))
		attr(U, 'gradient') <- (- coefU * (gradient(EvtInnovDist)(x) * r1 - gradient(EvtRefDist)(x) * r0 + r * (pT0 * (gradient(EvtRefDist)(x) * r0 + gradient(NoEvtRefDist)(x) * (1 - r0)) + pT1 * (gradient(EvtInnovDist)(x) * r1 + gradient(NoEvtInnovDist)(x) * (1 - r1)))))
		attr(U, 'hessian') <- (- coefU * (hessian(EvtInnovDist)(x) * r1 - hessian(EvtRefDist)(x) * r0 + r * (pT0 * (hessian(EvtRefDist)(x) * r0 + hessian(NoEvtRefDist)(x) * (1 - r0)) + pT1 * (hessian(EvtInnovDist)(x) * r1 + hessian(NoEvtInnovDist)(x) * (1 - r1)))))
		U
}
defineStartVal <- function(gridStartVal, utility, hessTol, ...) {
	Us <- utility(x = gridStartVal, ...)
	if (all(diff(Us) <= 0) | all(diff(Us) >= 0)) {
		startVal <- NA
		codeVal <- "Monotonic function"
		return(list(startVal, codeVal))
	}
	UsHessPos <- Us[abs(attr(Us, "hessian")) > hessTol & !(is.na(Us))]
	if (length(UsHessPos) == 0) {
		startVal <- NA
		codeVal <- "No minimum found"
		return(list(startVal, codeVal))
	}
	startVal <- gridStartVal[which(Us == min(Us[abs(attr(Us, "hessian")) > hessTol], na.rm = TRUE))][1]
	codeVal <- "Convex function"
	return(list(startVal, codeVal))
}
explicitThreshold <- function(EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist, r0, r1, r) {
	if ((methods::is(EvtRefDist, "compoundEvtRefDist") | methods::is(EvtInnovDist, "compoundEvtInnovDist")) & methods::is(NoEvtRefDist, "normalDist") & methods::is(NoEvtInnovDist, "normalDist") & r == 0) {
		a <- NoEvtRefDist@mu - NoEvtInnovDist@mu
		R <- (1 - r0) / (1 - r1)
		if (NoEvtRefDist@sd != NoEvtInnovDist@sd) {
			b <- NoEvtRefDist@sd / NoEvtInnovDist@sd
			thresholdEst <- (NoEvtInnovDist@mu * (b ^ 2 - 1) - a + b * sqrt(a ^ 2 + (b ^ 2 - 1) * NoEvtInnovDist@sd ^ 2 * log(b ^ 2 * R ^ 2))) / (b ^ 2 - 1)
			return(thresholdEst)
		}
		else {
			thresholdEst <- (NoEvtInnovDist@sd^2 * log(R ^ 2) + NoEvtRefDist@mu ^ 2 - NoEvtInnovDist@mu ^ 2) / (2 * a)
			return(thresholdEst)
		}	
	}
	if ((is(NoEvtRefDist, "compoundNoEvtRefDist") | is(NoEvtInnovDist, "compoundNoEvtInnovDist")) & is(EvtRefDist, "normalDist") & is(EvtInnovDist, "normalDist") & r == 0) {
		a <- EvtInnovDist@mu - EvtRefDist@mu
		R <- r1 / r0
		if (EvtInnovDist@sd != EvtRefDist@sd) {
			b <- EvtInnovDist@sd / EvtRefDist@sd
			thresholdEst <- (EvtRefDist@mu * (b ^ 2 - 1) - a + b * sqrt(a ^ 2 + (b ^ 2 - 1) * EvtRefDist@sd ^ 2 * log(b ^ 2 * R ^ 2))) / (b ^ 2 - 1)
			return(thresholdEst)
		}
		else {
			thresholdEst <- (EvtInnovDist@sd ^ 2 * log(R ^ 2) + EvtInnovDist@mu ^ 2 - EvtRefDist@mu ^ 2) / (2 * a)
			return(thresholdEst)
		}
	}
	return(NA)
}
maxUtility <- function(i, EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist, r0, r1, pT0, pT1, r, coefU, minMarker, maxMarker, utility, hessTol, pb) {
	if (!is.null(pb)) setTxtProgressBar(pb, i)
	thresholdEst <- explicitThreshold(EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist, r0, r1, r)
	if (!(is.na(thresholdEst))) return(thresholdEst)
	gridStartVal <- seq(minMarker, maxMarker, length.out = 1000)
	startVal <- defineStartVal(gridStartVal, utility, hessTol, r0, r1, r, pT0, pT1, coefU, EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist)
	if (startVal[[2]] %in% c("Monotonic function", "No minimum found")) {
		warning("Utility function is not convex.")
		return(NA)
	}
	thresholdEst <- stats::nlm(utility, startVal[[1]], r0 = r0, r1 = r1, r = r, pT0 = pT0, pT1 = pT1, coefU = coefU, EvtRefDist = EvtRefDist, NoEvtRefDist = NoEvtRefDist, EvtInnovDist = EvtInnovDist, NoEvtInnovDist = NoEvtInnovDist, hessian = FALSE, check.analyticals = TRUE)
	Uopt <- utility(thresholdEst$estimate, r0, r1, r, pT0, pT1, coefU, EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist)
	if (is.nan(attr(Uopt, "hessian"))) {
		warning("Hessian of the utility function is not defined at threshold estimate.")
		return(NA)
	}
	if (abs(attr(Uopt, "hessian")) < hessTol) {
		warning("Hessian of the utility function at threshold estimate is too close to zero to correspond to a minimum estimate.")
		return(NA)
	}
	return(thresholdEst$estimate)
}
markerBasedRisk <- function(EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist, r0, r1, Threshold, lowRef) {
	if (lowRef) {
		markerBasedRiskRef <- cdf(EvtRefDist)(Threshold) * r0 / (cdf(EvtRefDist)(Threshold) * r0 + cdf(NoEvtRefDist)(Threshold) * (1 - r0))
		markerBasedRiskInnov <- (1 - cdf(EvtInnovDist)(Threshold)) * r1 / ((1 - cdf(EvtInnovDist)(Threshold)) * r1 + (1 - cdf(NoEvtInnovDist)(Threshold)) * (1 - r1))
	}
	else {
		markerBasedRiskRef <- (1 - cdf(EvtRefDist)(Threshold)) * r0 / ((1 - cdf(EvtRefDist)(Threshold)) * r0 + (1 - cdf(NoEvtRefDist)(Threshold)) * (1 - r0))
		markerBasedRiskInnov <- cdf(EvtInnovDist)(Threshold) * r1 / (cdf(EvtInnovDist)(Threshold) * r1 + cdf(NoEvtInnovDist)(Threshold) * (1 - r1))
	}
	c(markerBasedRiskRef = markerBasedRiskRef, markerBasedRiskInnov = markerBasedRiskInnov)
}
defineUserDist0E <- function(EvtRefDist) {
	.optimalThreshold_defineUser_Env <- new.env()
	if (!(is(EvtRefDist@cdf, "MethodDefinition"))) {
		fCdf <- EvtRefDist@cdf
		formals(fCdf) <- alist(x =)
		EvtRefDist@cdf <- eval(parse(text = c("function(object) {", paste(paste(EvtRefDist@paraNames, EvtRefDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fCdf), "}")))
	}
	if (!(is(EvtRefDist@gradient, "MethodDefinition"))) {
		fGradient <- EvtRefDist@gradient
		formals(fGradient) <- alist(x =)
		EvtRefDist@gradient <- eval(parse(text = c("function(object) {",paste(paste(EvtRefDist@paraNames, EvtRefDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fGradient), "}")))
	}
	if (!(is(EvtRefDist@hessian, "MethodDefinition"))) {
		fHessian <- EvtRefDist@hessian
		formals(fHessian) <- alist(x =)
		EvtRefDist@hessian <- eval(parse(text = c("function(object) {",paste(paste(EvtRefDist@paraNames, EvtRefDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fHessian), "}")))
	}
	paraTypes <- rep("numeric", length(EvtRefDist@paraNames))
	names(paraTypes) <- EvtRefDist@paraNames
	methods::setClass("userDefinedDist0E", slots = paraTypes, where = .optimalThreshold_defineUser_Env)
	methods::setClass("fitUserDefinedDist0E", slots = c(x = "numeric", n = "numeric", paraNames = "character", mcmcList = "mcmc.list", mcmc = "logical"), where = .optimalThreshold_defineUser_Env)
	methods::setMethod("samplePosteriorDist", "fitUserDefinedDist0E",
		function(object, K, do.pb, seed) {
			listPara <- eval(parse(text = paste("lapply(object@mcmcList, function(mcmcList) mapply(function(",paste(object@paraNames, collapse = ","), ") methods::new('userDefinedDist0E', ",paste(paste(object@paraNames, object@paraNames, sep=" = "), collapse = ","), "),", paste(paste(object@paraNames, paste("as.vector(mcmcList[,'", object@paraNames, "'])", sep = ""), sep = " = "), collapse = ","), "))", sep = "")))
			listPara
		},
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("cdf", "userDefinedDist0E",
		EvtRefDist@cdf,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("gradient", "userDefinedDist0E",
		EvtRefDist@gradient,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("hessian", "userDefinedDist0E",
		EvtRefDist@hessian,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setIs("userDefinedDist0E", "allowedDist", where = .optimalThreshold_defineUser_Env)
	list(env = .optimalThreshold_defineUser_Env, object = methods::new("fitUserDefinedDist0E", x = EvtRefDist@x, n = EvtRefDist@n, paraNames = EvtRefDist@paraNames, mcmcList = EvtRefDist@mcmcList, mcmc = EvtRefDist@mcmc))
}
defineUserDist0Eb <- function(NoEvtRefDist) {
	.optimalThreshold_defineUser_Env <- new.env()
	if (!(is(NoEvtRefDist@cdf, "MethodDefinition"))) {
		fCdf <- NoEvtRefDist@cdf
		formals(fCdf) <- alist(x =)
		NoEvtRefDist@cdf <- eval(parse(text = c("function(object) {",paste(paste(NoEvtRefDist@paraNames, NoEvtRefDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fCdf), "}")))
	}
	if (!(is(NoEvtRefDist@gradient, "MethodDefinition"))) {
		fGradient <- NoEvtRefDist@gradient
		formals(fGradient) <- alist(x =)
		NoEvtRefDist@gradient <- eval(parse(text = c("function(object) {",paste(paste(NoEvtRefDist@paraNames, NoEvtRefDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fGradient), "}")))
	}
	if (!(is(NoEvtRefDist@hessian, "MethodDefinition"))) {
		fHessian <- NoEvtRefDist@hessian
		formals(fHessian) <- alist(x =)
		NoEvtRefDist@hessian <- eval(parse(text = c("function(object) {",paste(paste(NoEvtRefDist@paraNames, NoEvtRefDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fHessian), "}")))
	}
	paraTypes <- rep("numeric", length(NoEvtRefDist@paraNames))
	names(paraTypes) <- NoEvtRefDist@paraNames
	methods::setClass("userDefinedDist0Eb", slots = paraTypes, where = .optimalThreshold_defineUser_Env)
	methods::setClass("fitUserDefinedDist0Eb", slots = c(x = "numeric", n = "numeric", paraNames = "character", mcmcList = "mcmc.list", mcmc = "logical"), where = .optimalThreshold_defineUser_Env)
	methods::setMethod("samplePosteriorDist", "fitUserDefinedDist0Eb",
		function(object, K, do.pb, seed) {
			listPara <- eval(parse(text = paste("lapply(object@mcmcList, function(mcmcList) mapply(function(",paste(object@paraNames, collapse = ","), ") methods::new('userDefinedDist0Eb',", paste(paste(object@paraNames, object@paraNames, sep=" = "), collapse = ","), "),", paste(paste(object@paraNames, paste("as.vector(mcmcList[,'", object@paraNames, "'])", sep = ""), sep = " = "), collapse = ","), "))", sep = "")))
			listPara
		},
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("cdf", "userDefinedDist0Eb",
		NoEvtRefDist@cdf,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("gradient", "userDefinedDist0Eb",
		NoEvtRefDist@gradient,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("hessian", "userDefinedDist0Eb",
		NoEvtRefDist@hessian,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setIs("userDefinedDist0Eb", "allowedDist", where = .optimalThreshold_defineUser_Env)
	list(env = .optimalThreshold_defineUser_Env, object = methods::new("fitUserDefinedDist0Eb", x = NoEvtRefDist@x, n = NoEvtRefDist@n, paraNames = NoEvtRefDist@paraNames, mcmcList = NoEvtRefDist@mcmcList, mcmc = NoEvtRefDist@mcmc))
}
defineUserDist1E <- function(EvtInnovDist) {
	.optimalThreshold_defineUser_Env <- new.env()
	if (!(is(EvtInnovDist@cdf, "MethodDefinition"))) {
		fCdf <- EvtInnovDist@cdf
		formals(fCdf) <- alist(x =)
		EvtInnovDist@cdf <- eval(parse(text = c("function(object) {",paste(paste(EvtInnovDist@paraNames, EvtInnovDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fCdf), "}")))
	}
	if (!(is(EvtInnovDist@gradient, "MethodDefinition"))) {
		fGradient <- EvtInnovDist@gradient
		formals(fGradient) <- alist(x =)
		EvtInnovDist@gradient <- eval(parse(text = c("function(object) {",paste(paste(EvtInnovDist@paraNames, EvtInnovDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fGradient), "}")))
	}
	if (!(is(EvtInnovDist@hessian, "MethodDefinition"))) {
		fHessian <- EvtInnovDist@hessian
		formals(fHessian) <- alist(x =)
		EvtInnovDist@hessian <- eval(parse(text = c("function(object) {",paste(paste(EvtInnovDist@paraNames, EvtInnovDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fHessian), "}")))
	}
	paraTypes <- rep("numeric", length(EvtInnovDist@paraNames))
	names(paraTypes) <- EvtInnovDist@paraNames
	methods::setClass("userDefinedDist1E", slots = paraTypes, where = .optimalThreshold_defineUser_Env)
	methods::setClass("fitUserDefinedDist1E", slots = c(x = "numeric", n = "numeric", paraNames = "character", mcmcList = "mcmc.list", mcmc = "logical"), where = .optimalThreshold_defineUser_Env)
	methods::setMethod("samplePosteriorDist", "fitUserDefinedDist1E",
		function(object, K, do.pb, seed) {
			listPara <- eval(parse(text = paste("lapply(object@mcmcList,function(mcmcList) mapply(function(", paste(object@paraNames, collapse = ","), ") methods::new('userDefinedDist1E',", paste(paste(object@paraNames, object@paraNames, sep = " = "), collapse = ","), "),", paste(paste(object@paraNames, paste("as.vector(mcmcList[,'", object@paraNames, "'])", sep = ""), sep = " = "), collapse = ","), "))", sep = "")))
			listPara
		},
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("cdf", "userDefinedDist1E",
		EvtInnovDist@cdf,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("gradient", "userDefinedDist1E",
		EvtInnovDist@gradient,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("hessian", "userDefinedDist1E",
		EvtInnovDist@hessian,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setIs("userDefinedDist1E", "allowedDist", where = .optimalThreshold_defineUser_Env)
	list(env = .optimalThreshold_defineUser_Env, object = methods::new("fitUserDefinedDist1E", x = EvtInnovDist@x, n = EvtInnovDist@n, paraNames = EvtInnovDist@paraNames, mcmcList = EvtInnovDist@mcmcList, mcmc = EvtInnovDist@mcmc))
}
defineUserDist1Eb <- function(NoEvtInnovDist) {
	.optimalThreshold_defineUser_Env <- new.env()
	if (!(is(NoEvtInnovDist@cdf, "MethodDefinition"))) {
		fCdf <- NoEvtInnovDist@cdf
		formals(fCdf) <- alist(x =)
		NoEvtInnovDist@cdf <- eval(parse(text = c("function(object) {",paste(paste(NoEvtInnovDist@paraNames, NoEvtInnovDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fCdf), "}")))
	}
	if (!(is(NoEvtInnovDist@gradient, "MethodDefinition"))) {
		fGradient <- NoEvtInnovDist@gradient
		formals(fGradient) <- alist(x =)
		NoEvtInnovDist@gradient <- eval(parse(text = c("function(object) {",paste(paste(NoEvtInnovDist@paraNames, NoEvtInnovDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fGradient), "}")))
	}
	if (!(is(NoEvtInnovDist@hessian, "MethodDefinition"))) {
		fHessian <- NoEvtInnovDist@hessian
		formals(fHessian) <- alist(x =)
		NoEvtInnovDist@hessian <- eval(parse(text = c("function(object) {",paste(paste(NoEvtInnovDist@paraNames, NoEvtInnovDist@paraNames, sep = '<-object@'), collapse = ';'), ";", deparse(fHessian), "}")))
	}
	paraTypes <- rep("numeric", length(NoEvtInnovDist@paraNames))
	names(paraTypes) <- NoEvtInnovDist@paraNames
	methods::setClass("userDefinedDist1Eb", slots = paraTypes, where = .optimalThreshold_defineUser_Env)
	methods::setClass("fitUserDefinedDist1Eb", slots = c(x = "numeric", n = "numeric", paraNames = "character", mcmcList = "mcmc.list", mcmc = "logical"), where = .optimalThreshold_defineUser_Env)
	methods::setMethod("samplePosteriorDist", "fitUserDefinedDist1Eb",
		function(object, K, do.pb, seed) {
			listPara <- eval(parse(text = paste("lapply(object@mcmcList,function(mcmcList) mapply(function(",paste(object@paraNames, collapse = ","), ") methods::new('userDefinedDist1Eb',", paste(paste(object@paraNames, object@paraNames, sep = " = "), collapse = ","), "),", paste(paste(object@paraNames, paste("as.vector(mcmcList[,'", object@paraNames, "'])", sep = ""), sep = " = "), collapse = ","), "))", sep = "")))
			listPara
		},
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("cdf", "userDefinedDist1Eb",
		NoEvtInnovDist@cdf,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("gradient", "userDefinedDist1Eb",
		NoEvtInnovDist@gradient,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setMethod("hessian", "userDefinedDist1Eb",
		NoEvtInnovDist@hessian,
		where = .optimalThreshold_defineUser_Env
	)
	methods::setIs("userDefinedDist1Eb", "allowedDist", where = .optimalThreshold_defineUser_Env)
	list(env = .optimalThreshold_defineUser_Env, object = methods::new("fitUserDefinedDist1Eb", x = NoEvtInnovDist@x, n = NoEvtInnovDist@n, paraNames = NoEvtInnovDist@paraNames, mcmcList = NoEvtInnovDist@mcmcList, mcmc = NoEvtInnovDist@mcmc))
}
Mode <- function(x, ...) {
  densX <- stats::density(x, ...)
  return(densX$x[which.max(densX$y)])
}

#' Marker-by-treatment predictiveness curves plot
#'
#' This function plots the marker-by-treatment predictiveness curves for treatment selection markers, corresponding to the risk of event in each treatment arm in function of the cumulative distribution of the marker.
#'
#' @param x0E a numeric vector of the marker values for patients in the reference arm that developed the event.
#' @param x0Eb a numeric vector of the marker values for patients in the reference arm that did not develop the event.
#' @param x1E a numeric vector of the marker values for patients in the innovative arm that developed the event.
#' @param x1Eb a numeric vector of the marker values for patients in the innovative arm that did not develop the event.
#' @param ylab label of the Y-axis.
#' @param xlab label of the X-axis.
#' @param main title of the graph.
#' @return None
#' @details The function uses regression splines to plot the marker-by-treatment predictiveness curves. This graph may be used to check graphically the strength of the marker-by-treatment interaction, and to know whether low values of the marker are associated with a better response of the reference treatment (this information is needed in the \code{trtSelThresh} function).
#' @section References:
#' Janes, H, Pepe, MS, Bossuyt, PM, and Barlow, WE. Measuring the performance of markers for guiding treatment decisions. \emph{Annals of Internal Medicine}. 2011; 154(4): 253-259.
#' @seealso \code{\link[mgcv]{gam}} for more details about regression splines.
#' @examples
#' x0E <- rnorm(100, 2, 1)
#' x0Eb <- rnorm(100, 4, 1)
#' x1E <- rnorm(100, 4, 1)
#' x1Eb <- rnorm(100, 2, 1)
#' riskCurves(x0E, x0Eb, x1E, x1Eb)
#' @export riskCurves
riskCurves <- function(x0E, x0Eb, x1E, x1Eb, ylab = "Predicted risk of event", xlab = "Empirical cumulative distribution function of the marker", main = "Marker-by-treatment predictiveness curves") {
	opar <- graphics::par(no.readonly = TRUE)
	on.exit(graphics::par(opar))
	x <- c(x0E, x0Eb, x1E, x1Eb)
	Fx <- stats::ecdf(x)
	d0 <- data.frame(Evt = c(rep(1, length(x0E)), rep(0, length(x0Eb))), Marker = c(x0E, x0Eb))
	m0 <- mgcv::gam(Evt ~ s(Marker), family = "binomial", data = d0)
	d1 <- data.frame(Evt = c(rep(1, length(x1E)), rep(0, length(x1Eb))), Marker = c(x1E, x1Eb))
	m1 <- mgcv::gam(Evt ~ s(Marker), family = "binomial", data = d1)
	newd <- data.frame(Marker = seq(min(d0$Marker, d1$Marker, na.rm = TRUE), max(d0$Marker, d1$Marker, na.rm = TRUE), length.out = 20))
	graphics::par(lwd = 2)
	graphics::plot(mgcv::predict.gam(m0, newdata = newd, type = "response") ~ Fx(newd$Marker), type = "l", col = "blue", ylab = ylab, xlab = xlab, main = main, sub = "Blue: reference treatment; Green: innovative treatment", ylim = c(0,1), xlim = c(0,1))
	graphics::lines(mgcv::predict.gam(m1, newdata = newd, type = "response") ~ Fx(newd$Marker), col = "green")
}

#' Density curves plot
#'
#' This function plots the density curves of the marker values in each treatment arm (treatment selection marker), or for diseased and non-diseased patients (diagnostic marker).
#'
#' @param x0 a numeric vector containing the marker values of patients in the reference arm (treatment selection marker), or non-diseased patients (diagnostic marker).
#' @param x1 a numeric vector containing the marker values of patients in the innovative arm (treatment selection marker), or diseased patients (diagnostic marker).
#' @param type a character argument that specifies the type of the marker ("treatment selection" for treatment selection marker or "diagnostic" for diagnostic marker).
#' @param ylab label of the Y-axis.
#' @param xlab label of the X-axis.
#' @param main title of the graph.
#' @param col0 color of the density curve in the reference arm (treatment selection marker), or for non-diseased patients (diagnostic marker).
#' @param col1 color of the density curve in the innovative arm (treatment selection marker), or for diseased patients (diagnostic marker).
#' @param lty0 type of the line for the density curve in the reference arm (treatment selection marker), or for non-diseased patients (diagnostic marker).
#' @param lty1 type of the line for the density curve in the innovative arm (treatment selection marker), or for diseased patients (diagnostic marker).
#' @param pos.legend the x and y co-ordinates to be used to position the legend. They can be specified by keyword or in any way which is accepted by xy.coords.
#' @param ... other arguments to be passed to the plot function.
#' @return None
#' @details When assessing treatment selection markers and estimating their optimal threshold, it is necessary that the randomization constraint be respected. If the density curves of the marker are different when comparing the two treatment arms, then it is likely that the \code{trtSelThresh} function will provide a threshold that do not correspond to the true optimal threshold. 
#' When assessing diagnostic markers, it is necessary to define the decision rule (classically high values of the marker are associated with a worst outcome). This decision rule may be checked with the density curves of the marker for diseased and non-diseased patients.  
#' @examples
#' ### Plotting density curves for a treatment selection marker
#' # Data generation
#' x0E <- rnorm(100, 2, 1)
#' x0Eb <- rnorm(200, 4, 1)
#' x1E <- rnorm(100, 4, 1)
#' x1Eb <- rnorm(200, 2, 1)
#' densCurves(x0 = c(x0E, x0Eb), x1 = c(x1E, x1Eb), type = "treatment selection")
#'
#' ### Plotting density curves for a diagnostic marker
#' # Data generation
#' xE <- rnorm(30, 3, 1)
#' xEb <- rnorm(90, 1, 1)
#' densCurves(x0 = xEb, x1 = xE, type = "diagnostic")
#' @export densCurves
densCurves <- function(x0, x1, type, ylab = "Density", xlab = "Marker values", main = "Density curves", col0 = "blue", col1 = "green", lty0 = 1, lty1 = 1, pos.legend = "topright", ...) {
	opar <- graphics::par(no.readonly = TRUE)
	on.exit(graphics::par(opar))
	if (missing(type)) stop("You must specify the type of the marker under study: diagnostic or treatment selection.")
	type <- match.arg(type, c("diagnostic", "treatment selection"))
	if (type == "diagnostic") {
		legends <- c("With event", "Without event")
	}
	else {
		legends <- c("Innovative treatment", "Reference treatment")
	}
	dens0 <- stats::density(x0, na.rm = TRUE)
	dens1 <- stats::density(x1, na.rm = TRUE)
	graphics::plot(dens0$y ~ dens0$x, type = "l", ylab = ylab, xlab = xlab, main = main, col = col0, lty = lty0, ...)
	graphics::lines(dens1$y ~ dens1$x, col = col1, lty = lty1, ...)
	legend(pos.legend, col = c(col1, col0), lty = c(lty1, lty0), legend = legends)
}

#' Estimation of the optimal threshold of a treatment selection marker
#'
#' This function produces a sample of the posterior distribution of the optimal threshold of a treatment selection marker. The optimal threshold is defined as the marker value that maximized the utility of using the marker to decide between two treatment options (innovative and reference one). The utility function takes into account the efficacy of the treatment options as well as treatment induced toxicities.
#' The estimation of the utility function needs data from a clinical trial about the two treatment options, in which the success of a treatment is defined by the absence of an event of interest in a given post-treatment interval (binary data). For the time being, the package cannot estimate the optimal threshold in case of censored data about the occurrence of the event in the given post-treatment interval.
#' To calculate the utility function, the user needs to specify: 
#' \itemize{
#' \item the distribution of the marker in the four groups defined by the treatment option and the outcome; in fact only three distributions need to be specified, the fourth one being derived from the three others and the mean risks of event in the two treatment arms through the randomization constraint (the distribution of the marker being the same in both treatment arms; see the \code{\link[optimalThreshold]{fit}} function for more details),
#' \item and the mean risks of the event in the two treatment arms. The user must also specify: the cost of the innovative treatment relative to the cost of the event (see Details). 
#' }
#' The optimal threshold and its credible interval are calculated using a Monte Carlo approach. 
#'
#' @param EvtRefDist an object of class allowedFitDist that summarizes the distribution fitted to the marker values of patients that developed the event of interest in the reference arm. This class of objects is obtained using the \code{\link[optimalThreshold]{fit}} function.
#' @param NoEvtRefDist an object of class allowedFitDist that summarizes the distribution fitted to the marker values of patients that did not develop the event of interest in the reference arm. This class of objects is obtained using the \code{\link[optimalThreshold]{fit}} function.
#' @param EvtInnovDist an object of class allowedFitDist that summarizes the distribution fitted to the marker values of patients that developed the event of interest in the innovative arm. This class of objects is obtained using the \code{\link[optimalThreshold]{fit}} function.
#' @param NoEvtInnovDist an object of class allowedFitDist that summarizes the distribution fitted to the marker values of patients that did not develop the event of interest in the innovative arm. This class of objects is obtained using the \code{\link[optimalThreshold]{fit}} function.
#' @param mRiskRef an object of class mcmc.list provided by the user. It must be a sample of the posterior distribution of the mean risk of event in the reference treatment arm. If NULL, the function samples values in the posterior distribution of the mean risk of event in the reference arm assuming Jeffrey's prior (Beta(0.5,0.5)), and estimating the mean risk using the number of marker values specified in each treatment arm.
#' @param mRiskInnov an object of class mcmc.list provided by the user. It must be a sample of the posterior distribution of the mean risk of event in the innovative treatment arm. If NULL, the function samples values in the posterior distribution of the mean risk of event in the innovative arm assuming Jeffrey's prior (Beta(0.5,0.5)), and estimating the mean risk using the number of marker values specified in each treatment arm..
#' @param lowRef a logical value indicating whether low values of the marker are associated with low (TRUE) or high (FALSE) risk under the reference treatment arm.
#' @param toxRef a logical value indicating whether the reference treatment arm (TRUE) or the innovative treatment arm (FALSE) must be preferred at equal efficacy taking into account toxicity.
#' @param r a numeric value indicating the cost ratio between the most harmful treatment and the event (see Details).
#' @param le.MCMC length of the desired MCMC chain.
#' @param hessTol tolerance for the hessian value of the utility function at the optimal threshold.
#' @param plot a logical value indicating whether routine graphics must be produced.
#' @param progress.bar a character string indicating whether the user wishes to print a progress bar during the function process.
#' @param seed a numerical value used to fix the random seed.
#' @return Returns an object of class \code{trtSelOptThresh}.
#' @details When \code{toxRef==FALSE} then Janes et al. (2014) defined the costs of event and treatment as:
#' \tabular{lcc}{
#'	\tab Y=0 \tab Y=1 \cr
#' Z=0 \tab 0 \tab \eqn{C_Y} \cr
#' Z=1 \tab \eqn{C_Z} \tab \eqn{C_Z+C_Y} 
#'}
#' When \code{toxRef==TRUE} it is defined as :
#' \tabular{lcc}{
#'	\tab Y=0 \tab Y=1 \cr
#' Z=0 \tab \eqn{C_Z} \tab \eqn{C_Z+C_Y} \cr
#' Z=1 \tab 0 \tab \eqn{C_Y} 
#'}
#' According to the value of \code{toxRef}, the r ratio is simply \eqn{r=C_Z/C_Y}. The r ratio can also be indirectly specified by the absolute difference in risk of event between the two treatments above which a physician would recommend the use of the most harmful treatment. The inverse of the r ratio can also be interpreted as the number of patients for whom the physician is ready to give the most harmful treatment to prevent one additional case compared with the less harmful treatment.
#' @section References :
#' Blangero, Y, Rabilloud, M, Ecochard, R, and Subtil, F. A Bayesian method to estimate the optimal threshold of a marker used to select patients' treatment. \emph{Statistical Methods in Medical Research}. 2019.
#' @examples
#' #Simulating data from four gaussian distributions, 
#' #with mean risks equal to 0.5 in each arm:
#' x0E <- rnorm(250) # reference arm, event
#' x0Eb <- rnorm(250, 2) # reference arm, no event
#' x1E <- rnorm(250, 2) # innovative arm, event
#' x1Eb <- rnorm(250) # innovative arm, no event
#'
#' #When working with real data. You can check the randomization constraint using the 
#' #densCurves function:
#' densCurves(x0 = c(x0E, x0Eb), x1 = c(x1E, x1Eb), type = "treatment selection")
#'
#' #You can also use the riskCurves function to know if low values of the marker are associated
#' #with a better response under the reference treatment or not:
#' library(mgcv)
#' riskCurves(x0E, x0Eb, x1E, x1Eb)
#'
#' #Fit normal distributions on three groups. And let the last one (1E) be undefined (derived 
#' #indirectly using the randomization constraint):
#' fit0E <- fit(x0E, "norm")
#' fit0Eb <- fit(x0Eb, "norm")
#' fit1E <- fit(x1E, "undefined")
#' fit1Eb <- fit(x1Eb, "norm")
#'
#' #Apply the main function to estimate the optimal threshold:
#' # first case: the mean risks of event in the two treatment arms are left unspecified (are 
#' # determined by the number of marker measurements in the fit0E, fi0Eb, fit1E, fit1Eb) 
#' \donttest{
#' res <- trtSelThresh(fit0E, fit0Eb, fit1E, fit1Eb, 
#'                     lowRef = FALSE, toxRef = FALSE, r = 0.02, le.MCMC = 5000, plot = TRUE, 
#'                     progress.bar = "text")
#'
#' # second case: the mean risks of event in the two treatment arms are given through mcmc.lists 
#' # that correspond to their posterior distributions (see the fit man page for examples on how
#' # to generate posterior distributions manually)
#'
#' #You can summarize the results using the summary() function:
#' summary(res, method = "median")
#'
#' #You can extract the estimates and CI bounds of each indicator presented in the summary:
#' estimates(res, method = "median")
#' credibleIntervals(res)
#'
#' #Plot the decision curves (this function is time-consuming):
#' dCres <- decisionCurve(res, r = seq(0, 0.2, length.out = 6))
#'
#' #You can plot again the decision curves by applying the plot method to dCres, 
#' #this function is a lot faster than the previous one. It also has more options
#' #to customize the plot:
#' plot(dCres)
#' plot(dCres, which = 1)
#' plot(dCres, which = 2)
#' }
#' @export trtSelThresh
trtSelThresh <- function(EvtRefDist = NULL, NoEvtRefDist = NULL, EvtInnovDist = NULL, NoEvtInnovDist = NULL, mRiskRef = NULL, mRiskInnov = NULL, lowRef = TRUE, toxRef = TRUE, r = 0, le.MCMC = 1000, hessTol = 10 ^ (-6), plot = FALSE, progress.bar = NULL, seed = NULL) {
	if (!is.null(seed)) set.seed(seed)
	if (!(methods::is(EvtRefDist, "allowedFitDist"))) stop("The specified distribution for EvtRefDist is not supported or misspecified.")
	if (!(methods::is(NoEvtRefDist, "allowedFitDist"))) stop("The specified distribution for NoEvtRefDist is not supported or misspecified.")
	if (!(methods::is(EvtInnovDist, "allowedFitDist"))) stop("The specified distribution for EvtInnovDist is not supported or misspecified.")
	if (!(methods::is(NoEvtInnovDist, "allowedFitDist"))) stop("The specified distribution for NoEvtInnovDist is not supported or misspecified.")
	if (!(is.logical(lowRef))) stop("Argument lowRef must be logical.")
	if (!(is.logical(toxRef))) stop("Argument toxRef must be logical.")
	if (!(is.numeric(r))) stop("r must be numeric.")
	if (r < 0 | r >= 1) stop("r must be in [0;1[.")
	if (!(is.null(progress.bar))) {
		match.arg(progress.bar, c("text", "none"))
		if (progress.bar == "none") progress.bar <- NULL
	}
	if (sum(sapply(c(EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist), function(x) is(x, "undefined"))) != 1) stop("One distribution must be 'undefined'.")
	minMarker <- min(EvtRefDist@x, NoEvtRefDist@x, EvtInnovDist@x, NoEvtInnovDist@x)
	maxMarker <- max(EvtRefDist@x, NoEvtRefDist@x, EvtInnovDist@x, NoEvtInnovDist@x)
	pT0 <- length(c(EvtRefDist@x, NoEvtRefDist@x)) / length(c(EvtRefDist@x, NoEvtRefDist@x, EvtInnovDist@x, NoEvtInnovDist@x))
	pT1 <- length(c(EvtInnovDist@x, NoEvtInnovDist@x)) / length(c(EvtRefDist@x, NoEvtRefDist@x, EvtInnovDist@x, NoEvtInnovDist@x))
	le.MCMC <- as.integer(le.MCMC)
	if (le.MCMC <= 0) stop("Argument le.MCMC must be positive.")
	if (methods::is(EvtRefDist, "fitUserDefinedDist")) {
		res <- defineUserDist0E(EvtRefDist)
		.optimalThreshold_envEvtRefDist <- res$env
		EvtRefDist <- res$object}
	if (methods::is(NoEvtRefDist, "fitUserDefinedDist")) {
		res <- defineUserDist0Eb(NoEvtRefDist)
		.optimalThreshold_envNoEvtRefDist <- res$env
		NoEvtRefDist <- res$object}
	if (methods::is(EvtInnovDist, "fitUserDefinedDist")) {
		res <- defineUserDist1E(EvtInnovDist)
		.optimalThreshold_envEvtInnovDist <- res$env
		EvtInnovDist <- res$object}
	if (methods::is(NoEvtInnovDist, "fitUserDefinedDist")) {
		res <- defineUserDist1Eb(NoEvtInnovDist)
		.optimalThreshold_envNoEvtInnovDist <- res$env
		NoEvtInnovDist <- res$object}
	nchains <- vector(length = 1)
	if (EvtRefDist@mcmc) nchains <- length(EvtRefDist@ini)
	if (methods::is(EvtRefDist, "fitUserDefinedDist0E")) nchains <- length(EvtRefDist@mcmcList)
	if (NoEvtRefDist@mcmc) {
		if (nchains == FALSE) nchains <- length(NoEvtRefDist@ini)
		else { 
			if (nchains != length(NoEvtRefDist@ini)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (methods::is(NoEvtRefDist, "fitUserDefinedDist0Eb")) {
		if (nchains == FALSE) nchains <- length(NoEvtRefDist@mcmcList)
		else { 
			if (nchains != length(NoEvtRefDist@mcmcList)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (EvtInnovDist@mcmc) {
		if (nchains == FALSE) nchains <- length(EvtInnovDist@ini)
		else { 
			if (nchains != length(EvtInnovDist@ini)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (methods::is(EvtInnovDist, "fitUserDefinedDist1E")) {
		if (nchains == FALSE) nchains <- length(EvtInnovDist@mcmcList)
		else { 
			if (nchains != length(EvtInnovDist@mcmcList)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (NoEvtInnovDist@mcmc) {
		if (nchains == FALSE) nchains <- length(NoEvtInnovDist@ini)
		else { 
			if (nchains != length(NoEvtInnovDist@ini)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (methods::is(NoEvtInnovDist, "fitUserDefinedDist1Eb")) {
		if (nchains == FALSE) nchains <- length(NoEvtInnovDist@mcmcList)
		else { 
			if (nchains != length(NoEvtInnovDist@mcmcList)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (coda::is.mcmc.list(mRiskRef)) {
		if (nchains == FALSE) nchains <- length(mRiskRef)
		else {
			if (nchains != length(mRiskRef)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (coda::is.mcmc.list(mRiskInnov)) {
		if (nchains == FALSE) nchains <- length(mRiskInnov)
		else {
			if (nchains != length(mRiskInnov)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (nchains == FALSE) nchains <- 1
	le.MCMC <- as.integer(le.MCMC)
	do.pb <- interactive() &&  !is.null(progress.bar) && (le.MCMC*nchains) >= 100
	if (!(is.null(mRiskRef))) {
		if (!(coda::is.mcmc.list(mRiskRef))) stop("mRiskRef must be an object of class mcmc.list.")
		else {
			if (length(mRiskRef) != nchains) stop("Different number of MCMC chains were supplied for mRiskRef.")
			else {
				if (any(sapply(mRiskRef, function(i) coda::nvar(i)) != 1)) stop("Multiple variables supplied in a chain of 'mRiskRef' argument.")
				else { 
					mcmc.list.r0 <- mRiskRef
					mcmcChainr0 <- lapply(mRiskRef, function(i) as.vector(i))
				}
			}
		}
	}
	else {
		mcmcChainr0 <- replicate(nchains, rbeta(le.MCMC, 0.5 + EvtRefDist@n, 0.5 + NoEvtRefDist@n), simplify = FALSE)
	}
	if (!(is.null(mRiskInnov))) {
		if (!(coda::is.mcmc.list(mRiskInnov))) stop("mRiskInnov must be an object of class mcmc.list.")
		else {
			if (length(mRiskInnov) != nchains) stop("Different number of MCMC chains were supplied for mRiskInnov.")
			else {
				if (any(sapply(mRiskInnov, function(i) coda::nvar(i)) != 1)) stop("Multiple variables supplied in a chain of 'mRiskInnov' argument.")
				else {
					mcmc.list.r1 <- mRiskInnov
					mcmcChainr1 <- lapply(mRiskInnov, function(i) as.vector(i))
				}
			}
		}
	}
	else {
		mcmcChainr1 <- replicate(nchains, stats::rbeta(le.MCMC, 0.5 + EvtInnovDist@n, 0.5 + NoEvtInnovDist@n), simplify = FALSE)
	}
	if (class(EvtRefDist) == "undefined") {
		if (do.pb) cat("\n-----------------------------\nSampling in NoEvtRefDist\n-----------------------------\n")
		if (NoEvtRefDist@mcmc | methods::is(NoEvtRefDist, "fitUserDefinedDist0Eb")) mcmcChainNoEvtRef <- samplePosteriorDist(NoEvtRefDist, le.MCMC, do.pb, seed)
		else mcmcChainNoEvtRef <- samplePosteriorDist(NoEvtRefDist, le.MCMC, nchains)
		if (do.pb) cat("\n-----------------------------\nSampling in EvtInnovDist\n-----------------------------\n")
		if (EvtInnovDist@mcmc | methods::is(EvtInnovDist, "fitUserDefinedDist1E")) mcmcChainEvtInnov <- samplePosteriorDist(EvtInnovDist, le.MCMC, do.pb, seed)
		else mcmcChainEvtInnov <- samplePosteriorDist(EvtInnovDist, le.MCMC, nchains)
		if (do.pb) cat("\n-----------------------------\nSampling in NoEvtInnovDist\n-----------------------------\n")
		if (NoEvtInnovDist@mcmc | methods::is(NoEvtInnovDist, "fitUserDefinedDist1Eb")) mcmcChainNoEvtInnov <- samplePosteriorDist(NoEvtInnovDist, le.MCMC, do.pb, seed)
		else mcmcChainNoEvtInnov <- samplePosteriorDist(NoEvtInnovDist, le.MCMC, nchains)
		mcmcChainEvtRef <- sapply(1:nchains, function(i) mapply(compoundEvtRefDist, mcmcChainNoEvtRef[[i]], mcmcChainEvtInnov[[i]], mcmcChainNoEvtInnov[[i]], mcmcChainr0[[i]], mcmcChainr1[[i]], SIMPLIFY = FALSE), simplify = FALSE)
	}
	else {
		if (class(NoEvtRefDist) == "undefined") {
			if (do.pb) cat("\n-----------------------------\nSampling in EvtRefDist\n-----------------------------\n")
			if (EvtRefDist@mcmc | methods::is(EvtRefDist, "fitUserDefinedDist0E")) mcmcChainEvtRef <- samplePosteriorDist(EvtRefDist, le.MCMC, do.pb, seed)
			else mcmcChainEvtRef <- samplePosteriorDist(EvtRefDist, le.MCMC, nchains)
			if (do.pb) cat("\n-----------------------------\nSampling in EvtInnovDist\n-----------------------------\n")
			if (EvtInnovDist@mcmc | methods::is(EvtInnovDist, "fitUserDefinedDist1E")) mcmcChainEvtInnov <- samplePosteriorDist(EvtInnovDist, le.MCMC, do.pb, seed)
			else mcmcChainEvtInnov <- samplePosteriorDist(EvtInnovDist, le.MCMC, nchains)
			if (do.pb) cat("\n-----------------------------\nSampling in NoEvtInnovDist\n-----------------------------\n")
			if (NoEvtInnovDist@mcmc | methods::is(NoEvtInnovDist, "fitUserDefinedDist1Eb")) mcmcChainNoEvtInnov <- samplePosteriorDist(NoEvtInnovDist, le.MCMC, do.pb, seed)
			else mcmcChainNoEvtInnov <- samplePosteriorDist(NoEvtInnovDist, le.MCMC, nchains)
			mcmcChainNoEvtRef <- sapply(1:nchains, function(i) mapply(compoundNoEvtRefDist, mcmcChainEvtRef[[i]], mcmcChainEvtInnov[[i]], mcmcChainNoEvtInnov[[i]], mcmcChainr0[[i]], mcmcChainr1[[i]], SIMPLIFY = FALSE), simplify = FALSE)
		}
		else {
			if (class(EvtInnovDist) == "undefined") {
				if (do.pb) cat("\n-----------------------------\nSampling in EvtRefDist\n-----------------------------\n")
				if (EvtRefDist@mcmc | methods::is(EvtRefDist, "fitUserDefinedDist0E")) mcmcChainEvtRef <- samplePosteriorDist(EvtRefDist, le.MCMC, do.pb, seed)
				else mcmcChainEvtRef <- samplePosteriorDist(EvtRefDist, le.MCMC, nchains)
				if (do.pb) cat("\n-----------------------------\nSampling in NoEvtRefDist\n-----------------------------\n")
				if (NoEvtRefDist@mcmc | methods::is(NoEvtRefDist, "fitUserDefinedDist0Eb")) mcmcChainNoEvtRef <- samplePosteriorDist(NoEvtRefDist, le.MCMC, do.pb, seed)
				else mcmcChainNoEvtRef <- samplePosteriorDist(NoEvtRefDist, le.MCMC, nchains)
				if (do.pb) cat("\n-----------------------------\nSampling in NoEvtInnovDist\n-----------------------------\n")
				if (NoEvtInnovDist@mcmc | methods::is(NoEvtInnovDist, "fitUserDefinedDist1Eb")) mcmcChainNoEvtInnov <- samplePosteriorDist(NoEvtInnovDist, le.MCMC, do.pb, seed)
				else mcmcChainNoEvtInnov <- samplePosteriorDist(NoEvtInnovDist, le.MCMC, nchains)
				mcmcChainEvtInnov <- sapply(1:nchains, function(i) mapply(compoundEvtInnovDist, mcmcChainEvtRef[[i]], mcmcChainNoEvtRef[[i]], mcmcChainNoEvtInnov[[i]], mcmcChainr0[[i]], mcmcChainr1[[i]], SIMPLIFY = FALSE), simplify = FALSE)
			}
			else {
				if (class(NoEvtInnovDist) == "undefined") {
					if (do.pb) cat("\n-----------------------------\nSampling in EvtRefDist\n-----------------------------\n")
					if (EvtRefDist@mcmc | methods::is(EvtRefDist, "fitUserDefinedDist0E")) mcmcChainEvtRef <- samplePosteriorDist(EvtRefDist, le.MCMC, do.pb, seed)
					else mcmcChainEvtRef <- samplePosteriorDist(EvtRefDist, le.MCMC, nchains)
					if (do.pb) cat("\n-----------------------------\nSampling in NoEvtRefDist\n-----------------------------\n")
					if (NoEvtRefDist@mcmc | methods::is(NoEvtRefDist, "fitUserDefinedDist0Eb")) mcmcChainNoEvtRef <- samplePosteriorDist(NoEvtRefDist, le.MCMC, do.pb, seed)
					else mcmcChainNoEvtRef <- samplePosteriorDist(NoEvtRefDist, le.MCMC, nchains)
					if (do.pb) cat("\n-----------------------------\nSampling in EvtInnovDist\n-----------------------------\n")
					if (EvtInnovDist@mcmc | methods::is(EvtInnovDist, "fitUserDefinedDist1E")) mcmcChainEvtInnov <- samplePosteriorDist(EvtInnovDist, le.MCMC, do.pb, seed)
					else mcmcChainEvtInnov <- samplePosteriorDist(EvtInnovDist, le.MCMC, nchains)
					mcmcChainNoEvtInnov <- sapply(1:nchains, function(i) mapply(compoundNoEvtInnovDist, mcmcChainEvtRef[[i]], mcmcChainNoEvtRef[[i]], mcmcChainEvtInnov[[i]], mcmcChainr0[[i]], mcmcChainr1[[i]], SIMPLIFY = FALSE), simplify = FALSE)
				}
				else {
					stop("One distribution must be of class undefined to perform the analysis")
				}
			}
		}
	}
	if (lowRef == TRUE & toxRef == TRUE) {
		r <- (-r)
		coefU <- 1
	}
	else {
		if (lowRef == TRUE & toxRef == FALSE) {
			coefU <-1
		}
		else {
			if (lowRef == FALSE & toxRef == TRUE) {
				r <- (-r)
				coefU <- (-1)
			}
			else {
				coefU <- (-1)
			}
		}
	}
	if (do.pb) {
		cat("\nEstimating optimal threshold\n\n")
		pb <- utils::txtProgressBar(1, le.MCMC * nchains, initial = 1, style = 3, width = 50, char = "*")
	}
	else pb <- NULL
	mcmcChainThreshold <- sapply(1:nchains, function(i) mapply(maxUtility, ((i - 1) * le.MCMC + 1):(i * le.MCMC), mcmcChainEvtRef[[i]], mcmcChainNoEvtRef[[i]], mcmcChainEvtInnov[[i]], mcmcChainNoEvtInnov[[i]], mcmcChainr0[[i]], mcmcChainr1[[i]], MoreArgs = list(pT0 = pT0, pT1 = pT1, r = r, coefU = coefU, minMarker = minMarker, maxMarker = maxMarker, utility = utility, hessTol = hessTol, pb = pb), SIMPLIFY = TRUE), simplify = FALSE)
	if (do.pb) close(pb)
	listTabRes <- NULL
	if (coda::is.mcmc.list(mRiskRef)) listTabRes <- mRiskRef
	if (coda::is.mcmc.list(mRiskInnov)) listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), mRiskRef, mRiskInnov, SIMPLIFY = FALSE))
	if (class(EvtRefDist) == "undefined") {
		if (NoEvtRefDist@mcmc) {
			mcmc.list.DistNoEvtRef <- lapply(mcmcChainNoEvtRef, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) {nameS <- methods::slotNames(i); sapply(nameS,function(n) getElement(i, n))}))); coda::varnames(M) <- paste0(coda::varnames(M), "_NoEvtRef"); return(M)})
			if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistNoEvtRef)
			else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistNoEvtRef, SIMPLIFY = FALSE))
		}
		if (EvtInnovDist@mcmc) {
			mcmc.list.DistEvtInnov <- lapply(mcmcChainEvtInnov, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) sapply(methods::slotNames(i), function(n) getElement(i, n))))); coda::varnames(M) <- paste0(coda::varnames(M), "_EvtInnov"); return(M)})
			if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistEvtInnov)
			else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistEvtInnov, SIMPLIFY = FALSE))
		}
		if (NoEvtInnovDist@mcmc) {
			mcmc.list.DistNoEvtInnov <- lapply(mcmcChainNoEvtInnov, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) sapply(methods::slotNames(i), function(n) getElement(i, n))))); coda::varnames(M) <- paste0(coda::varnames(M), "_NoEvtInnov"); return(M)})
			if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistNoEvtInnov)
			else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistNoEvtInnov, SIMPLIFY = FALSE))
		}
	}
	else {
		if (class(NoEvtRefDist) == "undefined") {
			if (EvtRefDist@mcmc) {
				mcmc.list.DistEvtRef <- lapply(mcmcChainEvtRef, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) {nameS <- methods::slotNames(i); sapply(nameS, function(n) getElement(i, n))}))); coda::varnames(M) <- paste0(coda::varnames(M), "_EvtRef"); return(M)})
				if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistEvtRef)
				else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistEvtRef, SIMPLIFY = FALSE))
			}
			if (EvtInnovDist@mcmc) {
				mcmc.list.DistEvtInnov <- lapply(mcmcChainEvtInnov, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) sapply(methods::slotNames(i), function(n) getElement(i, n))))); coda::varnames(M) <- paste0(coda::varnames(M), "_EvtInnov"); return(M)})
				if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistEvtInnov)
				else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistEvtInnov, SIMPLIFY = FALSE))
			}
			if (NoEvtInnovDist@mcmc) {
				mcmc.list.DistNoEvtInnov <- lapply(mcmcChainNoEvtInnov, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) sapply(methods::slotNames(i), function(n) getElement(i, n))))); coda::varnames(M) <- paste0(coda::varnames(M), "_NoEvtInnov"); return(M)})
				if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistNoEvtInnov)
				else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistNoEvtInnov, SIMPLIFY = FALSE))
			}			
		}
		else {
			if(class(EvtInnovDist) == "undefined") {
				if (EvtRefDist@mcmc) {
					mcmc.list.DistEvtRef <- lapply(mcmcChainEvtRef, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) {nameS <- methods::slotNames(i); sapply(nameS, function(n) getElement(i, n))}))); coda::varnames(M) <- paste0(coda::varnames(M), "_EvtRef"); return(M)})
					if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistEvtRef)
					else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistEvtRef, SIMPLIFY = FALSE))
				}
				if (NoEvtRefDist@mcmc) {
					mcmc.list.DistNoEvtRef <- lapply(mcmcChainNoEvtRef, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) {nameS <- methods::slotNames(i); sapply(nameS, function(n) getElement(i, n))}))); coda::varnames(M) <- paste0(coda::varnames(M), "_NoEvtRef"); return(M)})
					if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistNoEvtRef)
					else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistNoEvtRef, SIMPLIFY = FALSE))
				}
				if (NoEvtInnovDist@mcmc) {
					mcmc.list.DistNoEvtInnov <- lapply(mcmcChainNoEvtInnov, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) sapply(methods::slotNames(i), function(n) getElement(i, n))))); coda::varnames(M) <- paste0(coda::varnames(M), "_NoEvtInnov"); return(M)})
					if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistNoEvtInnov)
					else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistNoEvtInnov, SIMPLIFY = FALSE))
				}	
			}
			else {
				if (EvtRefDist@mcmc) {
					mcmc.list.DistEvtRef <- lapply(mcmcChainEvtRef, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) {nameS <- methods::slotNames(i); sapply(nameS, function(n) getElement(i, n))}))); coda::varnames(M) <- paste0(coda::varnames(M), "_EvtRef"); return(M)})
					if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistEvtRef)
					else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistEvtRef, SIMPLIFY = FALSE))
				}
				if (NoEvtRefDist@mcmc) {
					mcmc.list.DistNoEvtRef <- lapply(mcmcChainNoEvtRef, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) {nameS <- methods::slotNames(i); sapply(nameS, function(n) getElement(i, n))}))); coda::varnames(M) <- paste0(coda::varnames(M), "_NoEvtRef"); return(M)})
					if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistNoEvtRef)
					else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistNoEvtRef, SIMPLIFY = FALSE))
				}
				if (EvtInnovDist@mcmc) {
					mcmc.list.DistEvtInnov <- lapply(mcmcChainEvtInnov, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) sapply(methods::slotNames(i), function(n) getElement(i, n))))); coda::varnames(M) <- paste0(coda::varnames(M), "_EvtInnov"); return(M)})
					if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistEvtInnov)
					else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistEvtInnov, SIMPLIFY = FALSE))
				}
			}
		}
	}
	mcmcMarkerBasedRisk <- sapply(1:nchains, function(i) unname(mapply(markerBasedRisk, mcmcChainEvtRef[[i]], mcmcChainNoEvtRef[[i]], mcmcChainEvtInnov[[i]], mcmcChainNoEvtInnov[[i]], mcmcChainr0[[i]], mcmcChainr1[[i]], mcmcChainThreshold[[i]], MoreArgs = list(lowRef = lowRef))), simplify = FALSE)
	mcmcChainThreshold <- unname(unlist(mcmcChainThreshold))
	mcmcMarkerBasedRiskRef <- unlist(lapply(1:nchains, function(i) mcmcMarkerBasedRisk[[i]][1, ]))
	mcmcMarkerBasedRiskInnov <- unlist(lapply(1:nchains, function(i) mcmcMarkerBasedRisk[[i]][2, ]))
	percentNA <- length(mcmcChainThreshold[is.na(mcmcChainThreshold)]) / (nchains * le.MCMC)
	if (plot == TRUE) {
		opar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(opar))
		if (any(sapply(c(EvtRefDist, NoEvtRefDist, EvtInnovDist, NoEvtInnovDist), function(x) x@mcmc))) {
			devAskNewPage(ask = TRUE)
			plot(listTabRes, ask = TRUE)
			coda::autocorr.plot(listTabRes, ask = TRUE)
			graphics::par(mfrow = c(1, 1), lwd = 2)
			graphics::hist(mcmcChainThreshold, yaxs = "i", main = "MCMC sample distribution of optimal threshold", col = "gray85", border = "darkgrey", xlab = "Optimal threshold estimate", freq = FALSE, breaks = seq(min(mcmcChainThreshold, na.rm = TRUE), max(mcmcChainThreshold, na.rm = TRUE), length.out = 20))
			colAlpha <- t(grDevices::col2rgb("gray85", alpha = TRUE))
			colAlpha[4] <- 160
			graphics::polygon(c(stats::density(mcmcChainThreshold, na.rm = TRUE)$x, rev(stats::density(mcmcChainThreshold, na.rm = TRUE)$x)), c(stats::density(mcmcChainThreshold, na.rm = TRUE)$y, rep(0, length(stats::density(mcmcChainThreshold, na.rm = TRUE)$x))), col = grDevices::rgb(colAlpha[1], colAlpha[2], colAlpha[3], colAlpha[4], maxColorValue = 255), border = NA)
			graphics::lines(stats::density(mcmcChainThreshold, na.rm = TRUE), col = "darkgrey", lwd = 2)
			graphics::box(bty = "L", lwd = 1)
			devAskNewPage(ask = FALSE)
		}
		else {
			graphics::par(mfrow = c(1, 1), lwd = 2)
			graphics::hist(mcmcChainThreshold, yaxs = "i", main = "MCMC sample distribution of optimal threshold", col = "gray85", border = "darkgrey", xlab = "Optimal threshold estimate", freq = FALSE, breaks = seq(min(mcmcChainThreshold, na.rm = TRUE), max(mcmcChainThreshold, na.rm = TRUE), length.out = 20))
			colAlpha <- t(grDevices::col2rgb("gray85", alpha = TRUE))
			colAlpha[4] <- 160
			graphics::polygon(c(stats::density(mcmcChainThreshold, na.rm = TRUE)$x, rev(stats::density(mcmcChainThreshold, na.rm = TRUE)$x)), c(stats::density(mcmcChainThreshold, na.rm = TRUE)$y, rep(0, length(stats::density(mcmcChainThreshold, na.rm = TRUE)$x))), col = grDevices::rgb(colAlpha[1], colAlpha[2], colAlpha[3], colAlpha[4], maxColorValue = 255), border = NA)
			graphics::lines(stats::density(mcmcChainThreshold, na.rm = TRUE), col = "darkgrey", lwd = 2)
			graphics::box(bty = "L", lwd = 1)
		}
	}
	paraNamesUserDefined <- vector("list", 4)
	names(paraNamesUserDefined) <- c("EvtRef", "NoEvtRef", "EvtInnov", "NoEvtInnov")
	cdfUserDefined <- vector("list", 4)
	names(cdfUserDefined) <- c("EvtRef", "NoEvtRef", "EvtInnov", "NoEvtInnov")
	gradientUserDefined <- vector("list", 4)
	names(gradientUserDefined) <- c("EvtRef", "NoEvtRef", "EvtInnov", "NoEvtInnov")
	hessianUserDefined <- vector("list", 4)
	names(hessianUserDefined) <- c("EvtRef", "NoEvtRef", "EvtInnov", "NoEvtInnov")
	if (any(methods::is(EvtRefDist, "fitUserDefinedDist0E"), methods::is(NoEvtRefDist, "fitUserDefinedDist0Eb"), methods::is(EvtInnovDist, "fitUserDefinedDist1E"), methods::is(NoEvtInnovDist, "fitUserDefinedDist1Eb"))) {
		if (methods::is(EvtRefDist, "fitUserDefinedDist0E")) {
			paraNamesUserDefined$EvtRef <- EvtRefDist@paraNames
			cdfUserDefined$EvtRef <- methods::getMethod(cdf, "userDefinedDist0E")
			gradientUserDefined$EvtRef <- methods::getMethod(gradient, "userDefinedDist0E")
			hessianUserDefined$EvtRef <- methods::getMethod(hessian, "userDefinedDist0E")
		}
		if (methods::is(NoEvtRefDist, "fitUserDefinedDist0Eb")) {
			paraNamesUserDefined$NoEvtRef <- NoEvtRefDist@paraNames
			cdfUserDefined$NoEvtRef <- methods::getMethod(cdf, "userDefinedDist0Eb")
			gradientUserDefined$NoEvtRef <- methods::getMethod(gradient, "userDefinedDist0Eb")
			hessianUserDefined$NoEvtRef <- methods::getMethod(hessian, "userDefinedDist0Eb")
		}
		if (methods::is(EvtInnovDist, "fitUserDefinedDist1E")) {
			paraNamesUserDefined$EvtInnov <- EvtInnovDist@paraNames
			cdfUserDefined$EvtInnov <- methods::getMethod(cdf, "userDefinedDist1E")
			gradientUserDefined$EvtInnov <- methods::getMethod(gradient, "userDefinedDist1E")
			hessianUserDefined$EvtInnov <- methods::getMethod(hessian, "userDefinedDist1E")
		}
		if (methods::is(NoEvtInnovDist, "fitUserDefinedDist1Eb")) {
			paraNamesUserDefined$NoEvtInnov <- NoEvtInnovDist@paraNames
			cdfUserDefined$NoEvtInnov <- methods::getMethod(cdf, "userDefinedDist1Eb")
			gradientUserDefined$NoEvtInnov <- methods::getMethod(gradient, "userDefinedDist1Eb")
			hessianUserDefined$NoEvtInnov <- methods::getMethod(hessian, "userDefinedDist1Eb")
		}
	}
	res <- methods::new("trtSelOptThresh", optThresh = mcmcChainThreshold, r0 = unlist(mcmcChainr0), r1 = unlist(mcmcChainr1), xEvtRef = EvtRefDist@x, xNoEvtRef = NoEvtRefDist@x, xEvtInnov = EvtInnovDist@x, xNoEvtInnov = NoEvtInnovDist@x, lowRef = lowRef, toxRef = toxRef, markerBasedRiskRef = mcmcMarkerBasedRiskRef, markerBasedRiskInnov = mcmcMarkerBasedRiskInnov, mcmcChainEvtRef = mcmcChainEvtRef, mcmcChainNoEvtRef = mcmcChainNoEvtRef, mcmcChainEvtInnov = mcmcChainEvtInnov, mcmcChainNoEvtInnov = mcmcChainNoEvtInnov, tabMCMCChain = listTabRes, paraNamesUserDefined = paraNamesUserDefined, cdfUserDefined = cdfUserDefined, gradientUserDefined = gradientUserDefined, hessianUserDefined = hessianUserDefined, percentNA = percentNA)
	if (any(methods::is(EvtRefDist, "fitUserDefinedDist0E"), methods::is(NoEvtRefDist, "fitUserDefinedDist0Eb"), methods::is(EvtInnovDist, "fitUserDefinedDist1E"), methods::is(NoEvtInnovDist, "fitUserDefinedDist1Eb"))) {	
		if (methods::is(EvtRefDist, "fitUserDefinedDist0E")) {
			methods::removeMethod(cdf, "userDefinedDist0E", where = .optimalThreshold_envEvtRefDist)
			methods::removeMethod(gradient, "userDefinedDist0E", where = .optimalThreshold_envEvtRefDist)
			methods::removeMethod(hessian, "userDefinedDist0E", where = .optimalThreshold_envEvtRefDist)
			methods::removeMethod(samplePosteriorDist, "fitUserDefinedDist0E", where = .optimalThreshold_envEvtRefDist)
			methods::removeClass("fitUserDefinedDist0E", where = .optimalThreshold_envEvtRefDist)
			methods::removeClass("userDefinedDist0E", where = .optimalThreshold_envEvtRefDist)			
		}
		if (methods::is(NoEvtRefDist, "fitUserDefinedDist0Eb")) {
			methods::removeMethod(cdf, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
			methods::removeMethod(gradient, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
			methods::removeMethod(hessian, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
			methods::removeMethod(samplePosteriorDist, "fitUserDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
			methods::removeClass("fitUserDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
			methods::removeClass("userDefinedDist0Eb", where = .optimalThreshold_envNoEvtRefDist)
		}
		if (methods::is(EvtInnovDist, "fitUserDefinedDist1E")) {
			methods::removeMethod(cdf, "userDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
			methods::removeMethod(gradient, "userDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
			methods::removeMethod(hessian, "userDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
			methods::removeMethod(samplePosteriorDist, "fitUserDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
			methods::removeClass("fitUserDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
			methods::removeClass("userDefinedDist1E", where = .optimalThreshold_envEvtInnovDist)
		}
		if (methods::is(NoEvtInnovDist, "fitUserDefinedDist1Eb")) {
			methods::removeMethod(cdf, "userDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
			methods::removeMethod(gradient, "userDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
			methods::removeMethod(hessian, "userDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
			methods::removeMethod(samplePosteriorDist, "fitUserDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
			methods::removeClass("fitUserDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
			methods::removeClass("userDefinedDist1Eb", where = .optimalThreshold_envNoEvtInnovDist)
		}
	}
	return(res)
}

explicitDiagThreshold <- function(EvtDist, NoEvtDist, p, r, coefU) {
	if (methods::is(EvtDist, "normalDist") & methods::is(NoEvtDist, "normalDist")) {
		a <- EvtDist@mu - NoEvtDist@mu
		W <- (r / (1 - r)) * ((1 - p) / p)
		if (EvtDist@sd != NoEvtDist@sd) {
			b <- EvtDist@sd / NoEvtDist@sd
			thresholdEst <- (NoEvtDist@mu*(b ^ 2 - 1) - a + b * sqrt(a ^ 2 + (b ^ 2 - 1) * NoEvtDist@sd ^ 2 * log(b ^ 2 * W ^ 2))) / (b ^ 2 - 1)
			return(thresholdEst)
		}
		else {
			thresholdEst <- (NoEvtDist@sd ^ 2 * log(W ^ 2) + EvtDist@mu ^ 2 - NoEvtDist@mu ^ 2) / (2 * a)
			return(thresholdEst)
		}	
	}
	return(NA)
}
maxDiagUtility <- function(i, EvtDist, NoEvtDist, prev, r, minMarker, maxMarker, utility, hessTol, coefU, pb) {
	if (!is.null(pb)) setTxtProgressBar(pb, i)
	thresholdEst <- explicitDiagThreshold(EvtDist, NoEvtDist, prev, r)
	if (!(is.na(thresholdEst))) return(thresholdEst)
	gridStartVal <- seq(minMarker, maxMarker, length.out = 1000)
	startVal <- defineStartVal(gridStartVal, utility, hessTol, r, prev, EvtDist, NoEvtDist, coefU)
	if (startVal[[2]] %in% c("Monotonic function", "No minimum found")) {
		warning("Utility function is not convex.")
		return(NA)
	}
	thresholdEst <- stats::nlm(utility, startVal[[1]], r = r, prev = prev, EvtDist = EvtDist, NoEvtDist = NoEvtDist, coefU = coefU, hessian = FALSE, check.analyticals = TRUE)
	Uopt <- utility(thresholdEst$estimate, r, prev, EvtDist, NoEvtDist, coefU)
	if (is.nan(attr(Uopt, "hessian"))) {
		warning("Hessian of the utility function is not defined at threshold estimate.")
		return(NA)
	}
	if (abs(attr(Uopt, "hessian")) < hessTol) {
		warning("Hessian of the utility function at threshold estimate is too close to zero to correspond to a minimum estimate.")
		return(NA)
	}
	return(thresholdEst$estimate)
}

#' Estimation of the optimal threshold of a diagnostic marker
#'
#' This function produces a sample of the posterior distribution of the optimal threshold of a diagnostic marker. The optimal threshold is defined as the marker value that maximized the utility of using the marker to make the diagnosis and treat the patient (treat or not the patient). The utility function takes into account the proportions of patients well classified and miss-classified (through the sensitivity and specificity), the prevalence of the disease in the target population, and the cost and benefits of treating wrongly or rightly the subject.
#' To calculate the utility function, the user needs to specify:
#' \itemize{
#' \item the distribution of the marker in the subject with and without the disease (see the \code{\link[optimalThreshold]{fit}} function)
#' \item the prevalence of the disease in the target population
#' \item the cost of treating subject without the disease and the benefit of treating a patient with the disease (see Details). 
#' }
#' The optimal threshold and its credible interval are calculated using a Monte Carlo approach. 
#'
#' @param EvtDist an object of class allowedFitDist that summarizes the distribution fitted to the marker values of patients with the disease of interest. This class of objects is obtained using the \code{\link[optimalThreshold]{fit}} function.
#' @param NoEvtDist an object of class allowedFitDist that summarizes the distribution fitted to the marker values of patients without the disease of interest. This class of objects is obtained using the \code{\link[optimalThreshold]{fit}} function.
#' @param p the prevalence of the disease in the target population.
#' @param r the risk threshold preference (see Details).
#' @param lowEvt logical argument that specifies whether low values of the marker are associated with the presence of the disease or not.
#' @param le.MCMC length of the desired MCMC chain.
#' @param hessTol a numeric value used in the optimization algorithm to control the tolerance threshold of the hessian value at the optimal threshold estimate.
#' @param plot a logical value indicating whether routine graphics must be produced.
#' @param progress.bar a character string indicating whether the user wishes to print a progress bar during the function process.
#' @param seed a numerical value used to fix the random seed.
#' @return Returns an object of class \code{diagOptThresh}.
#' @details The r value can be defined as the probability of disease above which a patient or a physician would accept the treatment. The value (1-r)/r can be interpreted as the NB/NC ratio, i.e. the number of subjects without the disease a physician would accept to treat wrongly to be able to detect and treat one diseased patient.
#' @section References :
#' Subtil, F, and Rabilloud. A Bayesian method to estimate the optimal threshold of a longitudinal marker. \emph{Biometrical Journal}. 2010; 52(3): 333-347.
#' @examples
#' #Simulating data from two gaussian distributions:
#' xE <- rnorm(100) # distribution of the marker in diseased patients
#' xEb <- rnorm(400, 2) # distribution of the marker in the subjects without the disease
#'
#' #When working with real data. You can check the decision rule (whether low or high 
#' #value of the marker are associated with the disease) using the densCurves function:
#' densCurves(x0 = xEb, x1 = xE, type = "diagnostic")
#'
#' #Fit normal distributions on the two groups:
#' fitE <- fit(xE, "norm")
#' fitEb <- fit(xEb, "norm")
#'
#' #Apply the main function to estimate the optimal threshold:
#' \donttest{
#' res <- diagThresh(fitE, fitEb, p = 0.2, r = 0.3, lowEvt = TRUE, le.MCMC = 5000, 
#'                   plot = TRUE, progress.bar = "text")
#'
#' #You can summarize the results using the summary() function:
#' summary(res,method = "median")
#'
#' #You can extract the estimates and CI bounds of each indicator presented in the summary:
#' estimates(res, method = "median")
#' credibleIntervals(res)
#'
#' #Plot the decision curves (this function is time-consuming):
#' dCres <- decisionCurve(res, r = seq(0, 0.5, length.out = 10))
#'
#' #You can plot again the decision curves by applying the plot method to dCres, 
#' #this function is a lot faster than the previous one. It also has more options
#' #to customize the plot:
#' plot(dCres)
#' }
#' @export diagThresh
diagThresh <- function(EvtDist = NULL, NoEvtDist = NULL, p, r, lowEvt = FALSE, le.MCMC = 1000, hessTol = 10 ^ (-6), plot = FALSE, progress.bar = NULL, seed = NULL) {
	if (!is.null(seed)) set.seed(seed)
	if (!(methods::is(EvtDist, "allowedFitDist"))) stop("The specified distribution for EvtDist is not supported or misspecified.")
	if (!(methods::is(NoEvtDist, "allowedFitDist"))) stop("The specified distribution for NoEvtDist is not supported or misspecified.")
	if (methods::is(EvtDist, "undefined")) stop("The specified distribution for EvtDist cannot be undefined.")
	if (methods::is(NoEvtDist, "undefined")) stop("The specified distribution for NoEvtDist cannot be undefined.")
	if (missing(p)) stop("Prevalence 'p' is missing.")
	if (!(is.numeric(p))) stop("p must be numeric.")
	if (p < 0 | p > 1) stop("p must be in [0;1].")
	if (missing(r)) stop("Risk threshold 'r' is missing.")
	if (!(is.numeric(r))) stop("r must be numeric.")
	if (r < 0 | r > 1) stop("r must be in [0;1].")
	if (!(is.null(progress.bar))) {
		match.arg(progress.bar, c("text", "none"))
		if (progress.bar == "none") progress.bar <- NULL
	}
	if (!(is.logical(lowEvt))) stop("lowEvt must be logical.")
	minMarker <- min(EvtDist@x, NoEvtDist@x)
	maxMarker <- max(EvtDist@x, NoEvtDist@x)
	if (lowEvt) coefU <- (-1)
	else coefU <- 1
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
	le.MCMC <- as.integer(le.MCMC)
	if (le.MCMC <= 0) stop("Argument le.MCMC must be positive.")
	if (methods::is(EvtDist, "fitUserDefinedDist")) {
		res <- defineUserDist0E(EvtDist)
		.optimalThreshold_envEvtDist <- res$env
		EvtDist <- res$object}
	if (methods::is(NoEvtDist, "fitUserDefinedDist")) {
		res <- defineUserDist0Eb(NoEvtDist)
		.optimalThreshold_envNoEvtDist <- res$env
		NoEvtDist <- res$object}
	nchains <- vector(length = 1)
	if (EvtDist@mcmc) nchains <- length(EvtDist@ini)
	if (methods::is(EvtDist, "fitUserDefinedDist0E")) nchains <- length(EvtDist@mcmcList)
	if (NoEvtDist@mcmc) {
		if (nchains == FALSE) nchains <- length(NoEvtDist@ini)
		else { 
			if (nchains != length(NoEvtDist@ini)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (methods::is(NoEvtDist, "fitUserDefinedDist0Eb")) {
		if (nchains == FALSE) nchains <- length(NoEvtDist@mcmcList)
		else { 
			if (nchains != length(NoEvtDist@mcmcList)) stop("Different number of MCMC chains were supplied.")
		}
	}
	if (nchains == FALSE) nchains <- 1
	do.pb <- interactive() &&  !is.null(progress.bar) && (le.MCMC*nchains) >= 100
	if (do.pb) cat("\n-----------------------------\nSampling in EvtDist\n-----------------------------\n")
	if (EvtDist@mcmc | methods::is(EvtDist, "fitUserDefinedDist0E")) mcmcChainEvt <- samplePosteriorDist(EvtDist, le.MCMC, do.pb, seed)
	else mcmcChainEvt <- samplePosteriorDist(EvtDist, le.MCMC, nchains)
	if (do.pb) cat("\n-----------------------------\nSampling in NoEvtDist\n-----------------------------\n")
	if (NoEvtDist@mcmc | methods::is(NoEvtDist, "fitUserDefinedDist0Eb")) mcmcChainNoEvt <- samplePosteriorDist(NoEvtDist, le.MCMC, do.pb, seed)
	else mcmcChainNoEvt <- samplePosteriorDist(NoEvtDist, le.MCMC, nchains)
	if (do.pb) {
		cat("\nEstimating optimal threshold\n\n")
		pb <- utils::txtProgressBar(1, le.MCMC * nchains, initial = 1, style = 3, width = 50, char = "*")
	}
	else pb <- NULL
	mcmcChainThreshold <- sapply(1:nchains, function(i) mapply(maxDiagUtility, ((i - 1) * le.MCMC + 1):(i * le.MCMC), mcmcChainEvt[[i]], mcmcChainNoEvt[[i]], MoreArgs = list(prev = p, r = r, minMarker = minMarker, maxMarker = maxMarker, utility = diagEB, hessTol = hessTol, coefU = coefU, pb = pb), SIMPLIFY = TRUE), simplify = FALSE)
	mcmcChainThreshold <- unname(unlist(mcmcChainThreshold))
	if (do.pb) close(pb)
	listTabRes <- NULL
	if (EvtDist@mcmc) {
		mcmc.list.DistEvt <- lapply(mcmcChainEvt, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) {nameS <- methods::slotNames(i); sapply(nameS, function(n) getElement(i, n))}))); coda::varnames(M) <- paste0(coda::varnames(M), "_Evt"); return(M)})
		listTabRes <- coda::mcmc.list(mcmc.list.DistEvt)
	}
	if (NoEvtDist@mcmc) {
		mcmc.list.DistNoEvt <- lapply(mcmcChainNoEvt, function(chain) {M <- coda::mcmc(t(sapply(chain, function(i) {nameS <- methods::slotNames(i); sapply(nameS, function(n) getElement(i, n))}))); coda::varnames(M) <- paste0(coda::varnames(M), "_NoEvt"); return(M)})
		if (is.null(listTabRes)) listTabRes <- coda::mcmc.list(mcmc.list.DistNoEvt)
		else listTabRes <- coda::mcmc.list(mapply(function(...) coda::mcmc(cbind(...)), listTabRes, mcmc.list.DistNoEvt, SIMPLIFY = FALSE))
	}
	percentNA <- length(mcmcChainThreshold[is.na(mcmcChainThreshold)]) / (nchains * le.MCMC)
	if (plot == TRUE) {
		opar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(opar))
		if (any(sapply(c(EvtDist, NoEvtDist), function(x) x@mcmc))) {
			devAskNewPage(ask = TRUE)
			plot(listTabRes, ask = TRUE)
			coda::autocorr.plot(listTabRes, ask = TRUE)
			graphics::par(mfrow = c(1, 1), lwd = 2)
			graphics::hist(mcmcChainThreshold, yaxs = "i", main = "MCMC sample distribution of optimal threshold", col = "gray85", border = "darkgrey", xlab = "Optimal threshold estimate", freq = FALSE, breaks = seq(min(mcmcChainThreshold, na.rm = TRUE), max(mcmcChainThreshold, na.rm = TRUE), length.out = 20))
			colAlpha <- t(grDevices::col2rgb("gray85", alpha = TRUE))
			colAlpha[4] <- 160
			graphics::polygon(c(stats::density(mcmcChainThreshold, na.rm = TRUE)$x, rev(stats::density(mcmcChainThreshold, na.rm = TRUE)$x)), c(stats::density(mcmcChainThreshold, na.rm = TRUE)$y, rep(0, length(stats::density(mcmcChainThreshold, na.rm = TRUE)$x))), col = grDevices::rgb(colAlpha[1], colAlpha[2], colAlpha[3], colAlpha[4], maxColorValue = 255), border = NA)
			graphics::lines(stats::density(mcmcChainThreshold, na.rm = TRUE), col = "darkgrey", lwd = 2)
			graphics::box(bty = "L", lwd = 1)
			devAskNewPage(ask = FALSE)
		}
		else {
			graphics::par(mfrow = c(1, 1),lwd = 2)
			graphics::hist(mcmcChainThreshold, yaxs = "i", main = "MCMC sample distribution of optimal threshold", col = "gray85", border = "darkgrey", xlab = "Optimal threshold estimate", freq = FALSE, breaks = seq(min(mcmcChainThreshold, na.rm = TRUE), max(mcmcChainThreshold, na.rm = TRUE), length.out = 20))
			colAlpha <- t(grDevices::col2rgb("gray85", alpha = TRUE))
			colAlpha[4] <- 160
			graphics::polygon(c(stats::density(mcmcChainThreshold, na.rm = TRUE)$x, rev(stats::density(mcmcChainThreshold, na.rm = TRUE)$x)), c(stats::density(mcmcChainThreshold, na.rm = TRUE)$y, rep(0, length(stats::density(mcmcChainThreshold, na.rm = TRUE)$x))), col = grDevices::rgb(colAlpha[1], colAlpha[2], colAlpha[3], colAlpha[4], maxColorValue = 255), border = NA)
			graphics::lines(stats::density(mcmcChainThreshold, na.rm = TRUE), col = "darkgrey", lwd = 2)
			graphics::box(bty = "L", lwd = 1)
		}
	}
	paraNamesUserDefined <- vector("list", 2)
	names(paraNamesUserDefined) <- c("Evt", "NoEvt")
	cdfUserDefined <- vector("list", 2)
	names(cdfUserDefined) <- c("Evt", "NoEvt")
	gradientUserDefined <- vector("list", 2)
	names(gradientUserDefined) <- c("Evt", "NoEvt")
	hessianUserDefined <- vector("list", 2)
	names(hessianUserDefined) <- c("Evt", "NoEvt")
	if (any(methods::is(EvtDist, "fitUserDefinedDist0E"), methods::is(NoEvtDist, "fitUserDefinedDist0Eb"))) {
		if (methods::is(EvtDist, "fitUserDefinedDist0E")) {
			paraNamesUserDefined$EvtRef <- EvtDist@paraNames
			cdfUserDefined$EvtRef <- methods::getMethod(cdf, "userDefinedDist0E")
			gradientUserDefined$EvtRef <- methods::getMethod(gradient, "userDefinedDist0E")
			hessianUserDefined$EvtRef <- methods::getMethod(hessian, "userDefinedDist0E")
		}
		if (methods::is(NoEvtDist,"fitUserDefinedDist0Eb")) {
			paraNamesUserDefined$NoEvtRef <- NoEvtDist@paraNames
			cdfUserDefined$NoEvtRef <- methods::getMethod(cdf, "userDefinedDist0Eb")
			gradientUserDefined$NoEvtRef <- methods::getMethod(gradient, "userDefinedDist0Eb")
			hessianUserDefined$NoEvtRef <- methods::getMethod(hessian, "userDefinedDist0Eb")
		}
	}
	res <- methods::new("diagOptThresh", optThresh = mcmcChainThreshold, p = p, r = r, N = EvtDist@n + NoEvtDist@n, xEvt = EvtDist@x, xNoEvt = NoEvtDist@x, lowEvt = lowEvt, mcmcChainEvt = mcmcChainEvt, mcmcChainNoEvt = mcmcChainNoEvt, tabMCMCChain = listTabRes, paraNamesUserDefined = paraNamesUserDefined, cdfUserDefined = cdfUserDefined, gradientUserDefined = gradientUserDefined, hessianUserDefined = hessianUserDefined, percentNA = percentNA)
	if (any(methods::is(EvtDist, "fitUserDefinedDist0E"),methods::is(NoEvtDist, "fitUserDefinedDist0Eb"))) {	
		if (methods::is(EvtDist, "fitUserDefinedDist0E")) {
			methods::removeMethod(cdf, "userDefinedDist0E", where = .optimalThreshold_envEvtDist)
			methods::removeMethod(gradient, "userDefinedDist0E", where = .optimalThreshold_envEvtDist)
			methods::removeMethod(hessian, "userDefinedDist0E", where = .optimalThreshold_envEvtDist)
			methods::removeMethod(samplePosteriorDist, "fitUserDefinedDist0E", where = .optimalThreshold_envEvtDist)
			methods::removeClass("fitUserDefinedDist0E", where = .optimalThreshold_envEvtDist)
			methods::removeClass("userDefinedDist0E", where = .optimalThreshold_envEvtDist)			
		}
		if (methods::is(NoEvtDist, "fitUserDefinedDist0Eb")) {
			methods::removeMethod(cdf, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
			methods::removeMethod(gradient, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
			methods::removeMethod(hessian, "userDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
			methods::removeMethod(samplePosteriorDist, "fitUserDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
			methods::removeClass("fitUserDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
			methods::removeClass("userDefinedDist0Eb", where = .optimalThreshold_envNoEvtDist)
		}
	}
	return(res)
}
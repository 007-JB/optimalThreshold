#################################################################################################
###                                  ClassUnionsDist.R                                        ###
#################################################################################################

#' An S4 union class to represent allowed dist object.
#'
#' This S4 union class describes what types of distribution are supported by the \code{optimalThreshold} package.
#' @details Five theoretical univariate and continuous distributions are supported by the package internally (normal, log-normal, gamma, scaled t, and logistic). However, it is possible for the user to define its own distribution type using the \code{fit} function, and passing the result to one of the two main functions: \code{trtSelThresh}, and \code{diagThresh}.
#' @seealso \code{\link[optimalThreshold]{fit}}, \code{\link[optimalThreshold]{trtSelThresh}}, and \code{\link[optimalThreshold]{diagThresh}}.
#' @name allowedDist-class
#' @aliases allowedDist
#' @exportClass allowedDist
setClassUnion("allowedDist", c("normalDist", "logNormalDist", "gammaDist", "studentDist", "logisticDist"))

#' An S4 union class to represent allowed fitDist objects.
#'
#' This S4 union class describes what types of distribution may be fitted by the user in the \code{optimalThreshold} package.
#' @details Five theoretical types of distribution fit are supported by the package internally (normal, log-normal, gamma, scaled t, and logistic). However, it is possible for the user to fit a personalized distribution using the \code{fit} function, and passing the result to one of the two main functions: \code{trtSelThresh}, and \code{diagThresh}. The 'undefined' type is used to indicate which distribution must be expressed as a function of the three others when estimating the optimal threshold of a treatment selection marker (see References for more details).
#' @section References :
#' Blangero, Y, Rabilloud, M, Ecochard, R, and Subtil, F. A Bayesian method to estimate the optimal threshold of a marker used to select patients' treatment. \emph{Statistical Methods in Medical Research}. 2019.
#' Subtil, F, and Rabilloud, M. A Bayesian method to estimate the optimal threshold of a longitudinal marker. \emph{Biometrical Journal}. 2010.
#' @seealso \code{\link[optimalThreshold]{fit}}, \code{\link[optimalThreshold]{trtSelThresh}}, and \code{\link[optimalThreshold]{diagThresh}}.
#' @name allowedFitDist-class
#' @aliases allowedFitDist
#' @exportClass allowedFitDist
setClassUnion("allowedFitDist", c("fitNormalDist", "fitLogNormalDist", "fitGammaDist", "fitStudentDist", "fitLogisticDist", "undefined", "fitUserDefinedDist"))
#################################################################################################
###                                  ClassCompoundDist.R                                      ###
#################################################################################################

#' An S4 class to represent a compound distribution.
#'
#' This S4 class describes the 'compound' distribution that is fitted to the marker values of patients that developed the event of interest in the reference treatment arm when the type of distribution is set to 'undefined' in the \code{fit} function. 
#' @slot NoEvtRefDist This slot is an object that describes the distribution of the marker for patients that did not develop the event in the reference arm. This object must be an "allowedDist" class object.
#' @slot EvtInnovDist This slot is an object that describes the distribution of the marker for patients that developed the event in the innovative arm. This object must be an "allowedDist" class object.
#' @slot NoEvtInnovDist This slot is an object that describes the distribution of the marker for patients that did not develop the event in the innovative arm. This object must be an "allowedDist" class object.
#' @slot r0 Mean risk of event occurrence in the reference arm. Numeric argument.
#' @slot r1 Mean risk of event occurrence in the innovative arm. Numeric argument.
#' @details You never have to create this class manually. This class is created internally when an undefined distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit an undefined distribution.
#' @name compoundEvtRefDist-class
#' @aliases compoundEvtRefDist
#' @exportClass compoundEvtRefDist
setClass("compoundEvtRefDist", slots = c(NoEvtRefDist = "allowedDist", EvtInnovDist = "allowedDist", NoEvtInnovDist = "allowedDist", r0 = "numeric", r1 = "numeric"))

#' An S4 class to represent a compound distribution.
#'
#' This S4 class describes the 'compound' distribution that is fitted to the marker values of patients that did not develop the event of interest in the reference treatment arm when the type of distribution is set to 'undefined' in the \code{fit} function. 
#' @slot EvtRefDist This slot is an object that describes the distribution of the marker for patients that developed the event in the reference arm. This object must be an "allowedDist" class object.
#' @slot EvtInnovDist This slot is an object that describes the distribution of the marker for patients that developed the event in the innovative arm. This object must be an "allowedDist" class object.
#' @slot NoEvtInnovDist This slot is an object that describes the distribution of the marker for patients that did not develop the event in the innovative arm. This object must be an "allowedDist" class object.
#' @slot r0 Mean risk of event occurrence in the reference arm. Numeric argument.
#' @slot r1 Mean risk of event occurrence in the innovative arm. Numeric argument.
#' @details You never have to create this class manually. This class is created internally when an undefined distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit an undefined distribution.
#' @name compoundNoEvtRefDist-class
#' @aliases compoundNoEvtRefDist
#' @exportClass compoundNoEvtRefDist
setClass("compoundNoEvtRefDist", slots = c(EvtRefDist = "allowedDist", EvtInnovDist = "allowedDist", NoEvtInnovDist = "allowedDist", r0 = "numeric", r1 = "numeric"))

#' An S4 class to represent a compound distribution.
#'
#' This S4 class describes the 'compound' distribution that is fitted to the marker values of patients that developed the event of interest in the innovative treatment arm when the type of distribution is set to 'undefined' in the \code{fit} function. 
#' @slot EvtRefDist This slot is an object that describes the distribution of the marker for patients that developed the event in the reference arm. This object must be an "allowedDist" class object.
#' @slot NoEvtRefDist This slot is an object that describes the distribution of the marker for patients that did not develop the event in the reference arm. This object must be an "allowedDist" class object.
#' @slot NoEvtInnovDist This slot is an object that describes the distribution of the marker for patients that did not develop the event in the innovative arm. This object must be an "allowedDist" class object.
#' @slot r0 Mean risk of event occurrence in the reference arm. Numeric argument.
#' @slot r1 Mean risk of event occurrence in the innovative arm. Numeric argument.
#' @details You never have to create this class manually. This class is created internally when an undefined distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit an undefined distribution.
#' @name compoundEvtInnovDist-class
#' @aliases compoundEvtInnovDist
#' @exportClass compoundEvtInnovDist
setClass("compoundEvtInnovDist", slots = c(EvtRefDist = "allowedDist", NoEvtRefDist = "allowedDist", NoEvtInnovDist = "allowedDist", r0 = "numeric", r1 = "numeric"))

#' An S4 class to represent a compound distribution.
#'
#' This S4 class describes the 'compound' distribution that is fitted to the marker values of patients that did not develop the event of interest in the innovative treatment arm when the type of distribution is set to 'undefined' in the \code{fit} function. 
#' @slot EvtRefDist This slot is an object that describes the distribution of the marker for patients that developed the event in the reference arm. This object must be an "allowedDist" class object.
#' @slot NoEvtRefDist This slot is an object that describes the distribution of the marker for patients that did not develop the event in the reference arm. This object must be an "allowedDist" class object.
#' @slot EvtInnovDist This slot is an object that describes the distribution of the marker for patients that developed the event in the innovative arm. This object must be an "allowedDist" class object.
#' @slot r0 Mean risk of event occurrence in the reference arm. Numeric argument.
#' @slot r1 Mean risk of event occurrence in the innovative arm. Numeric argument.
#' @details You never have to create this class manually. This class is created internally when an undefined distribution is fitted on the marker values.
#' @seealso \code{\link[optimalThreshold]{fit}} for more details on how to fit an undefined distribution.
#' @name compoundNoEvtInnovDist-class
#' @aliases compoundNoEvtInnovDist
#' @exportClass compoundNoEvtInnovDist
setClass("compoundNoEvtInnovDist", slots = c(EvtRefDist = "allowedDist", NoEvtRefDist = "allowedDist", EvtInnovDist = "allowedDist", r0 = "numeric", r1 = "numeric"))

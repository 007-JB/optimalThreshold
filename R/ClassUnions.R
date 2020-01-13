#################################################################################################
###                                  ClassUnions.R                                            ###
#################################################################################################
setOldClass("mcmc.list") 

#' An S4 union class to merge mcmc.list and NULL types of object.
#'
#' @name mcmc.listOrNull-class
#' @aliases mcmc.listOrNull
#' @exportClass mcmc.listOrNull
setClassUnion("mcmc.listOrNull", c("mcmc.list", "NULL"))
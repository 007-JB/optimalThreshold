#################################################################################################
###                                  ARS.R                                                    ###
#################################################################################################

fInArs <- function(shape, xVec, scale) {
	n <- length(xVec)
	- n * log(gamma(shape)) - n * shape * log(scale) + (shape - 1) * sum(log(xVec))
}

fPrimaInArs <- function(shape, xVec, scale) {
	n <- length(xVec)
	- n * digamma(shape) - n * log(scale) + sum(log(xVec))
}

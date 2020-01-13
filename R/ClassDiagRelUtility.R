#################################################################################################
###                                  ClassDiagRelUtility.R                                    ###
#################################################################################################
#' An S4 class to sum up the results from the decisionCurve methods.
#'
#' @slot U This slot is a matrix of the marker-based utility according to \code{r}.
#' @slot UNoTreat This slot is a matrix of the utility of the 'No Treat' strategy according to \code{r}.
#' @slot UTreatAll This slot is a matrix of the utility of the 'Treat All' strategy according to \code{r}.
#' @slot r Risk threshold preference. Numeric argument.
#' @details You never have to create this class manually. This class is created internally when the \code{decisionCurve} method is applied to a 'diagOptThresh' object.
#' @seealso \code{\link[optimalThreshold]{decisionCurve}} for more details on how to plot the decision curves.
#' @name diagRelUtility-class
#' @aliases diagRelUtility
#' @exportClass diagRelUtility
setClass("diagRelUtility", slots = c(U = "matrix", UNoTreat = "matrix", UTreatAll = "matrix", r = "numeric"))

#' Plot the decision curves of a diagnostic marker
#'
#' @name plot.diagRelUtility
#' @aliases plot,diagRelUtility-method
#' @param x a \code{diagRelUtility} object.
#' @param y unused parameter.
#' @param main an overall title for the plot.
#' @param lty the line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank", "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses 'invisible lines' (i.e., does not draw them). 
#' @param lwd the line width, a \emph{positive} number, defaulting to 1. The interpretation is device-specific, and some devices do not implement line widths less than one. (See the help on the device for details of the interpretation). 
#' @param xlim the x limits of the plot. 
#' @param ylim the x limits of the plot. 
#' @param ylab a label for the y axis.
#' @param xlab a label for the x axis
#' @param col.U color of the utility curve for the marker-based strategy. 
#' @param col.UNoTreat color of the utility curve for the "No Treat" strategy.
#' @param col.UTreatAll color of the utility curve for the "Treat All" strategy.
#' @param ... other graphical parameters.
#' @return None
#' @exportMethod plot
setMethod("plot", "diagRelUtility",
	function(x, y, main = "Decision curves", lty = 1, lwd = 1, xlim = range(x@r), ylim = c(min(x@U, x@UNoTreat, x@UTreatAll), max(x@U, x@UNoTreat, x@UTreatAll)), ylab = "Expected benefit", xlab = "r", col.U = "black", col.UNoTreat = "blue", col.UTreatAll = "green", ...) {
		opar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(opar))
		graphics::par(mfrow = c(1, 1), lwd = 1, mar = c(4.1, 4.1, 2, 0.5))
		Uc <- apply(x@U, 2, stats::median, na.rm = TRUE)
		UNoTreat <- apply(x@UNoTreat, 1, stats::median, na.rm = TRUE)
		UTreatAll <- apply(x@UTreatAll, 1, stats::median, na.rm = TRUE)
		graphics::plot(Uc ~ x@r, ylim = ylim, col = col.U, type = "l", xlab = xlab, ylab = ylab, main = main, xlim = xlim, lty = lty, lwd = lwd, ...)
		graphics::lines(UNoTreat ~ x@r, col = col.UNoTreat, lty = lty, lwd = lwd, ...)
		graphics::lines(UTreatAll ~ x@r, col = col.UTreatAll, lty = lty, lwd = lwd, ...)
		graphics::legend("topright", col = c(col.U, col.UNoTreat, col.UTreatAll), lty = c(lty, lty, lty), legend = c("Marker under study", "No Treat", "Treat All"))
	}
)
#################################################################################################
###                                  ClassTrtSelRelUtility.R                                  ###
#################################################################################################
#' An S4 class to the results from the decisionCurve methods.
#'
#' @slot RU This slot is a matrix of the relative utility according to the r ratios.
#' @slot U This slot is a matrix of the marker-based utility according to the r ratios.
#' @slot UT0 This slot is a matrix of the reference treatment utility according to the r ratios.
#' @slot UT1 This slot is a matrix of the innovative treatment utility according to the r ratios.
#' @slot Up This slot is a matrix of perfect marker utility according to the r ratios.
#' @slot r Ratio of treatment/event costs. Numeric argument.
#' @details You never have to create this class manually. This class is created internally when the \code{decisionCurve} method is applied to an 'optThresh' object.
#' @seealso \code{\link[optimalThreshold]{decisionCurve}} for more details on how to plot the decision curves.
#' @name trtSelRelUtility-class
#' @aliases trtSelRelUtility
#' @exportClass trtSelRelUtility
setClass("trtSelRelUtility", slots = c(RU = "matrix", U = "matrix", UT0 = "matrix", UT1 = "matrix", Up = "matrix", r = "numeric"))

#' Plot the decision curves of a treatment selection marker
#'
#' @name plot.trtSelRelUtility
#' @aliases plot,trtSelRelUtility-method
#' @param x a \code{trtSelRelUtility} object.
#' @param y unused parameter.
#' @param which indicates which graph should be plotted. Default is both graphs.
#' @param alpha alpha risk for the confidence intervals.
#' @param conf.int a logical value indicating whether the confidence intervals should be plotted for the relative utility curve.
#' @param main1 an overall title for the first plot.
#' @param main2 an overall title for the second plot.
#' @param lty the line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank", "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses 'invisible lines' (i.e., does not draw them).
#' @param lwd the line width, a \emph{positive} number, defaulting to 1. The interpretation is device-specific, and some devices do not implement line widths less than one. (See the help on the device for details of the interpretation). 
#' @param xlim the x limits of the plots. 
#' @param ylim1 the y limits of the first plot. 
#' @param ylim2 the y limits of the second plot.
#' @param ylab1 a label for the y axis of the first plot.
#' @param ylab2 a label for the y axis of the second plot.
#' @param xlab a label for the x axis of the plots.
#' @param col.U color of the utility curve for the marker-based strategy. 
#' @param col.Up color of the utility curve for the perfect marker-based strategy. 
#' @param col.UT0 color of the utility curve for the "Treat All with the reference treatment" strategy. 
#' @param col.UT1 color of the utility curve for the "Treat All with the innovative treatment" strategy. 
#' @param col.RU color of the relative utility curve. 
#' @param col.conf.int color of the confidence intervals. 
#' @param add a logical value indicating whether the relative utility curve should superimpose with an existing graph. Only works when \code{which} = 2.
#' @param legend1 a logical value indicating whether a legend should be added to the first plot.
#' @param ... other graphical parameters.
#' @return None
#' @exportMethod plot
setMethod("plot", "trtSelRelUtility",
	function(x, y, which = c(1, 2), alpha = 0.05, conf.int = TRUE, main1 = "Decision curves (unscaled)", main2 = "Decision curve (scaled)", lty = 1, lwd = 1, xlim = range(x@r), ylim1 = c(min(x@U, x@UT0, x@UT1, x@Up), max(x@U, x@UT0, x@UT1, x@Up)), ylim2 = c(0, 1), ylab1 = "Utility", ylab2 = "Relative utility", xlab = "r ratio", col.U = "black", col.Up = "red", col.UT0 = "blue", col.UT1 = "green", col.RU = "black", col.conf.int = "black", add = FALSE, legend1 = TRUE, ...) {
		opar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(opar))
		if (!(is.numeric(which)) || any(which > 2) || any(which < 1)) stop("which must be in 1:2.")
		if (alpha < 0 || alpha > 1) stop("alpha must be in [0;1].")
		if (length(which) == 2) graphics::par(mfrow = c(1, 2), lwd = 1, mar = c(4.1, 4.1, 2, 0.5))
		else graphics::par(mfrow=c(1, 1), lwd = 1, mar = c(4.1, 4.1, 2, 0.5))
		show <- rep(FALSE, 2)
		show[which] <- TRUE
		if (show[1]) {
			Uc <- apply(x@U, 2, stats::median, na.rm = TRUE)
			Up <- apply(x@Up, 2, stats::median, na.rm = TRUE)
			if (length(x@UT0) == 1) UT0 <- rep(x@UT0, length(x@r))
			else UT0 <- apply(x@UT0, 2, stats::median, na.rm = TRUE)
			if (length(x@UT1) == 1) UT1 <- rep(x@UT1, length(x@r))
			else UT1 <- apply(x@UT1, 2, stats::median, na.rm = TRUE)
			graphics::plot(Uc ~ x@r, ylim = ylim1, col = col.U, type = "l", xlab = xlab, ylab = ylab1, main = main1, xlim = xlim, lty = lty, lwd = lwd, ...)
			graphics::lines(Up ~ x@r, col = col.Up, lty = lty, lwd = lwd, ...)
			graphics::lines(UT0 ~ x@r, col = col.UT0, lty = lty, lwd = lwd, ...)
			graphics::lines(UT1 ~ x@r, col = col.UT1, lty = lty, lwd = lwd, ...)
			graphics::legend("topright", col = c(col.Up, col.U, col.UT0, col.UT1), lty = c(lty, lty, lty, lty), legend = c("Perfect marker", "Marker under study", "All under reference", "All under innovative"))
		}
		if (show[2]) {
			if (conf.int) {
				resRU <- apply(x@RU, 2, function(z) stats::quantile(z, probs = c(alpha / 2, 0.5, 1 - alpha / 2), na.rm = TRUE))
				if (add) {
					graphics::arrows(x0 = x@r, y0 = resRU[1, ], x1 = x@r, y1 = resRU[3,], length = 0.05, angle = 90, code = 3, col = col.conf.int, lwd = lwd)
					graphics::lines(resRU[2, ] ~ x@r, col = col.RU, lty = lty, lwd = lwd, ...)
				}
				else {
					graphics::plot(resRU[2,] ~ x@r, type = "n", main = main2, ylim = ylim2, ylab = ylab2, xlab = xlab, xlim = xlim, ...)
					graphics::abline(h = 0, lty = 2, col = "lightgrey")
					graphics::abline(h = 1, lty = 2, col = "lightgrey")
					for(i in seq(0.1, 0.9, by = 0.1)) graphics::abline(h = i, col = "lightgrey")
					graphics::arrows(x0 = x@r, y0 = resRU[1, ], x1 = x@r, y1 = resRU[3, ], length = 0.05, angle = 90, code = 3, col = col.conf.int, lwd = lwd)
					graphics::lines(resRU[2, ] ~ x@r, col = col.RU, lty = lty, lwd = lwd, ...)
				}
			}
			else {
				resRU <- apply(x@RU, 2, function(z) stats::median(z, na.rm = TRUE))
				if (add) {
					graphics::lines(resRU ~ x@r, col = col.RU, lty = lty, lwd = lwd, ...)
				}
				else {
					graphics::plot(resRU ~ x@r, type = "n", main = main2, ylim = ylim2, ylab = ylab2, xlab = xlab, xlim = xlim, ...)
					graphics::abline(h = 0, lty = 2, col = "lightgrey")
					graphics::abline(h = 1, lty = 2, col = "lightgrey")
					for(i in seq(0.1, 0.9, by = 0.1)) graphics::abline(h = i, col = "lightgrey")
					graphics::lines(resRU ~ x@r, col = col.RU, lty = lty, lwd = lwd, ...)
				}
			}
		}
	}
)
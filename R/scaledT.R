#################################################################################################
###                                  scaledT.R                                                ###
#################################################################################################
pt.scaled <- function(q, df, mean = 0, sd = 1, ncp, log.p = FALSE) {
    stats::pt((q - mean) / sd, df, ncp = ncp, log.p = log.p)
}

dt.scaled <- function(x, df, mean = 0, sd = 1, ncp, log = FALSE) {
    if (!log) stats::dt((x - mean) / sd, df, ncp = ncp, log = FALSE) / sd
    else stats::dt((x - mean) / sd, df, ncp = ncp, log = TRUE) - log(sd)
}

rt.scaled <- function(n, df, mean = 0, sd = 1, ncp) {
    mean + sd * stats::rt(n, df, ncp = ncp)
}






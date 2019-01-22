#' Compute summary statistics of per-slice FD measurements for one individual
#'
#' Function based on
#' https://github.com/UK-Digital-Heart-Project/AutoFD/blob/master/pft_JC_FDStatistics.m
#'
#' @param data [vector] of FD measurements at [NrSlices] for one individual
#' @param discard [logical] Set True if first and last slice should be exclueded
#' @param NaN.val [character] Name/vector of names for possible non-numeric
#' values for slice measurements, for instance explanation for technical failure
#' of FD measurement.
#' @return named [vactor] with i) the number of slices used to compute the
#' summary statistics (SlicesUsed), i) the MeanGlobalFD, iii) the MeanApicalFD,
#' iv) the MaxApicalFD, v) the MeanBasalFD and the vi) MaxBasalFD.

summaryStatistics <- function(data, discard=FALSE,
                      NaN.val=c("Sparse myocardium", "Meagre blood pool",
                                "FD measure failed")) {
    data <- as.matrix(data)
    # index of last non-NA value
    n <- which(!is.na(data))[length(which(!is.na(data)))]
    q <- floor(n/2)
    r <- n %% 2

    if (!discard) {
        SlicesUsed <- n
        LG <- 1           # Lower global index
        UG <- n           # Upper global index
        LB <- 1           # Lower basal index
        UB <- q           # Upper basal index
        LA <- 1 + q + r   # Lower apical index
        UA <- n           # Upper apical index
    } else {
        SlicesUsed <- n - 2
        LG <- 2           # Lower global index
        UG <- n - 1       # Upper global index
        LB <- 2           # Lower basal index
        UB <- q           # Upper basal index
        LA <- 1 + q + r   # Lower apical index
        UA <- n - 1       # Upper apical index
    }
    data[data %in% NaN.val] <- NA
    data <- as.numeric(data)

    # Trim the data if necessary, but in any case, assign the global, basal and
    # apical arrays
    GlobalFD <- data[LG:UG]
    BasalFD  <- data[LB:UB]
    ApicalFD <- data[LA:UA]

    MeanGlobalFD <- ifelse(all(is.na(GlobalFD)), NA, mean(GlobalFD, na.rm=TRUE))
    MeanApicalFD <- ifelse(all(is.na(ApicalFD)), NA, mean(ApicalFD, na.rm=TRUE))
    MaxApicalFD <- ifelse(all(is.na(ApicalFD)), NA, max(ApicalFD, na.rm=TRUE))
    MeanBasalFD <- ifelse(all(is.na(BasalFD)), NA, mean(BasalFD, na.rm=TRUE))
    MaxBasalFD <- ifelse(all(is.na(BasalFD)), NA, max(BasalFD, na.rm=TRUE))

    # Calculate the statistics
    Statistics <- c(SlicesUsed=SlicesUsed, MeanGlobalFD=MeanGlobalFD,
                    MeanApicalFD=MeanApicalFD, MaxApicalFD=MaxApicalFD,
                    MeanBasalFD=MeanBasalFD, MaxBasalFD=MaxBasalFD)

    return(Statistics)
}

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
#' @param sections [character] determines how many sections to compute summary
#' statistics for; one of BA (basal and apical summary statistics) and BMA
#' (basal, mid and apical section summary statistics).
#' @return named [vector] with i) the number of slices used to compute the
#' summary statistics (SlicesUsed); if sections == BA: ii) the MeanGlobalFD,
#' iii) the MeanBasalFD, iv) the MaxBasalFD, v) the MeanApicalFD,
#' vi) the MaxApicalFD; if sections == MBA: ii) the MeanGlobalFD,
#' iii) the MeanBasalFD, iv) the MaxBasalFD, v) the MeanMidFD, vi) the MaxMidFD,
#' vii) the MeanApicalFD, viii) the MaxApicalFD.

summaryStatistics <- function(data, discard=FALSE, sections="BA",
                      NaN.val=c("Sparse myocardium", "Meagre blood pool",
                                "FD measure failed")) {
    data <- as.matrix(data)
    # index of last non-NA value
    n <- which(!is.na(data))[length(which(!is.na(data)))]
    p <- floor(n/3)
    q <- floor(n/2)
    r <- n %% 2
    s <- n %% 3

    if (!discard) {
        SlicesUsed <- n
        LG <- 1           # Lower global index
        UG <- n           # Upper global index
        if (sections == "BA") {
            LB <- 1           # Lower basal index
            UB <- q           # Upper basal index
            LA <- 1 + q + r   # Lower apical index
            UA <- n           # Upper apical index
        } else if (sections == "BMA") {
            if (s == 2) {
                LB <- 1           # Lower basal index
                UB <- 1 + p       # Upper basal index
                LM <- 2 + p       # Lower mid index
                UM <- 2*p + 1     # Upper mid index
                LA <- 2*p + 2     # Lower apical index
                UA <- n           # Upper apical index
            } else if (s == 1) {
                LB <- 1
                UB <- p
                LM <- 1 + p
                UM <- 2*p
                LA <- 2*p + 1
                UA <- n
            } else {
                LB <- 1
                UB <- p
                LM <- 1 + p
                UM <- 2*p
                LA <- 2*p + 1
                UA <- n
            }
        } else {
            stop("Unknown section division chosen, Has to be one of BA
                 (Basal-Apical) or BMA (Basal-Mid-Apical)")
        }
    } else {
        SlicesUsed <- n - 2
        LG <- 2           # Lower global index
        UG <- n - 1       # Upper global index
        if (sections == "BA") {
            LB <- 2           # Lower basal index
            UB <- q           # Upper basal index
            LA <- 1 + q + r   # Lower apical index
            UA <- n - 1       # Upper apical index
        } else if (sections == "BMA") {
            s <- floor((n-2)/3)
            LB <- 2           # Lower basal index
            UB <- 2 + s       # Upper basal index
            LM <- 2 + s + 1   # Lower mid index
            UM <- 2 + 2*s     # Upper mid index
            LA <- 2+ 2*s + 1  # Lower apical index
            UA <- n - 1       # Upper apical index
        } else {
            stop("Unknown section division chosen, Has to be one of BA
                 (Basal-Apical) or BMA (Basal-Mid-Apical)")
        }
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

    if (sections == "BA") {
        Statistics <- c(SlicesUsed=SlicesUsed, MeanGlobalFD=MeanGlobalFD,
                        MeanBasalFD=MeanBasalFD, MaxBasalFD=MaxBasalFD,
                        MeanApicalFD=MeanApicalFD, MaxApicalFD=MaxApicalFD)
    } else {
        MidFD <- data[LM:UM]
        MeanMidFD <- ifelse(all(is.na(MidFD)), NA, mean(MidFD, na.rm=TRUE))
        MaxMidFD <- ifelse(all(is.na(MidFD)), NA, max(MidFD, na.rm=TRUE))
        Statistics <- c(SlicesUsed=SlicesUsed, MeanGlobalFD=MeanGlobalFD,
                        MeanBasalFD=MeanBasalFD, MaxBasalFD=MaxBasalFD,
                        MeanMidFD=MeanMidFD, MaxMidFD=MaxMidFD,
                        MeanApicalFD=MeanApicalFD, MaxApicalFD=MaxApicalFD)
    }

    return(Statistics)
}

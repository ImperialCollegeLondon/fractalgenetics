#' Remove related samples while keeping maximum number of samples in cohort
#'
#' @param ID_study [N] vector with N samples IDs of study cohort
#' @param relatedness relatedness data.frame as obtained from ukbgene -rel
#' containing pairwise relatedness estimates for individuals in the 500 UKB
#' genotypes. Has to have original header i.e. ID1 ID2 HetHet IBS0 Kinship
#' @return named list with related2filter containing IDs of related samples in
#' the study cohort that should be removed and related2keep, containing the IDs
#' of samples in the study cohort that have related individuals in UKB but those
#' are not in the current study - these IDs are just provided for information
#' and nothing needs to be done.

smartRelatednessFiter <- function(ID_study, relatedness) {
    relatedness_study <-
        relatedness[unique(union(which(relatedness$ID1 %in% ID_study),
                             which(relatedness$ID2 %in% ID_study))),]

    # clear individuals whose relatives is not part of ID_study
    id_class <- rowSums(cbind(relatedness_study$ID1 %in% ID_study,
            relatedness_study$ID2 %in% ID_study))
    both <- relatedness_study[id_class == 2,]

    single_tmp <- relatedness_study[id_class == 1,]
    single_in_both <- rowSums(cbind(single_tmp$ID1 %in% c(both$ID1, both$ID2),
            single_tmp$ID2 %in% c(both$ID1, both$ID2)))
    single <- single_tmp[single_in_both == 0,]

    # check
    id_both <- unique(c(both$ID1, both$ID2))
    id_single <- unique(c(single$ID1, single$ID2))
    id_single_in_both <- unique(c(single_tmp$ID1[single_in_both == 1],
                            single_tmp$ID2[single_in_both == 1]))
    id_not_single_nor_both <-
        id_single_in_both[!(id_single_in_both %in%
                        unique(c(both$ID1, both$ID2,single$ID1, single$ID2)))]
    id_all <- unique(c(relatedness_study$ID1, relatedness_study$ID2))

    id_combined <- length(id_both) + length(id_single) +
        length(id_not_single_nor_both)
    if (id_combined != length(id_all)) {
        stop("Not all IDs accounted for")
    }
    tokeep <- ID_study[ID_study %in% id_single]
    tofilter <- ID_study[ID_study %in% id_both]
    return(list(related2filter=tofilter, related2keep=tokeep))
}

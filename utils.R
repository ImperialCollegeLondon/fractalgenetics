#' Remove related samples while keeping maximum number of samples in cohort
#'
#' @param IDs [N] vector with N samples IDs of study cohort
#' @param relatedness relatedness data.frame as obtained from ukbgene -rel
#' containing pairwise relatedness estimates for individuals in the 500 UKB
#' genotypes. Has to have original header i.e. ID1 ID2 HetHet IBS0 Kinship
#' @return named list with related2filter, related2decide and related2keep.
#' related2filter contains IDs of related samples in the study cohort that
#' should be removed. related2decide contains a dataframe with pairs of related
#' that are only related to one another. One of each pair should be removed -
#' if no other reason is given, simply choose on column as samples to keep and
#' the other as samples to filter. related2keep, contains the IDs
#' of samples in the study cohort that have related individuals in UKB but those
#' are not in the current study - these IDs are just provided for information
#' and nothing needs to be done.

smartRelatednessFilter <- function(IDs, relatedness) {
    relatedness_study <-
        relatedness[unique(union(which(relatedness$ID1 %in% ID),
                             which(relatedness$ID2 %in% ID))),]
    ID_study <- ID[ID %in% c(relatedness_study$ID1, relatedness_study$ID2)]
    # clear individuals whose relatives is not part of ID_study
    id_class <- rowSums(cbind(relatedness_study$ID1 %in% ID_study,
            relatedness_study$ID2 %in% ID_study))
    both <- relatedness_study[id_class == 2,]

    single_tmp <- relatedness_study[id_class == 1,]
    single_in_both <- rowSums(cbind(single_tmp$ID1 %in% c(both$ID1, both$ID2),
            single_tmp$ID2 %in% c(both$ID1, both$ID2)))

    # Those are the ones to keep, as relative is not in study cohort
    single <- single_tmp[single_in_both == 0,]

    # check that ID numbers match up
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

    # Further selection of samples with relatives in cohort
    duplicate_IDs <- c(both$ID1, both$ID2)[duplicated(c(both$ID1, both$ID2))]

    # From these, either sample of a pair can be kept in the analysis
    singlets <- both[rowSums(cbind(both$ID1 %in% duplicate_IDs,
                                   both$ID2 %in% duplicate_IDs)) == 0,]

    # Check complex relative scenarios: investigate duplicate IDs
    duplicates <- both[rowSums(cbind(both$ID1 %in% duplicate_IDs,
                              both$ID2 %in% duplicate_IDs)) != 0,]

    # Get all related samples per individual
    rel_per_ID <- lapply(duplicate_IDs, function(x) {
        tmp <- both[rowSums(cbind(both$ID1 %in% x,  both$ID2 %in% x)) != 0,1:2]
        rel <- unique(unlist(tmp))
        return(rel)
    })

    # List of individuals that share at least 2 relatives
    rel_ge_2  <- lapply(seq_along(rel_per_ID), function(x) {
        tmp <- sapply(seq_along(rel_per_ID), function(y) {
            if (x >= y) {
                return(FALSE)
            }
            if (length(intersect(rel_per_ID[[x]],
                                 rel_per_ID[[y]])) == length(rel_per_ID[[x]])) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        })
        return(unique(unlist(rel_per_ID[which(tmp)])))
    })
    rel_ge_2 <- unique(unlist(rel_ge_2))

    # List of trios
    rel_trios_rm <- duplicate_IDs[!(duplicate_IDs %in% rel_ge_2)]
    rel_trios <- unlist(both[rowSums(cbind(both$ID1 %in% rel_trios_rm,
                              both$ID2 %in% rel_trios_rm)) != 0,1:2])
    tokeep_trios <- rel_trios[!(rel_trios %in% rel_trios_rm)]

    # Check that all IDs with relatives in cohort are accounted for
    id_singlets <- unique(c(singlets$ID1, singlets$ID2))

    id_both_combined <-  length(rel_ge_2) + length(id_singlets) +
        length(rel_trios_rm) + length(tokeep_trios)

    if (id_both_combined != length(id_both)) {
        stop("Not all IDs of samples with relatives in study accounted for")
    }

    tokeep <- ID_study[ID_study %in% c(id_single, tokeep_trios)]
    tofilter <- ID_study[ID_study %in% c(rel_ge_2, rel_trios_rm)]
    if(length(tokeep) + length(tofilter) + length(id_singlets)
       != length(ID_study)) {
        stop("Not all IDs returned")
    }

    return(list(related2filter=tofilter, related2decide=singlets,
        related2keep=tokeep))
}



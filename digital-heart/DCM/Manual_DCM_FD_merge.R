# From slack conversation with Antonio (september 18, 2019):
# If there are repeats I would pick the ED_edited version.
# All done manually rather than attempting to correct auto segmentations

both <- read.table("~/data/digital-heart/phenotype/FD/DCM.csv",
                      sep=",", stringsAsFactors = FALSE, header=TRUE)

ed_only <- read.table("~/data/digital-heart/phenotype/FD/DCM_ED_edited.csv",
                      sep=",", stringsAsFactors = FALSE, header=TRUE)

collection <- read.table("~/data/digital-heart/phenotype/FD/20191002_DCM_FD_Collection.csv",
                      sep=",", stringsAsFactors = FALSE, header=TRUE, quote="")

both$Folder <- gsub("([A-Z0-9]*)_.*", "\\1", both$Folder)
ed_only$Folder <- gsub("([A-Z0-9]*)_.*", "\\1", ed_only$Folder)

duplicated_both <- data.frame(ID=both$Folder[duplicated(both$Folder)],
                              Type='duplicated in DCM.csv',
                              stringsAsFactors = FALSE)
duplicated_across <- data.frame(ID=both$Folder[both$Folder %in% ed_only$Folder],
                                Type='duplicated across DCM.csv and DCM_ED_edited.csv',
                                stringsAsFactors = FALSE)
empty_both <- data.frame(ID=both$Folder[apply(both[,10:29], 1,
                                function(x) all(is.na(x) | x == "  "))],
                         Type='empty rows in DCM.csv',
                         stringsAsFactors = FALSE)
exclude_both <- rbind(duplicated_both, duplicated_across, empty_both)

# no duplicates in ed_only
duplicated_ed_only <- ed_only$Folder[duplicated(ed_only$Folder)]
empty_ed_only <- ed_only$Folder[apply(ed_only[,10:29], 1,
                                      function(x)  all(is.na(x) | x == "  "))]

# remove duplicate and empty IDs and combine both data.frames
combined <- rbind(both[!both$Folder %in% exclude_both$ID,], ed_only)

write.table(combined, file="~/data/digital-heart/phenotype/FD/DCM_FD_all.csv",
            sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

# merge exclude IDs with original data
exclude <- merge(exclude_both, both, by=1, all.y=FALSE)
write.table(exclude,
            file="~/data/digital-heart/phenotype/FD/DCM_FD_exclude.csv",
            sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

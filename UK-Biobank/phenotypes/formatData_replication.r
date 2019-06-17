#################
## libraries ####
#################
optparse <- modules::import_package('optparse')


#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-u", "--ukbdir"), action="store", dest="rawdir",
               type="character", help="Path to ukb directory with decrypted ukb
               key.html file [default: %default].", default=NULL),
    make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- optparse$parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$rawdir <- "~/data/ukbb/ukb-hrt/rawdata"
}
################
## analysis ####
################

## LV volume measurements ####
lvv <- data.table::fread(paste(args$rawdir,
                         "/190402_clinical_measures_vent_26k_qced_for_gwas.csv",
                         sep=""),
                         data.table=FALSE, stringsAsFactors=FALSE,
                         na.strings=c("NA", "NaN"))
hr <- data.table::fread(paste(args$rawdir, "/20181124_HR_from_CMR.csv", sep=""),
                        data.table=FALSE, stringsAsFactors=FALSE)

lvv <- lvv[,c(1, 15:24)]
colnames(lvv) <- gsub("(.*) \\(.*", "\\1", colnames(lvv))
colnames(lvv)[1] <- "ID"

lvv$SV <- lvv$LVEDV - lvv$LVESV
rownames(lvv) <- lvv[, 1]

lvv <- dplyr::select(lvv, ID, LVEDV, LVESV, LVEF, LVM, SV)
write.table(lvv, paste(args$rawdir, "/190402_cardiac_phenotypes.csv", sep=""),
            quote=FALSE, col.names=TRUE, row.names=FALSE)

###############################
### Libraries and functions ###
###############################
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('optparse')

###########
## data ###
###########

## command line arguments ####
option_list <- list(
    optparse$make_option(c("-o", "--outdir"), action="store", dest="outdir",
               type="character", help="Path to output directory
               [default: %default].", default=NULL),
    optparse$make_option(c("-p", "--pheno"), action="store", dest="phenofile",
               type="character", help="Path to replication bgenie fd phenotype
               file; [default:%default].", default=NULL),
    optparse$make_option(c("-cov", "--covariates"), action="store",
               dest="covfile", type="character", help="Path to replication
               bgenie fd covariates file; [default:%default].", default=NULL),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out [default: %default]."),
    optparse$make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the
               script [default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    directory <- "homes/hannah/data/ukbb/ukb-hrt"
    args$outdir <-"/homes/hannah/data/ukbb/ukb-hrt/grs"
    args$cohort <- "ukb"
    args$phenofile <- file.path(directory, "phenotypes",
                                "190402_fractal_dimension_26k",
                                "FD_summary_bgenie.txt")
    args$covfile <- file.path(directory, "phenotypes",
                              "190402_fractal_dimension_26k",
                              "FD_covariates_bgenie.txt")
    args$verbose <- TRUE
}

## Read files ####
if (args$verbose) message("Read files with phenotypes and covariates")
pheno <- data.table::fread(args$phenofile, data.table=FALSE,
                               stringsAsFactors=FALSE)
covs <- data.table::fread(args$covfile, data.table=FALSE, stringsAsFactors=FALSE)

###############
## analysis ###
###############

if (verbose) message("Format results for grs with PRSice")
## ukb discovery cohort
pheno$IID <- samples$ID_1
pheno <- pheno[,c(ncol(pheno),1:(ncol(pheno)-1))]

covs$IID <- samples$ID_1
covs <- covs[,c(ncol(covs), 1:(ncol(covs)-1))]
if (args$cohort == "ukb") {
    colnames(covs)[1:6] <- c("IID", "sex", "age", "weight", "bmi", "height")
}
if (args$cohort == "dh") {
    covs <- covs[,c(ncol(covs), 1:3, 5, 4, 6:(ncol(covs)-1))]
}

## write phenotype and summary stats file: one per pheno/region
write.table(pheno, file=paste(args$outdir, "/prsice_pheno_replication_",
                              args$cohort, ".txt", sep=""),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)
                })

## write covariates for replication
write.table(covs, file=file.path(args$outdir,
                                 paste("prsice_covariates_replication_",
                                       args$cohort, "ukb.txt", sep="")),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

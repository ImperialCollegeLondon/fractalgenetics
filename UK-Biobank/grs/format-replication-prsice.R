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
    optparse$make_option(c("-coh", "--cohort"), action="store",
               dest="cohort", type="character", help="Name of replication
               cohort; [default:%default].", default=NULL),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out [default: %default]."),
    optparse$make_option(c("-s", "--samples"), action="store",
               dest="samplesfile", type="character", help="Path to fam samples
               file; [default:%default].", default=NULL),
    optparse$make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the
               script [default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    directory <- "/homes/hannah/data/ukbb/ukb-hrt"
    args$outdir <-"/homes/hannah/data/ukbb/ukb-hrt/grs"
    args$cohort <- "ukb"
    args$phenofile <- file.path(directory, "phenotypes",
                                "190402_fractal_dimension_26k",
                                "FD_summary_EUnorel.csv")
    args$covfile <- file.path(directory, "phenotypes",
                              "190402_fractal_dimension_26k",
                              "FD_covariates_EUnorel.csv")
    args$samplesfile <- file.path(directory, "rawdata",
                              "ukb18545_cal_chr1_v2_s488282.fam")
    args$verbose <- TRUE
}

## Read files ####
if (args$verbose) message("Read files with phenotypes and covariates")
pheno <- data.table::fread(args$phenofile, data.table=FALSE,
                               stringsAsFactors=FALSE)
covs <- data.table::fread(args$covfile, data.table=FALSE, stringsAsFactors=FALSE)

samples <- data.table::fread(args$samplesfile, data.table=FALSE, skip=2,
                             stringsAsFactors=FALSE)[,1:2]
colnames(samples) <- c("FID", "IID")

###############
## analysis ###
###############

if (args$verbose) message("Format results for grs with PRSice")
## ukb discovery cohort
colnames(pheno)[1] <- "IID"

if (args$cohort == "ukb") {
    colnames(covs)[1:6] <- c("IID", "sex", "age", "weight", "bmi", "height")
}
if (args$cohort == "dh") {
    covs <- covs[,c(ncol(covs), 1:3, 5, 4, 6:(ncol(covs)-1))]
}

all_pheno <- merge(samples, pheno, by="IID", all.x=TRUE, sort=FALSE)
all_pheno <- all_pheno[match(samples$IID, all_pheno[,1]),-2]

all_covs <- merge(samples, covs, by="IID", all.x=TRUE, sort=FALSE)
all_covs <- all_covs[match(samples$IID, all_covs[,1]),-2]

## write phenotype and summary stats file: one per pheno/region
prsice_single <- sapply(1:3, function(x) {
    nn <- colnames(pheno)[x+1]
    message("Write PRSice output for ", nn)
    write.table(pheno[,c(1,x+1)],
                file=paste(args$outdir, "/prsice_pheno_replication_", nn, "_",
                           args$cohort, ".txt", sep=""),
                sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)
})

## write covariates for replication
write.table(covs, file=paste(args$outdir, "/prsice_covariates_replication_",
                             args$cohort, ".txt", sep=""),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

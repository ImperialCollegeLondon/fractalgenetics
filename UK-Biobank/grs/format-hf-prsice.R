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
    optparse$make_option(c("-n", "--name"), action="store", dest="name",
               type="character", help="Name of analysis; has to be the same as
               in naming bgenie files [default: %default].", default=NULL),
    optparse$make_option(c("-hf", "--hf"), action="store", dest="hf",
               type="character", help="Path to hf id file; [default:%default].",
               default=NULL),
    optparse$make_option(c("-cad", "--cad"), action="store", dest="cad",
               type="character", help="Path to cad id file; [default:%default].",
               default=NULL),
    optparse$make_option(c("-icm", "--icm"), action="store", dest="icm",
               type="character", help="Path to icm id file; [default:%default].",
               default=NULL),
    optparse$make_option(c("-nicm", "--nicm"), action="store", dest="nicm",
               type="character", help="Path to nicm id file; [default:%default].",
               default=NULL),
    optparse$make_option(c("-sznicm", "--sznicm"), action="store", dest="sznicm",
               type="character", help="Path to sz_nicm id file;
               [default:%default].",
               default=NULL),
    optparse$make_option(c("-anicm", "--aragamnicm"), action="store",
                         dest="aragamnicm",
               type="character", help="Path to aragam nicm id file;
               [default:%default].",
               default=NULL),
    optparse$make_option(c("-id", "--idoverview"), action="store",
               dest="qcid", type="character", help="Path to .rds heart failure
               overview file; [default:%default].", default=NULL),
    optparse$make_option(c("-cov", "--covariates"), action="store",
               dest="covfile", type="character", help="Path to heart failure
               covariates file; [default:%default].", default=NULL),
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
    directory <-"/homes/hannah/data/ukbb/ukb-hrt"
    args <- list()
    args$outdir <- file.path(directory, "grs")
    args$samplesfile <- file.path(directory, "rawdata",
                              "ukb40616_imp_chr1_v3_s487317.sample")
    args$covfile <- file.path(directory, "heart_failure_phenotypes",
                              "heart_failure_covariates.csv")
    args$hf <- file.path(directory, "heart_failure_phenotypes", "hf_eid.csv")
    args$nicm <- file.path(directory, "heart_failure_phenotypes", "nicm_eid.csv")
    args$icm <- file.path(directory, "heart_failure_phenotypes", "icm_eid.csv")
    args$sznicm <- file.path(directory, "heart_failure_phenotypes",
                              "sz_nicm_eid.csv")
    args$aragamnicm <- file.path(directory, "heart_failure_phenotypes",
                                  "aragam_nicm_eid.csv")
    args$cad <- file.path(directory, "heart_failure_phenotypes", "cad_eid.csv")
    args$qcid <- file.path(directory, "heart_failure_phenotypes",
                           "heart_failure_samples_id.rds")
    args$verbose <- TRUE
}

## Read files ####
## IDs of samples in ukb with hear failure phenotypes (application 40616)
if (verbose) message("Read all sample ids of heart failure phenotypes")
hf <- data.table::fread(args$hf, data.table=FALSE, stringsAsFactors=FALSE)
icm <- data.table::fread(args$icm, data.table=FALSE, stringsAsFactors=FALSE)
nicm <- data.table::fread(args$nicm, data.table=FALSE, stringsAsFactors=FALSE)
sz_nicm <- data.table::fread(args$sznicm, data.table=FALSE,
                             stringsAsFactors=FALSE)
aragam_nicm <- data.table::fread(args$aragamnicm, data.table=FALSE,
                                 stringsAsFactors=FALSE)
cad <- data.table::fread(args$cad, data.table=FALSE, stringsAsFactors=FALSE)

if (args$verbose) message("Read files with heart failure covariates")
covs_hf_norelated <- data.table::fread(args$covfile, sep=",", data.table=FALSE,
                                       stringsAsFactors=FALSE)
overview_id <- readRDS(args$qcid)

if (args$verbose) message("Read genotype sample ids")
## heart failure phenotypes (application 40616)
samples <- data.table::fread(args$samplesfile, data.table=FALSE, skip=2,
                             stringsAsFactors=FALSE,
                             col.names=c("ID_1", "ID_2", "missing", "sex"))

###############
## analysis ###
###############

if (verbose) message("Format results for grs with PRSice")
## construct and write heart failure phenotypes and covariate files ####
failures <- list(hf=hf, cad=cad, icm=icm, nicm=nicm, sz_nicm=sz_nicm,
                 aragam_nicm=aragam_nicm)

formated <- lapply(seq_along(failures), function(x) {
        tmp <- failures[[x]]
        tmp$dummy <- rep(1, nrow(tmp))
        all <- merge(samples, tmp, by=1, all.x=TRUE, sort=FALSE)
        all <- all[match(samples$ID_1, all[,1]),ncol(all)]
        all[is.na(all)] <- 0
        return(all)
})
failures_df <- data.frame(IID=samples$ID_1, do.call(cbind, formated))
colnames(failures_df)[-1] <- names(failures)
failures_prsice <- failures_df

# set IIDs with missing covariates/related/non-white to NA
failures_prsice[!failures_prsice$IID %in% covs_hf_norelated$eid, -1] <- NA
write.table(failures_prsice, file=file.path(args$outdir,
                                 "prsice_heart_failures_ukb.txt"),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

# match covs with geno samples file
covs_all <- merge(samples[,1], covs_hf_norelated, by=1, all.x=TRUE, sort=FALSE)
covs_all <- covs_all[match(samples$ID_1, covs_all[,1]),]
colnames(covs_all)[1] <- "IID"
write.table(covs_all, file=file.path(args$outdir,
                                 "prsice_covariates_heart_failures_ukb.txt"),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

## qc breakdown ####
phenotyped <- length(overview_id$all)
length_pheno <- sapply(failures, nrow)
length_pheno <- c(phenotyped, length_pheno)
names(length_pheno)[1] <-  "all"

qc <- data.frame(sapply(overview_id, function(ids) {
                        tmp <- failures_df[failures_df$IID  %in% ids,]
                        sapply(tmp, function (x) sum(x != 0))
}))
qc <- cbind(length_pheno, qc)
colnames(qc)[1:2] <- c("phenotyped", "genotyped")
write.table(qc,
            file=file.path(args$outdir, "heart_failures_ukb_overview.csv"),
            sep=",", quote=FALSE, col.names=TRUE, row.names=TRUE)

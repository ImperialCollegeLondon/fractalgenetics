###############################
### Libraries and functions ###
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('optparse')

###########
## data ###
###########


## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path to directory with bgenie association
               results [default: %default].", default=NULL),
    optparse$make_option(c("-n", "--name"), action="store", dest="name",
               type="character", help="Name of analysis; has to be the same as
               in naming bgenie files [default: %default].", default=NULL),
    optparse$make_option(c("-p", "--pheno"), action="store", dest="phenofile",
               type="character", help="Path to fd phenotype file;
               [default:%default].", default=NULL),
    optparse$make_option(c("-cov", "--covariates"), action="store",
               dest="covfile", type="character", help="Path to fd covariates file;
               [default:%default].", default=NULL),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out [default: %default]."),
    optparse$make_option(c("-i", "--interpolate"), action="store",
               dest="interpolate",
               type="integer", help="Number of slices to interpolate to
               [default: %default].", default=9),
    optparse$make_option(c("--debug"), action="store_true",
                        dest="debug", default=FALSE, type="logical",
                        help="If set, predefined arguments are used to test the
                        script [default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$directory <-"/homes/hannah/data/ukbb/ukb-hrt/gwas/180628_fractal_dimension"
    args$basedirectory <-"/homes/hannah/data/ukbb/ukb-hrt"
    args$samples18545 <- "~/data/ukbb/ukb-hrt/rawdata/ukb18545_imp_chr1_v3_s487378.sample"
    args$samples40616 <- "~/data/ukbb/ukb-hrt/rawdata/ukb40616_imp_chr1_v3_s487317.sample"
    args$samplesDH <- "~/data/digital-heart/genotype/imputation/combined/genotypes/gencall.combined.clean.related.chr1.sample"
    args$name <-"slices"
    args$valueMeff <- NULL
    args$discoveryfile <- "/homes/hannah/data/ukbb/ukb-hrt/phenotypes/180628_fractal_dimension/FD_summary_bgenie.txt"
    args$cov_discoveryfile <- "/homes/hannah/data/ukbb/ukb-hrt/phenotypes/180628_fractal_dimension/FD_covariates_bgenie.txt"
    args$replicationfile <- "/homes/hannah/data/ukbb/ukb-hrt/phenotypes/190402_fractal_dimension_26k/FD_summary_bgenie.txt"
    args$cov_replicationfile <- "/homes/hannah/data/ukbb/ukb-hrt/phenotypes/190402_fractal_dimension_26k/FD_covariates_bgenie.txt"
    args$dhfile <- "/homes/hannah/data/digital-heart/phenotype/FD/FD_summary_bgenie.txt"
    args$cov_dhfile <- "/homes/hannah/data/digital-heart/phenotype/FD/FD_covariates_bgenie.txt"
    args$verbose <- TRUE
    args$interpolate <- 9
}

directory <- args$directory
basedirectory <- args$basedirectory
name <- args$name

discoveryfile <- args$discoveryfile
cov_discoveryfile <- args$cov_discoveryfile

replicationfile <- args$replicationfile
cov_replicationfile <- args$cov_replicationfile

dhfile <- args$dhfile
cov_dhfile <- args$cov_dhfile

samplesfile18545 <- args$samples18545
samplesfile40616 <- args$samples40616
samplesfileDH <- args$samplesDH
verbose <- args$verbose

## Read files ####
if (verbose) message("Read files with phenotypes and covariates")
## discovery data
discovery <- data.table::fread(discoveryfile, data.table=FALSE,
                               stringsAsFactors=FALSE)
covs_discovery <- data.table::fread(cov_discoveryfile, data.table=FALSE,
                               stringsAsFactors=FALSE)

## replication UKB
replication <- data.table::fread(replicationfile, data.table=FALSE,
                               stringsAsFactors=FALSE)
covs_replication <- data.table::fread(cov_replicationfile, data.table=FALSE,
                               stringsAsFactors=FALSE)

## replication digital heart
dh <- data.table::fread(dhfile, data.table=FALSE,
                               stringsAsFactors=FALSE)
covs_dh <- data.table::fread(cov_dhfile, data.table=FALSE,
                               stringsAsFactors=FALSE)

if (verbose) message("Read genotype sample ids")
# ukbb genotype samples via ukbgene imp -c1 -m
## ukb discovery and replication FD association (application 18545)
samples18545 <- data.table::fread(samplesfile18545, data.table=FALSE, skip=2,
                             stringsAsFactors=FALSE,
                             col.names=c("ID_1", "ID_2", "missing", "sex"))

## digital-heart replication FD association
samplesDH <- data.table::fread(samplesfileDH, data.table=FALSE, skip=2,
                             stringsAsFactors=FALSE)[,1:4]
colnames(samplesDH) <- c("ID_1", "ID_2", "missing", "sex")

## heart failure phenotypes (application 40616)
samples40616 <- data.table::fread(samplesfile40616, data.table=FALSE, skip=2,
                             stringsAsFactors=FALSE,
                             col.names=c("ID_1", "ID_2", "missing", "sex"))

## summary statistics of FD GWAS in ukb discovery cohort (application 18545)
if (verbose) message("Read files with association results")
genomewide <- data.table::fread(file=paste(directory, "/bgenie_", name,
                                           "_lm_st_genomewide.csv", sep=""),
                                sep=" ", data.table=FALSE,
                                stringsAsFactors=FALSE)

## all directly genotyped SNPs
if (verbose) message("Read file with rsids of directly genotyped snps")
genotyped <- data.table::fread(file=file.path(basedirectory, "genotypes",
                                              "ukb_cal_genome.bim"),
                                sep="\t", data.table=FALSE,
                                stringsAsFactors=FALSE)

## IDs of samples in ukb with hear failure phenotypes (application 40616)
if (verbose) message("Read all sample ids of heart failure phenotypes")
hf <- data.table::fread(file.path(basedirectory, "heart_phenotypes",
                                  "hf_eid.csv"),
                        data.table=FALSE, stringsAsFactors=FALSE)
nicm <- data.table::fread(file.path(basedirectory, "heart_phenotypes",
                                  "nicm_eid.csv"),
                        data.table=FALSE, stringsAsFactors=FALSE)
sz_nicm <- data.table::fread(file.path(basedirectory, "heart_phenotypes",
                                  "sz_nicm_eid.csv"),
                        data.table=FALSE, stringsAsFactors=FALSE)
aragam_nicm <- data.table::fread(file.path(basedirectory, "heart_phenotypes",
                                  "aragam_nicm_eid.csv"),
                        data.table=FALSE, stringsAsFactors=FALSE)
cad <- data.table::fread(file.path(basedirectory, "heart_phenotypes",
                                  "cad_eid.csv"),
                        data.table=FALSE, stringsAsFactors=FALSE)

covs_hf <- data.table::fread(file.path(basedirectory, "heart_phenotypes",
                                  "extra_fields.csv"),
                        sep=",", data.table=FALSE, stringsAsFactors=FALSE)

###############
## analysis ###
###############

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))
index_beta <- which(grepl("beta", colnames(genomewide)))

if (verbose) message("Format results for grs with PRSice")
## ukb discovery cohort
discovery$IID <- samples18545$ID_1
discovery <- discovery[,c(ncol(discovery),1:(ncol(discovery)-1))]

covs_discovery$IID <- samples18545$ID_1
covs_discovery <- covs_discovery[,c(ncol(covs_discovery),
                                    1:(ncol(covs_discovery)-1))]
colnames(covs_discovery)[1:6] <- c("IID", "sex", "age", "weight", "bmi", "height")

## ukb replication cohort
replication$IID <- samples18545$ID_1
replication <- replication[,c(ncol(replication),1:(ncol(replication)-1))]

covs_replication$IID <- samples18545$ID_1
covs_replication <- covs_replication[,c(ncol(covs_replication),
                                    1:(ncol(covs_replication)-1))]
colnames(covs_replication)[1:6] <- c("IID", "sex", "age", "weight", "bmi", "height")

## digital heart replication cohort
dh$IID <- samplesDH$ID_1
dh <- dh[,c(ncol(dh),1:(ncol(dh)-1))]

covs_dh$IID <- samplesDH$ID_1
covs_dh <- covs_dh[,c(ncol(covs_dh), 1:3, 5, 4, 6:(ncol(covs_dh)-1))]

## filter genomewide summary results for directly genotyped snps
genomewide <- dplyr:: filter(genomewide, rsid %in% genotyped[,2])
## write phenotype and summary stats file: one per pheno/region
prsice_single <- sapply(seq_along(index_beta), function(x) {
                   b <- index_beta[x]
                   p <- index_logp[x]
                   tmp <- genomewide[,1:5]
                   tmp$beta <- genomewide[,b]
                   tmp$p <- 10^-(genomewide[,p])
                   nn <- paste(name, "_",
                               gsub("_beta", "", colnames(genomewide)[b]),
                               sep="")
                   message("Write PRSice output for ", nn)
                   outdir <- file.path(basedirectory, "grs")
                   if(!dir.exists(outdir)) dir.create(outdir)
                   write.table(tmp, file=paste(outdir, "/prsice_association_",
                                               nn, ".txt", sep=""),
                               sep=" ", quote=FALSE, col.names=TRUE,
                               row.names=FALSE)
                   write.table(discovery[,c(1,x+1)],
                               file=paste(outdir, "/prsice_pheno_",
                                               nn, "_discovery_ukb.txt", sep=""),
                               sep=" ", quote=FALSE,
                               col.names=c("IID", "pheno"),
                               row.names=FALSE)
                   write.table(replication[,c(1,x+1)],
                               file=paste(outdir, "/prsice_pheno_",
                                               nn, "_replication_ukb.txt", sep=""),
                               sep=" ", quote=FALSE,
                               col.names=c("IID", "pheno"),
                               row.names=FALSE)
                   write.table(dh[,c(1,x+1)],
                               file=paste(outdir, "/prsice_pheno_",
                                               nn, "_replication_dh.txt", sep=""),
                               sep=" ", quote=FALSE,
                               col.names=c("IID", "pheno"),
                               row.names=FALSE)
                })

## write covariates for ukb discovery and replication and dh replication
write.table(covs_discovery, file=file.path(basedirectory, "grs",
                                 "prsice_covariates_summary_discovery_ukb.txt"),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(covs_replication, file=file.path(basedirectory, "grs",
                                 "prsice_covariates_summary_replication_ukb.txt"),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

write.table(covs_dh, file=file.path(basedirectory, "grs",
                                 "prsice_covariates_summary_replication_dh.txt"),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)




## construct and write heart failure phenofiles
failures <- list(hf=hf, cad=cad, nicm=nicm, sz_nicm=sz_nicm,
                 aragam_nicm=aragam_nicm)
length_pheno <- sapply(failures, nrow)

formated <- lapply(seq_along(failures), function(x) {
        tmp <- failures[[x]]
        tmp$dummy <- rep(1, nrow(tmp))
        all <- merge(samples40616, tmp, by=1, all.x=TRUE, sort=FALSE)
        all <- all[match(samples40616$ID_1, all[,1]),ncol(all)]
        all[is.na(all)] <- 0
        return(all)
})
length_geno <- sapply(formated, function (x) sum(x != 0))

failures_df <- data.frame(IID=samples40616$ID_1, do.call(cbind, formated))
colnames(failures_df)[-1] <- names(failures)

write.table(failures_df, file=file.path(basedirectory, "grs",
                                 "prsice_hffailures_ukb.txt"),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

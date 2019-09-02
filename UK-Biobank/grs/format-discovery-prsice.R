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
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path to ukb parent directory
               [default: %default].", default=NULL),
    optparse$make_option(c("-p", "--pheno"), action="store", dest="phenofile",
               type="character", help="Path to discovery fd phenotype file;
               [default:%default].", default=NULL),
    optparse$make_option(c("-cov", "--covariates"), action="store",
               dest="covfile", type="character", help="Path to discovery fd
               covariates file; [default:%default].", default=NULL),
    optparse$make_option(c("-s", "--samples"), action="store",
               dest="samplesfile", type="character", help="Path to bgen samples
               file; [default:%default].", default=NULL),
    optparse$make_option(c("-gwas", "--gwas"), action="store",
               dest="gwasfile", type="character", help="Path to discovery fd
               genomwide bgenie results file; [default:%default].",
               default=NULL),
    optparse$make_option(c("-g", "--genotypedids"), action="store",
               dest="genotypes", type="character", help="Path to .bim file with
               all directly genotyped markers; [default:%default].",
               default=NULL),
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
    args$directory <-"/homes/hannah/data/ukbb/ukb-hrt"
    args$genotypes <-file.path(args$directory, "genotypes", "ukb_cal_genome.bim")
    args$samplesfile <- file.path(args$directory, "rawdata",
                              "ukb18545_imp_chr1_v3_s487378.sample")
    args$phenofile <- file.path(args$directory, "phenotypes",
                                "180628_fractal_dimension",
                                "FD_summary_bgenie.txt")
    args$covfile <- file.path(args$directory, "phenotypes",
                              "180628_fractal_dimension",
                              "FD_covariates_bgenie.txt")
    args$gwasfile <- file.path(args$directory, "gwas",
                              "180628_fractal_dimension",
                              "bgenie_summary_lm_st_genomewide.csv")
    args$verbose <- TRUE
}

## Read files ####
if (args$verbose) message("Read files with phenotypes and covariates")
## discovery data
pheno <- data.table::fread(args$phenofile, data.table=FALSE,
                           stringsAsFactors=FALSE)
covs <- data.table::fread(args$covfile, data.table=FALSE,
                          stringsAsFactors=FALSE)


if (args$verbose) message("Read genotype sample ids")
# ukbb genotype samples via ukbgene imp -c1 -m
## ukb discovery and replication FD association (application 18545)
samples <- data.table::fread(args$samplesfile, data.table=FALSE, skip=2,
                             stringsAsFactors=FALSE,
                             col.names=c("ID_1", "ID_2", "missing", "sex"))

## summary statistics of FD GWAS in ukb discovery cohort (application 18545)
if (args$verbose) message("Read files with association results")
genomewide <- data.table::fread(file=paste(args$gwasfile,
                                sep=" ", data.table=FALSE,
                                stringsAsFactors=FALSE)

## all directly genotyped SNPs
if (args$verbose) message("Read file with rsids of directly genotyped snps")
genotyped <- data.table::fread(file=args$genotyped, sep="\t", data.table=FALSE,
                                stringsAsFactors=FALSE)

###############
## analysis ###
###############

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))
index_beta <- which(grepl("beta", colnames(genomewide)))

if (args$verbose) message("Format results for grs with PRSice")
## ukb discovery cohort
pheno$IID <- samples$ID_1
pheno <- pheno[,c(ncol(pheno),1:(ncol(pheno)-1))]

covs$IID <- samples$ID_1
covs <- covs[,c(ncol(covs), 1:(ncol(covs)-1))]
colnames(covs)[1:6] <- c("IID", "sex", "age", "weight", "bmi", "height")

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
                   outdir <- file.path(args$directory, "grs")
                   if(!dir.exists(outdir)) dir.create(outdir)
                   write.table(tmp, file=paste(outdir, "/prsice_association_",
                                               nn, ".txt", sep=""),
                               sep=" ", quote=FALSE, col.names=TRUE,
                               row.names=FALSE)
                   write.table(pheno[,c(1,x+1)],
                               file=paste(outdir, "/prsice_pheno_",
                                          nn, "_discovery_ukb.txt", sep=""),
                               sep=" ", quote=FALSE,
                               col.names=c("IID", "pheno"),
                               row.names=FALSE)
                })

## write covariates for ukb discovery
write.table(covs,
            file=file.path(args$directory, "grs",
                           "prsice_covariates_discovery_ukb.txt"),
            sep=" ", quote=FALSE, col.names=TRUE, row.names=FALSE)

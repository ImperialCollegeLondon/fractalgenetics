###############################
### Libraries and functions ###
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('optparse')
bgenie <- modules::import("bgenieResults")
plots <- modules::import("plots")
meta <- modules::import("metaanalysis")
ldfilter <- modules::import("LDfilter")

stPlots <- function(trait_index, bgenie_result, directory, is.negLog=FALSE,
                         Meff=1, ymin=1, ymax="max", name=NULL) {
    if (!is.null(name)) name <- paste("_", name, sep="")
    if (!is.negLog) {
        pvalues <- sapply(bgenie_result[,trait_index],
                          function(x) min(x * valueMeff, 1))
        which.sig <- pvalues < ymin
    } else {
        pvalues <- sapply(bgenie_result[,trait_index],
                          function(x) max(x - log10(valueMeff), 0))
        which.sig <- pvalues > -log10(ymin)
    }
    p_manhattan <- plots$manhattan(data.frame(bgenie_result[which.sig, 1:3],
                                              P=pvalues[which.sig]),
                                   title=colnames(bgenie_result)[trait_index],
                                   size.x.labels=12, size.y.labels=12,
                                   is.negLog=is.negLog,
                                   color=c("#fc8d59", "#b30000"), chr='chr',
                                   snp="rsid", bp="pos", max.y=ymax)
    ggplot2::ggsave(plot=p_manhattan, height=4, width=10,
                    file=paste(directory,"/bgenie", name, "_lm_",
                               colnames(bgenie_result)[trait_index],
                               "_manhattanplot.pdf", sep=""))
    p_qq <- plots$qqplot(bgenie_result[,trait_index], size.text=14,
                      size.title=14, is.negLog=is.negLog, raster=FALSE)
    ggplot2::ggsave(plot=p_qq, height=7, width=7,
                  file=paste(directory,"/bgenie", name, "_lm_",
                         colnames(bgenie_result)[trait_index],"_qqplot.png",
                    sep=""))
}

###############
## analysis ###
###############


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
               type="character", help="Path to fd phenotype file; if mEffective
               is not provided, --pheno can be used to estimate --mEffective. If
               neither are set, default mEffective = NrTraits_tested
               [default:%default].", default=NULL),
    optparse$make_option(c("-m", "--mEffective"), action="store",
               dest="valueMeff",
               type="double", help="Estimated number of effective test; based
               on Galwey (2009) Genetic Epidemiology [default:
               NrTraits_tested].", default=NULL),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out [default: %default]."),
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$directory <-"/homes/hannah/data/digital-heart/association/FD"
    args$name <-"slices"
    args$valueMeff <- NULL
    args$phenofile <- "/homes/hannah/data/digital-heart/association/FD_slices.csv"
    args$verbose <- TRUE
}

directory <- args$directory
name <- args$name
phenofile <- args$phenofile
verbose <- args$verbose
tags_prefix <- "European_ukb_imp_chr"
tags_suffix <- "_v3_maf0.001_500kb_r0.6.tags.list"
tags_dir <- "~/data/ukbb/ukb-hrt/tags"

genomewide <- lapply(1:22, bgenie$readBgenieOutput, directory=directory,
                     name=paste("bgenie_", name, "_lm_st_chr", sep=""),
                     maf=0.001, biallelicOnly=FALSE)
genomewide <- do.call(rbind, genomewide)

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))
index_beta <- which(grepl("beta", colnames(genomewide)))
index_t <- which(grepl("_t", colnames(genomewide)))

## write results ####
if (verbose) message("Write file with genome-wide association results")
write.table(genomewide,
            file=paste(directory, "/bgenie_", name, "_lm_st_genomewide.csv",
                       sep=""),
            sep=",",quote=FALSE, col.names=TRUE, row.names=FALSE)

## Meta-analysis single-trait summary statistics ####
if (verbose) message("Estimate pseudo-multitrait association results")
multitrait_res <- meta$pseudoMultitrait(genomewide[, index_t])
multitrait <- cbind(genomewide[,index_snp], multitrait_res$multitrait_p)
colnames(multitrait) <- c("CHR", "SNP","BP", "P")

if (verbose) message("Write pseudo-multitrait association results")
write.table(multitrait,
            file=paste(directory, "/bgenie_", name, "_lm_pseudomt_all.csv",
                       sep=""),
            sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
title="Pseudomulti-trait all FD slices"
if (verbose) message("Filter pseudo-multitrait association results for LD")
multitrait_ld <- ldfilter$filterSigLoci4LD(gwas=cbind(multitrait[,1:3],
                                                      genomewide[,4:7],
                                                      P=multitrait$P),
                                           tags_dir=tags_dir,
                                           tags_prefix=tags_prefix,
                                           tags_suffix=tags_suffix,
                                           matchBy="SNP", tags_ID="SNP",
                                           tagsSplit="|")
ldfilter$writeSig(multitrait_ld, threshold=5*10^(-8), directory=directory,
               name=paste("Pseudomultitrait", name, sep=""))

if (verbose) message("Plot pseudo-multitrait association results")
ymin <- 10^(-1)
ymax <- max(c(-log10(min(multitrait$P)), max(genomewide[,index_logp])))
sig <- multitrait[multitrait$P < ymin, ]
p_manhattan_mt <- plots$manhattan(d=sig,
                    title=title,
                    size.x.labels=12, size.y.labels=12, is.negLog=FALSE,
                    color=c("#fc8d59", "#b30000"), max.y=ymax)

ggplot2::ggsave(plot=p_manhattan_mt, height=4, width=10,
       file=paste(directory,"/bgenie_", name,
                  "_lm_pseudomt_manhattanplot.pdf", sep=""))

p_qq_mt <- plots$qqplot(multitrait$P, size.text=14, size.title=14)
ggplot2::ggsave(plot=p_qq_mt, height=7, width=7,
       file=paste(directory,"/bgenie_", name,
                  "_lm_pseudomt_qqplot.pdf", sep=""))


## per trait qq and manhattan plots ####
# effective number of tests to adjust single-trait assocation p-values by
if (!is.null(phenofile)) {
    pheno <- data.table::fread(phenofile, data.table=FALSE,
                               stringsAsFactors=FALSE)
    valueMeff <- meta$Teff(as.matrix(pheno[,-1]))
} else if (is.null(args$valueMeff)) {
    valueMeff <- length(index_logp)
} else {
    valueMeff <- args$valueMeff
}

# plot single-trait GWAS
if (verbose) message("Plot single-trait association results")
plots_perTrait <- sapply(index_logp, stPlots, bgenie_result=genomewide, ymin=0.1,
                         is.negLog=TRUE, directory=directory, Meff=valueMeff,
                         name=name,  ymax=ymax)

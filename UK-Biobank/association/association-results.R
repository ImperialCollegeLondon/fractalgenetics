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
        which.min <- pvalues < ymin
        which.sig <- pvalues < 5e-8
    } else {
        pvalues <- sapply(bgenie_result[,trait_index],
                          function(x) max(x - log10(valueMeff), 0))
        which.min <- pvalues > -log10(ymin)
        which.sig <- pvalues > -log10(5e-8)
    }
    highlight <- bgenie_result$rsid[which.sig]
    p_manhattan <- plots$manhattan(data.frame(bgenie_result[which.min, 1:3],
                                              P=pvalues[which.min]),
                                   title=colnames(bgenie_result)[trait_index],
                                   size.x.text=12, size.y.text=12,
                                   is.negLog=is.negLog,
                                   color=c('#878787', '#4d4d4d'), chr='chr',
                                   highlight=highlight, colorHighlight='#35978f',
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
    args$directory <-"/homes/hannah/data/ukbb/ukb-hrt/gwas"
    args$name <-"slices"
    args$valueMeff <- NULL
    args$phenofile <- "/homes/hannah/data/ukbb/ukb-hrt/phenotypes/FD_slices_EUnorel.csv"
    args$verbose <- TRUE
    args$interpolate <- 9
}

directory <- args$directory
name <- args$name
phenofile <- args$phenofile
verbose <- args$verbose
tags_prefix <- "European_ukb_imp_chr"
tags_suffix <- "_v3_maf0.001_250kb_r0.6.tags.list"
tags_dir <- "~/data/ukbb/ukb-hrt/tags/180628_fractal_dimension"

## Read results ####
if (verbose) message("Read files with genome-wide association results")
genomewide <- lapply(1:22, bgenie$readBgenieOutput, directory=directory,
                     name=paste("bgenie_", name, "_lm_st_chr", sep=""),
                     maf=0.001, biallelicOnly=FALSE)
genomewide <- do.call(rbind, genomewide)

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))
index_beta <- which(grepl("beta", colnames(genomewide)))
index_t <- which(grepl("_t", colnames(genomewide)))
ymax <- 'max'

## Write combined results ####
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
if (name == 'slices') {
    title="Pseudomulti-trait all FD slices"
} else {
    title="Pseudomulti-trait mean basal, mid and apical FD"
}

if (verbose) message("Plot pseudo-multitrait association results")
ymin <- 10^(-1)
ymax <- max(c(-log10(min(multitrait$P)), max(genomewide[,index_logp])))
sig <- multitrait[multitrait$P < ymin, ]
p_manhattan_mt <- plots$manhattan(d=sig,
                    title=title,
                    size.x.text=13, size.y.text=13,
                    size.x.title=15, size.y.title=15, is.negLog=FALSE,
                    color=c('#878787', '#4d4d4d'),
                    colorHighlight="#fc8d59",
                    highlight=multitrait$SNP[multitrait$P < 5e-8],
                    max.y=ymax)

ggplot2::ggsave(plot=p_manhattan_mt, height=4, width=11,
       file=paste(directory,"/bgenie_", name,
                  "_lm_pseudomt_manhattanplot.pdf", sep=""))

p_qq_mt <- plots$qqplot(multitrait$P, size.text=14, size.title=14)
ggplot2::ggsave(plot=p_qq_mt, height=7, width=7,
       file=paste(directory,"/bgenie_", name,
                  "_lm_pseudomt_qqplot.pdf", sep=""))

if (verbose) message("Filter pseudo-multitrait association results for LD")
multitrait_ld <-
    ldfilter$filterSigLoci4LD(gwas=cbind(multitrait[,1:3], P=multitrait$P,
                                         genomewide[,c(4:7, index_beta)]),
                              tags_dir=tags_dir,
                              tags_prefix=tags_prefix,
                              tags_suffix=tags_suffix,
                              matchBy="SNP", tags_ID="SNP",
                              tagsSplit="|")
ldfilter$writeSig(multitrait_ld, threshold=5*10^(-8), directory=directory,
               name="Pseudomultitrait_Slices")
sig <- genomewide[genomewide$rsid %in% multitrait_ld$sig_no_ld,]
write.table(sig, paste(directory, "/", name, "_sig5e08_ldFiltered.txt",
                              sep=""),
            col.names=TRUE, row.names=FALSE, quote=FALSE)

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

sig <- apply(as.matrix(genomewide[, index_logp]), 1,
             function(x, index) any(x > -log10(5e-8/valueMeff)))
sigAssociations <- genomewide[sig,]
write.table(sigAssociations, paste(directory, "/bgenie_", name,
                                   "_lm_st_significant.csv",  sep=""),
            quote=FALSE, col.names=TRUE, row.names=FALSE)


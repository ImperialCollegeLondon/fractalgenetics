###############################
### Libraries and functions ###
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

clustermq <- modules::import_package('clustermq')
plots <- modules::import("plots")
utils <- modules::import("utils")



stPlots <- function(trait_index, bgenie_result, directory, is.negLog=FALSE,
                         Meff=1) {
    if (!is.negLog) {
        pvalues <- sapply(bgenie_result[,trait_index],
                          function(x) min(x * valueMeff, 1))
        which.sig <- pvalue < Ymin
    } else {
        pvalues <- sapply(bgenie_result[,trait_index],
                          function(x) max(x - log10(valueMeff), 0))
        which.sig <- pvalue > -log10(Ymin)
    }
    p_manhattan <- plots$manhattan(data.frame(bgenie_result[which.sig, 1:3],
                                              P=pvalues[which.sig]),
                                   title=colnames(bgenie_result)[trait_index],
                                   size.x.labels=12, size.y.labels=12,
                                   is.negLog=is.negLog,
                                   color=c("#fc8d59", "#b30000"), chr='chr',
                                   snp="rsid", bp="pos")
    ggplot2::ggsave(plot=p_manhattan, height=4, width=10,
                    file=paste(directory,"/bgenie_lm_",
                               colnames(bgenie_result)[trait_index],
                               "_manhattanplot.pdf", sep=""))
    p_qq <- qp$qqplot(bgenie_result[which.sig,trait_index], size.text=14,
                      size.title=14)
    ggplot2::ggsave(plot=p_qq, height=7, width=7,
                    file=paste(directory,"/bgenie_lm_",
                         colnames(bgenie_result)[trait_index],"_qqplot.pdf",
                    sep=""))
}

###############
## analysis ###
###############

## command line arguments ####
option_list <- list(
    make_option(c("-d", "--directory"), action="store", dest="directory",
               type="character", help="Path to directory with bgenie association
               results [default: %default].", default=NULL),
    make_option(c("-n", "--name"), action="store", dest="name",
               type="character", help="Name of analysis; has to be the same as
               in naming bgenie files [default: %default].", default=NULL),
    make_option(c("-p", "--pheno"), action="store", dest="phenofile",
               type="character", help="Path to fd phenotype file; if mEffective
               is not provided, --pheno can be used to estimate --mEffective. If
               neither are set, default mEffective = NrTraits_tested
               [default:%default].", default=NULL),
    make_option(c("-m", "--mEffective"), action="store", dest="valueMeff",
               type="double", help="Estimated number of effective test; based
               on Galwey (2009) Genetic Epidemiology [default:
               NrTraits_tested].", default=NULL))

args <- optparse$parse_args(OptionParser(option_list=option_list))

if (FALSE) {
    args <- list()
    args$directory <-"~/data/ukbb/ukb-hrt/gwas"
    args$name <-"summary_EDV"
    args$valueMeff <- NULL
    args$phenofile <- "~/data/ukbb/ukb-hrt/phenotypes/FD_phenotypes_EUnorel.csv"
}

directory <- args$directory
name <- args$name
phenofile <- args$phenofile

## read bgenie output files per chromosome ####
genomewide <- lapply(1:22, function(chr) {
    chrFile <- paste(directory, "/bgenie_", name, "_lm_st_chr", chr,".gz",
                     sep="")
    if (!file.exists(chrFile)) {
        message("Bgenie association output for chr", chr, " cannot be found, ",
        "skip to next chromosome")
        return(NULL)
    }
    readString <- paste("zcat", chrFile)
    tmp <- data.table::fread(readString, sep=" ", stringsAsFactors=FALSE,
                             data.table=FALSE, header=TRUE)
    tmp <- tmp[!tmp$af %in% c(0,1),]
    return(tmp)
})
genomewide <- do.call(rbind, genomewide)

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))
index_beta <- which(grepl("beta", colnames(genomewide)))
index_t <- which(grepl("_t", colnames(genomewide)))

## write results ####
write.table(genomewide,
            file=paste(directory, "/bgenie_", name, "lm_st_genomewide.csv",
                       sep=""),
sep=",",quote=FALSE, col.names=TRUE, row.names=FALSE)

## Meta-analysis single-trait summary statistics ####
multitrait <- utils$pseudoMultitrait(genomewide[, index_t])
multitrait <- cbind(genomewide[,index_snp], multitrait)
colnames(multitrait) <- c("CHR", "SNP","BP", "P")
write.table(multitrait,
            file=paste(directory, "/bgenie_", name, "lm_pseudomt_all.csv",
                       sep=""),
            sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

ymin = 10^(-1)
sig <- multitrait[multitrait$P < ymin, ]
p_manhattan_mt <- plots$manhattan(d=sig,
                    title="Pseudomulti-trait all",
                    size.x.labels=12, size.y.labels=12, is.negLog=FALSE,
                    color=c("#fc8d59", "#b30000"))

ggplot2::ggsave(plot=p_manhattan_mt, height=4, width=10,
       file=paste(directory,"/bgenie_lm_", name,
                  "_lm_pseudomt_manhattanplot.pdf", sep=""))

p_qq_mt <- plots$qqplot(multitrait$P, size.text=14, size.title=14)
ggplot2::ggsave(plot=p_qq_mt, height=7, width=7,
       file=paste(directory,"/bgenie_", name,
                  "_lm_pseudomt_qqplot.png", sep=""))

if (grepl("summary", name)) {
    multitrait_maxApicalBasal <-
        utils$pseudoMultitrait(genomewide[, index_t][,4:5])
    multitrait_maxApicalBasal <- cbind(genomewide[,index_snp],
        multitrait_maxApicalBasal)
    colnames(multitrait_maxApicalBasal) <- c("CHR", "SNP","BP", "P")
    write.table(multitrait_maxApicalBasal,
                file=paste(directory, "/bgenie_", name,
                           "_lm_pseudomt_maxApicalBasal.csv", sep=""),
                sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

     sig <- multitrait_maxApicalBasal[multitrait_maxApicalBasal$P < ymin, ]
     p_manhattan_mt <- plots$manhattan(d=sig,
                    title="Pseudomulti-trait max Basal and Apical FD",
                    size.x.labels=12, size.y.labels=12, is.negLog=FALSE,
                    color=c("#fc8d59", "#b30000"))
     ggplot2::ggsave(plot=p_manhattan_mt, height=4, width=10,
            file=paste(directory,"/bgenie_", name,
                 "_lm_pseudomt_manhattanplot_maxBasalApical.pdf", sep=""))

     p_qq_mt <- plots$qqplot(sig$P, size.text=14, size.title=14)
     ggplot2::ggsave(plot=p_qq_mt, height=7, width=7,
       file=paste(directory,"/bgenie_", name, "_lm_pseudomt_qqplot.pdf",
                  sep=""))
}


## per trait qq and manhattan plots ####
# effective number of tests to adjust single-trait assocation p-values by
if (!is.null(phenofile)) {
    pheno <- data.table::fread(phenofile, data.table=FALSE,
                               stringsAsFactors=FALSE)
    valueMeff <- utils$Teff(as.matrix(pheno[,-1]))
} else if (is.null(args$valueMeff)) {
    valueMeff <- length(index_logp)
} else {
    valueMeff <- args$valueMeff
}

# plot single-trait GWAS
plots_perTrait <- sapply(index_logp, stPlots, bgenie_result=genomewide,
                         is.negLog=TRUE, directory=directory, Meff=valueMeff)




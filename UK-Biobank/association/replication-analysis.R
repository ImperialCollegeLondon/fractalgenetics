###############################
### Libraries and functions ###
###############################

options(import.path=c("/homes/hannah/projects/GWAS",
                      "/homes/hannah/GWAS/analysis/fd",
                      "/homes/hannah/projects"))
options(bitmapType = 'cairo', device = 'pdf')

modules::import_package('ggplot2', attach=TRUE)
smooth <- modules::import('utils/smoothAddR2')
optparse <- modules::import_package('optparse')
ldfilter <- modules::import("LDfilter")

###############
## analysis ###
###############


## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path to directory with bgenie association
               results [default: %default].", default=NULL),
    optparse$make_option(c("-nl", "--numberloci"), action="store", dest="nloci",
               type="numeric", help="Number of significant, independent loci
               in discovery [default:%default].", default=1),
    optparse$make_option(c("-sl", "--sigloci"), action="store", dest="sloci",
               type="character", help="Path to file with replication association
               results filtered for significance in discovery
               [default:%default].", default=NULL),
    optparse$make_option(c("--debug"), action="store_true",
                        dest="debug", default=FALSE, type="logical",
                        help="If set, predefined arguments are used to test the
                        script [default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$directory <-"/homes/hannah/data/ukbb/ukb-hrt/gwas/190402_fractal_dimension_26k/"
    args$sloci <-"/homes/hannah/data/ukbb/ukb-hrt/gwas/190402_fractal_dimension_26k/"
    args$nloci <- 16
    args$discoverydir <- "/homes/hannah/data/ukbb/ukb-hrt/gwas/180628_fractal_dimension"
}

tags_prefix <- "European_ukb_imp_chr"
tags_suffix <- "_v3_maf0.001_250kb_r0.6.tags.list"
tags_dir <- "~/data/ukbb/ukb-hrt/tags/180628_fractal_dimension"

replication_of_sig_loci <- data.table::fread(args$sloci,
                                 data.table=FALSE, stringsAsFactors=FALSE)
replication_sig_ld <- ldfilter$filterSigLoci4LD(gwas=replication_of_sig_loci,
                                                tags_dir=tags_dir,
                                                tags_prefix=tags_prefix,
                                                tags_suffix=tags_suffix,
                                                threshold = 1,
                                                matchBy="SNP", tags_ID="SNP",
                                                tagsSplit="|")

write.table(replication_sig_ld$sig,
            file.path(args$directory,
                      "Replication_association_significant_in_discovery.txt"),
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(replication_sig_ld$sig_wo_ld,
            file.path(args$directory,
                      "Replication_association_significant_in_discovery_ldFiltered.txt"),
            sep="\t", col.names=TRUE, row.names=FALSE)

slices_rep <- read.table(paste(args$dir,
                               "/bgenie_slices_lm_st_genomewide_replication.csv",
                               sep=""),
                         sep=",", stringsAsFactors=FALSE, header=TRUE)

## ld-filtered, significant genome-wide association results ukb ####
slices_dsv <- read.table(paste(args$discoverydir,
                               "/Pseudomultitrait_slices_sig5e08_ldFiltered.txt",
                               sep=""),
                         sep=",", stringsAsFactors=FALSE, header=TRUE)
# LD filter misses these two SNPs, manually remove
slices_dsv <- slices_dsv[!slices_dsv$rsid %in% c("rs12214483", "rs117953218"),]

## format ukb betas and p-values per slice ####
slices_dsv_beta <- cbind(rsid=slices_dsv$rsid,
                          slices_dsv[,grepl('beta', colnames(slices_dsv))])
colnames(slices_dsv_beta) <- gsub("_beta", "", colnames(slices_dsv_beta))
beta <- reshape2::melt(slices_dsv_beta, id.var='rsid', value.name='beta',
                       variable.name='slice')

slices_dsv_logp <- cbind(rsid=slices_dsv$rsid,
                          slices_dsv[, grepl('log10p', colnames(slices_dsv))])
colnames(slices_dsv_logp) <- gsub("\\.log10p", "", colnames(slices_dsv_logp))
logp <- reshape2::melt(slices_dsv_logp, id.var='rsid', value.name='logp',
                       variable.name='slice')

dsv <- cbind(beta, logp=logp$logp)
dsv$rsid <- as.character(dsv$rsid)

## discovery significant slice pvalues and betas ####
dsv <- dsv[dsv$rsid %in% slices_rep$rsid, ]
dsv_sig_slices <- dplyr::filter(dsv, logp  > -log10(5e-8))

slices2sample <- as.data.frame(table(dsv_sig_slices$slice),
                               stringsAsFactors=FALSE)
colnames(slices2sample) <- c("slice", "freq")
slices2sample <- slices2sample[slices2sample$freq != 0,]

observedBetas <- dsv_sig_slices$beta

## betas and pvalues of significant snps and slice in discovery ####
rep_beta <- slices_rep[,c(2, which(grepl('beta', colnames(slices_rep))))]
colnames(rep_beta) <- gsub("_beta", "", colnames(rep_beta))
rep_logp <- slices_rep[,c(2, which(grepl('log10p', colnames(slices_rep))))]
colnames(rep_logp) <- gsub("\\.log10p", "", colnames(rep_logp))

rep_sig_slices <- data.frame(t(apply(as.matrix(dsv_sig_slices), 1, function(x) {
    pos <- which(rep_beta$rsid %in%  x[1])
    beta_tmp <- rep_beta[pos, colnames(rep_beta) %in% x[2]]
    logp_tmp <- rep_logp[pos, colnames(rep_logp) %in% x[2]]
    return(rbind(x[2], beta_tmp, logp_tmp))
})), stringsAsFactors=FALSE)
colnames(rep_sig_slices) <- c("slice", "beta", "logp")

rep_sig_slices$rsid <- dsv_sig_slices$rsid
rep_sig_slices <- dplyr::select(rep_sig_slices, rsid, slice, beta, logp)
rep_sig_slices$beta <- as.numeric(rep_sig_slices$beta)
rep_sig_slices$logp <- as.numeric(rep_sig_slices$logp)

## concordance of effect sizes ####
nrObservations <- length(sign(rep_sig_slices$beta * observedBetas))
observedConcordance <- sum(sign(rep_sig_slices$beta * observedBetas))

## Empirical concordance: matched for number of sig snps and slices ####
nrsig <- nrow(dsv_sig_slices)
nrtotal <- nrow(slices_rep)
draws <- 100000
seed <- 101

set.seed(seed)
testConcordance <- sapply(1:draws, function(dummy) {
    randomSnps <- rep_beta[sample(nrtotal, nrsig),
                          colnames(rep_beta) %in% slices2sample$slice]
    randomBetas <- unlist(sapply(1:nrow(slices2sample), function(x) {
        pos <- colnames(randomSnps) == slices2sample$slice[x]
        sample(randomSnps[,pos], slices2sample$freq[x])
    }))
    return(sum(sign(randomBetas * observedBetas)))
})

empiricalConcordance <-
    length(which(testConcordance >= observedConcordance))/draws
if (empiricalConcordance == 0) {
    empiricalConcordance <- 1/draws
}
concordance <- data.frame(observedConcordance=
                            (nrObservations - observedConcordance)/2,
                          empiricalConcordance=empiricalConcordance,
                          nrObservations=nrObservations)
write.table(concordance,
            paste(args$dir, "/Replication_concordance_summary.txt", sep=""),
            quote=FALSE, col.names=TRUE, row.names=FALSE)

## plot beta concordance discovery and replication ####
slices <- cbind(dsv_sig_slices, rep_sig_slices[,3:4])
colnames(slices) <- c('rsid', 'slices', 'discovery_beta', 'discovery_logp',
                      'replication_beta', 'replication_logp')

sig_adjust <- round(0.05/args$nloci,3)
slices$sig <- factor(as.numeric(slices$replication_logp > -log10(sig_adjust)),
                        labels=c(expression(p >= sig_adjust),
                                 expression(p < sig_adjust)))
slices$concordance <- factor(-sign(slices$discovery_beta * slices$replication_beta),
                             labels=c('yes', 'no'))

write.table(slices, paste(args$dir, "/Replication_concordance.txt",
                          sep=""),
            quote=FALSE, col.names=TRUE, row.names=FALSE)

max_y <- max(abs(slices$replication_beta))
max_x <- max(abs(slices$discovery_beta))
p <- ggplot(data=slices, aes(x=discovery_beta, y=replication_beta))
p <- p + geom_point(aes(color=concordance, shape=sig)) +
        smooth$stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE,
                                xpos=max_x - 1/5*max_x,
                                ypos=max_y + 1/10*max_y, vjust=0, color="black") +
        xlim(c(-max_x - 1/5*max_x, max_x + 1/5*max_x)) +
        ylim(c(-max_y - 1/5*max_y, max_y + 1/5*max_y)) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0) +
        xlab(expression(hat(beta)[discovery])) +
        ylab(expression(hat(beta)[replication])) +
        scale_color_manual(values=c('black', '#969696'), guide=FALSE) +
        scale_shape_manual(values=c(20, 17), name='Replication',
                           labels=c(bquote(p>=.(sig_adjust)),
                                    bquote(p<.(sig_adjust)))) +
        theme_bw()
ggsave(plot=p, paste(args$dir, "/Replication_concordance.pdf", sep=""),
       height=5, width=6.5)



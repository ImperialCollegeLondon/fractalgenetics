###############################
### Libraries and functions ###
###############################
options(import.path=c("/homes/hannah/projects/GWAS",
                      "/homes/hannah/projects"))
options(bitmapType = 'cairo', device = 'pdf')


modules::import_package('ggplot2', attach=TRUE)
smooth <- modules::import('utils/smoothAddR2')
optparse <- modules::import_package('optparse')
bgenie <- modules::import("bgenieResults")

###############
## analysis ###
###############

## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path to directory with digital-heart
               bgenie association results [default: %default].", default=NULL),
    optparse$make_option(c("-ukb", "--ukbdir"), action="store",
               dest="ukbdir",
               type="character", help="Path to directory with ukbb significant
               association results [default: %default].", default=NULL),
    optparse$make_option(c("-n", "--name"), action="store", dest="name",
               type="character", help="Name of analysis; has to be the same as
               in naming bgenie files [default: %default].", default=NULL),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out [default: %default]."),
    optparse$make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the script 
               [default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$directory <- "~/data/digital-heart/gwas/FD"
    args$ukbdir <- "~/data/ukbb/ukb-hrt/gwas/180628_fractal_dimension"
    args$name <- 'slices'
    args$verbose <- TRUE
}
directory <- args$directory
ukbdir <- args$ukbdir
name <- args$name
verbose <- args$verbose

## genome-wide association results digital-heart ####
slices_dh <- lapply(1:22, bgenie$readBgenieOutput, directory=directory,
                     name=paste("bgenie_", name, "_lm_st_chr", sep=""),
                     maf=0.001, biallelicOnly=FALSE)
slices_dh <- do.call(rbind, slices_dh)


if (verbose) message("Write file with genome-wide association results")
write.table(slices_dh,
            file=paste(directory, "/bgenie_", name, "_lm_st_genomewide.csv",
                       sep=""),
            sep=",",quote=FALSE, col.names=TRUE, row.names=FALSE)


## ld-filtered, significant genome-wide association results ukb ####
slices_ukb <- read.table(paste(ukbdir,
                               "/Pseudomultitrait_slices_sig5e08_ldFiltered.txt",
                               sep=""),
                         sep=",", stringsAsFactors=FALSE, header=TRUE)
# LD filter misses these two SNPs, manually remove
slices_ukb <- slices_ukb[!slices_ukb$rsid %in% c("rs12214483", "rs117953218"),]

## format ukb betas and p-values per slice ####
slices_ukb_beta <- cbind(rsid=slices_ukb$rsid,
                          slices_ukb[,grepl('beta', colnames(slices_ukb))])
colnames(slices_ukb_beta) <- gsub("_beta", "", colnames(slices_ukb_beta))
beta <- reshape2::melt(slices_ukb_beta, id.var='rsid', value.name='beta',
                       variable.name='slice')

slices_ukb_logp <- cbind(rsid=slices_ukb$rsid,
                          slices_ukb[, grepl('log10p', colnames(slices_ukb))])
colnames(slices_ukb_logp) <- gsub("\\.log10p", "", colnames(slices_ukb_logp))
logp <- reshape2::melt(slices_ukb_logp, id.var='rsid', value.name='logp',
                       variable.name='slice')

ukb <- cbind(beta, logp=logp$logp)
ukb$rsid <- as.character(ukb$rsid)

## ukb significant slice pvalues and betas ####
ukb <- ukb[ukb$rsid %in% slices_dh$rsid, ]
ukb_sig_slices <- dplyr::filter(ukb, logp  > -log10(5e-8))

slices2sample <- as.data.frame(table(ukb_sig_slices$slice),
                               stringsAsFactors=FALSE)
colnames(slices2sample) <- c("slice", "freq")
slices2sample <- slices2sample[slices2sample$freq != 0,]

observedBetas <- ukb_sig_slices$beta

## betas and pvalues of significant snps and slice in ukb ####
dh_beta <- slices_dh[,c(2, which(grepl('beta', colnames(slices_dh))))]
colnames(dh_beta) <- gsub("_beta", "", colnames(dh_beta))
dh_logp <- slices_dh[,c(2, which(grepl('log10p', colnames(slices_dh))))]
colnames(dh_logp) <- gsub("-log10p", "", colnames(dh_logp))

dh_sig_slices <- data.frame(t(apply(ukb_sig_slices, 1, function(x) {
    pos <- which(dh_beta$rsid %in%  x[1])
    beta_tmp <- dh_beta[pos, colnames(dh_beta) %in% x[2]]
    logp_tmp <- dh_logp[pos, colnames(dh_logp) %in% x[2]]
    return(rbind(x[2], beta_tmp, logp_tmp))
})), stringsAsFactors=FALSE)
colnames(dh_sig_slices) <- c("slice", "beta", "logp")
dh_sig_slices$rsid <- ukb_sig_slices$rsid
dh_sig_slices <- dplyr::select(dh_sig_slices, rsid, slice, beta, logp)
dh_sig_slices$beta <- as.numeric(dh_sig_slices$beta)
dh_sig_slices$logp <- as.numeric(dh_sig_slices$logp)

## concordance of effect sizes ####
nrObservations <- length(sign(dh_sig_slices$beta * observedBetas))
observedConcordance <- sum(sign(dh_sig_slices$beta * observedBetas))

## Empirical concordance: matched for number of sig snps and slices ####
nrsig <- nrow(ukb_sig_slices)
nrtotal <- nrow(slices_dh)
draws <- 100000
seed <- 101

set.seed(seed)
testConcordance <- sapply(1:draws, function(dummy) {
    randomSnps <- dh_beta[sample(nrtotal, nrsig),
                          colnames(dh_beta) %in% slices2sample$slice]
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
            paste(directory, "/", name, "_concordance_summary.txt", sep=""),
            quote=FALSE, col.names=TRUE, row.names=FALSE)

## plot beta concordance ukb and digital heart ####
slices <- cbind(ukb_sig_slices, dh_sig_slices[,3:4])
colnames(slices) <- c('rsid', 'slices', 'ukbb_beta', 'ukbb_logp', 'dh_beta',
                      'dh_logp')

slices$dh_p <- 10^-(slices$dh_logp)
slices$ukbb_p <- 10^-(slices$ukbb_logp)

sig_adjust <- round(0.05/args$nloci,3)
slices$sig <- factor(as.numeric(slices$dh_logp > -log10(sig_adjust)),
                        labels=c(expression(p >= sig_adjust),
                                 expression(p < sig_adjust)))
slices$concordance <- factor(-sign(slices$ukbb_beta * slices$dh_beta),
                             labels=c('yes', 'no'))

write.table(slices, paste(directory, "/", name, "_concordance.txt",
                          sep=""),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

max_y <- max(abs(slices$dh_beta))
max_x <- max(abs(slices$ukbb_beta))
p <- ggplot(data=slices, aes(x=ukbb_beta, y=dh_beta))
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

ggsave(plot=p, paste(directory, "/", name, "_concordance.pdf", sep=""),
       height=5, width=6.5)








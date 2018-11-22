###############################
### Libraries and functions ###
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

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
               about analyses are printed to standard out [default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (FALSE) {
    args <- list()
    args$directory <- "~/data/digital-hea"
    args$ukbdir <- "~/data/ukbb/ukb-hrt/gwas"
    args$name <- 'volumes'
    args$verbose <- TRUE
}
directory <- args$directory
ukbdir <- args$ukbdir
name <- args$name
verbose <- args$verbose

## genome-wide association results volumes ####
lvv <- data.table::fread(paste(ukbdir, "/bgenie_", name, "_lm_st_genomewide.csv",
                               sep=""), data.table=FALSE, stringsAsFactors=FALSE)

## ld-filtered, significant genome-wide association results ukb ####
slices_sig <- read.table(paste(ukbdir, "/Slices_sig5e08_ldFiltered.txt", sep=""),
                         sep=",", stringsAsFactors=FALSE, header=TRUE)

## format slice betas and p-values per slice ####
slices_beta <- cbind(rsid=slices_sig$rsid,
                          slices_sig[,grepl('beta', colnames(slices_sig))])
colnames(slices_beta) <- gsub("_beta", "", colnames(slices_beta))
beta <- reshape2::melt(slices_beta, id.var='rsid', value.name='beta',
                       variable.name='slice')

slices_se <- cbind(rsid=slices_sig$rsid,
                          slices_sig[,grepl('se', colnames(slices_sig))])
colnames(slices_se) <- gsub("_se", "", colnames(slices_sig))
se <- reshape2::melt(slices_se, id.var='rsid', value.name='se',
                       variable.name='slice')

slices_logp <- cbind(rsid=slices_sig$rsid,
                          slices_sig[, grepl('log10p', colnames(slices_sig))])
colnames(slices_logp) <- gsub("\\.log10p", "", colnames(slices_logp))
logp <- reshape2::melt(slices_logp, id.var='rsid', value.name='logp',
                       variable.name='slice')

slices <- cbind(beta, logp=logp$logp, se=se$se)
slices$rsid <- as.character(slices$rsid)

## slices significant slice pvalues and betas ####
slices <- slices[slices$rsid %in% lvv$rsid, ]
sig_per_slice <- dplyr::filter(slices, logp  > -log10(5e-8))

## betas and pvalues of significant snps and slice in ukb ####
lvv_beta <- lvv[,c(2, which(grepl('beta', colnames(lvv))))]
colnames(lvv_beta) <- gsub("_beta", "", colnames(lvv_beta))
lvv_logp <- lvv[,c(2, which(grepl('log10p', colnames(lvv))))]
colnames(lvv_logp) <- gsub("-log10p", "", colnames(lvv_logp))
lvv_se <- lvv[,c(2, which(grepl('se', colnames(lvv))))]
colnames(lvv_se) <- gsub("_se", "", colnames(lvv_se))

lvv_sig <- do.call(rbind, apply(sig_per_slice, 1, function(x) {
    pos <- which(lvv_beta$rsid %in%  x[[1]])
    beta_tmp <- lvv_beta[pos, -1]
    names(beta_tmp) <- paste(names(beta_tmp), "_beta", sep="")
    logp_tmp <- lvv_logp[pos, -1]
    names(logp_tmp) <- paste(names(logp_tmp), "_logp", sep="")
    se_tmp <- lvv_se[pos, -1]
    names(se_tmp) <- paste(names(se_tmp), "_se", sep="")
    tmp <- data.frame(beta_tmp, logp_tmp, se_tmp, stringsAsFactors=FALSE)
    return(tmp)
}))

## plot betas slices and volumes ####
combined <- cbind(sig_per_slice, lvv_sig)
colnames(slices)[1:5] <- c('rsid', 'slices', 'slices_beta', 'slices_logp',
                           'slices_se')


p <- ggplot(data=combined)
p <- p + geom_point(aes(x=slices_beta, y=CO_beta)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    xlab(expression(hat(beta)[slices])) +
    ylab(expression(hat(beta)[CO])) +
    geom_errorbar(aes(ymin=CO_beta - CO_se, ymax=CO_beta + CO_se, x=slices_beta)) +
    geom_errorbarh(aes(xmin=slices_beta-slices_se, xmax=slices_beta + slices_se,
                       y=CO_beta)) +
    facet_wrap(~slices) +
    theme_bw()
ggsave(plot=p, paste(directory, "/", name, "_MR_CO.pdf", sep=""),
       height=5, width=6.5)

p <- ggplot(data=combined)
p <- p + geom_point(aes(x=slices_beta, y=SV_beta)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    xlab(expression(hat(beta)[slices])) +
    ylab(expression(hat(beta)[SV])) +
    geom_errorbar(aes(ymin=SV_beta - SV_se, ymax=SV_beta + SV_se, x=slices_beta)) +
    geom_errorbarh(aes(xmin=slices_beta-slices_se, xmax=slices_beta + slices_se,
                       y=SV_beta)) +
    facet_wrap(~slices) +
    theme_bw()
ggsave(plot=p, paste(directory, "/", name, "_MR_SV.pdf", sep=""),
       height=5, width=6.5)


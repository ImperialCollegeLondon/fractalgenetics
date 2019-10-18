###############################
### Libraries and functions ###
###############################
options(bitmapType = 'cairo', device = 'pdf')

modules::import_package('dplyr', attach=TRUE)
modules::import_package('ggplot2', attach=TRUE)
optparse <- modules::import_package('optparse')


############
## data  ###
############

## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path ukbb root data rootdirectory
                [default: %default].", default=NULL),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out ",
               "[default: %default]."),
    optparse$make_option(c("--debug"), action="store_true",
                dest="debug", default=FALSE, type="logical",
                help="If set, predefined arguments are used to test the script",
                "[default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$directory <- "~/data/ukbb/ukb-hrt"
    args$verbose <- TRUE
    args$Teff <- 6.6
}
directory <- args$directory
Teff <- args$Teff
verbose <- args$verbose

## ld-filtered, significant genome-wide association results ukb ####
slices_sig <- read.table(paste(directory,
                               "/gwas/Pseudomultitrait_slices_sig5e08_ldFiltered.txt",
                               sep=""),
                         sep=",", stringsAsFactors=FALSE, header=TRUE)
# LD filter misses these two SNPs, manually remove
slices_sig <- slices_sig[!slices_sig$rsid %in% c("rs12214483", "rs117953218"),]

slices_sig$SNPID <- paste(slices_sig$chr, ":", slices_sig$pos, "_",
                          slices_sig$a_0, "_", slices_sig$a_1, sep="")


###############
## analysis ###
###############

## format slice betas and p-values per slice ####
slices_beta <- cbind(rsid=slices_sig$rsid, chr=slices_sig$chr,
                     pos=slices_sig$pos,
                     a_0=slices_sig$a_0,
                     a_1=slices_sig$a_1, af=slices_sig$af, samplesize=18097,
                     slices_sig[,grepl('beta', colnames(slices_sig))])
colnames(slices_beta) <- gsub("_beta", "", colnames(slices_beta))
beta <- reshape2::melt(slices_beta, id.var=c('rsid', 'chr', 'pos',
                                             'a_0', 'a_1', 'af',
                                             'samplesize'),
                                             value.name='beta',
                       variable.name='slice')

slices_se <- cbind(rsid=slices_sig$rsid, a_0=slices_sig$a_0,
                   a_1=slices_sig$a_1,
                   slices_sig[,grepl('se', colnames(slices_sig))])
colnames(slices_se) <- gsub("_se", "", colnames(slices_se))
se <- reshape2::melt(slices_se, id.var=c('rsid', 'a_0', 'a_1'),
                     value.name='se',
                     variable.name='slice')

slices_logp <- cbind(rsid=slices_sig$rsid, a_0=slices_sig$a_0,
                     a_1=slices_sig$a_1,
                     slices_sig[, grepl('log10p', colnames(slices_sig))])
colnames(slices_logp) <- gsub("\\.log10p", "", colnames(slices_logp))
logp <- reshape2::melt(slices_logp, id.var=c('rsid', 'a_0', 'a_1'),
                       value.name='logp',
                       variable.name='slice')

slices <- cbind(beta, p=10^(-logp$logp), se=se$se)
slices$rsid <- as.character(slices$rsid)
slices$a_0 <- as.character(slices$a_0)
slices$a_1 <- as.character(slices$a_1)
slices$sig <- 2
slices$sig[slices$p > 5*10^(-8)/Teff] <- 1
slices$sig <- factor(slices$sig, levels=c(1,2))
slices$name <- sapply(1:nrow(slices), function(x) {
    paste(slices$pos[x], slices$rsid[x], sep="-")
})
slices$slice_numeric <- gsub("Slice_", "", slices$slice)

## Effect size estimates per variant across slices ####
# Create a custom color scale
myColors <- RColorBrewer::brewer.pal(3,"Set1")[1:2]
names(myColors) <- levels(slices$sig)
colScale <- scale_colour_manual(name = "GWAS", values = myColors,
                                labels = c(expression(p<8%*%10^-9),
                                           expression(p>=8%*%10^-9)))

# generate legend
p_legend <- ggplot(slices, aes(x=slice_numeric, y=beta, color=sig))
p_legend <- p_legend + geom_point() +
    colScale +
    theme_bw()
sigLegend <- cowplot::get_legend(p_legend)

slices_per_chr <- split(slices, f=factor(slices$chr,levels=unique(slices$chr)))
perChr_beta <- lapply(slices_per_chr, function(chr) {
    chr <- chr[order(chr$pos),]
    slices_per_pos <- split(chr,
                            f=factor(chr$pos, levels=unique(chr$pos)))
    p_pos <- lapply(slices_per_pos, function(pos) {
        pos$sig <- factor(as.numeric(pos$sig), levels=c(2,1))
        p_beta <- ggplot(pos, aes(x=slice_numeric, y=beta, color=sig))
        p_beta <- p_beta + geom_point() +
            facet_wrap(~name) +
            ylim(c(-0.16,0.16)) +
            geom_hline(yintercept=0, color='grey') +
        colScale +
            xlab('Slice') +
            ylab(expression(beta~estimates)) +
            theme_bw() + theme(strip.background=element_rect(fill='white'),
                               legend.position = "none")
        return(p_beta)
    })
    null_plots <- 3 - length(unique(chr$pos))
    if (null_plots == 0) {
        p <- cowplot::plot_grid(plotlist=p_pos, ncol=3)
    }
    if (null_plots == 1) {
        null_list <- vector('list', length=null_plots)
        p <- cowplot::plot_grid(p_pos[[1]], p_pos[[2]], NULL, ncol=3)
    }
    if (null_plots == 2) {
        null_list <- vector('list', length=null_plots)
        p <- cowplot::plot_grid(p_pos[[1]], NULL, NULL, ncol=3)
    }
    return(p)
})

p_plots <- cowplot::plot_grid(plotlist=perChr_beta, nrow=length(slices_per_chr),
                        align='v', axis='l', hjust=0, vjust=1,
                        labels=paste('chr', names(slices_per_chr), sep=""))
p_all <- cowplot::plot_grid(sigLegend, p_plots, rel_heights=c(1,20), nrow=2)
ggsave(plot=p_all, paste(directory, "/gwas/Distribution_slices_beta.pdf", sep=""),
       height=15, width=6)

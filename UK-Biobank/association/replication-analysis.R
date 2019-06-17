###############################
### Libraries and functions ###
###############################

options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

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
    args$directory <-"/homes/hannah/data/ukbb/ukb-hrt/gwas/190402_fractal_dimension_26k"
    args$sloci <-"/homes/hannah/data/ukbb/ukb-hrt/gwas/190402_fractal_dimension_26k/"
    args$nloci <- 16
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
                                                threshold = 0.05/args$nloci,
                                                matchBy="SNP", tags_ID="SNP",
                                                tagsSplit="|")

write.table(replication_sig_ld$sig,
            file.path(args$directory,
                      "Replication_association_significant_in_discovery.txt"),
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(replication_sig_ld$sig_wo_ld,
            file.path(args$directory,
                      "Replication_association_significant_in_discovery_ldFiltered.txt"),
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

###############################
## Libraries and functions ####
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('optparse')
plots <- modules::import("plots")


################
## analysis ####
################

## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path to directory with plink association
               results [default: %default].", default=NULL),
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
    args$directory <-"/homes/hannah/data/genotype/control"
    args$name <-"HVOL.gencall.combined.clean.related_batch"
    args$verbose <- TRUE
}
directory <- args$directory
name <- args$name
verbose <- args$verbose

## read plink results ####
hvol <- data.table::fread(paste(directory, "/", name, ".assoc.linear", sep=""),
                          data.table=FALSE, stringsAsFactors=FALSE)
hvol <- hvol[!is.na(hvol$P),]
hvol$P[hvol$P ==1] <- 0.99999999

sig <- hvol[hvol$P < 5e-8,]
write.table(sig, paste(directory, "/", name, ".sig.txt", sep=""),
           col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(sig$SNP, paste(directory, "/", name, ".sig.txt", sep=""),
           col.names=FALSE, row.names=FALSE, quote=FALSE)

## plot results ####
if (verbose) message("Plot association results")
p_manhattan <- plots$manhattan(d=hvol,
                                  title="HVOL GWAS batch",
                                  size.x.labels=12, size.y.labels=12,
                                  is.negLog=FALSE, color=c("#fc8d59", "#b30000"))

ggplot2::ggsave(plot=p_manhattan, height=4, width=10,
                file=paste(directory,"/plink_", name, "_lm_manhattanplot.pdf",
                           sep=""))

p_qq <- plots$qqplot(hvol$P, size.text=14, size.title=14)
ggplot2::ggsave(plot=p_qq, height=7, width=7,
                file=paste(directory,"/plink_", name, "_lm_qqplot.pdf", sep=""))

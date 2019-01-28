###############################
## Libraries and functions ####
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('optparse')
bgenie <- modules::import("bgenieResults")
plots <- modules::import("plots")


################
## analysis ####
################

## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path to directory with bgenie association
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

## read bgenie results ####
genomewide <- lapply(1:22, bgenie$readBgenieOutput, directory=directory,
                     name=paste("bgenie_", name, "_lm_st_chr", sep=""),
                     maf=0.01, biallelicOnly=FALSE)
genomewide <- do.call(rbind, genomewide)

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))
index_beta <- which(grepl("beta", colnames(genomewide)))
index_t <- which(grepl("_t", colnames(genomewide)))

## write combined results ####
if (verbose) message("Write file with genome-wide association results")
write.table(genomewide,
            file=paste(directory, "/bgenie_", name, "_lm_st_genomewide.csv",
                       sep=""),
            sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

## plot results ####
if (verbose) message("Plot association results")
p_manhattan <- plots$manhattan(d=genomewide,
                                  chr='chr', bp='pos', p='batch-log10p',
                                  title="HVOL GWAS batch",
                                  size.x.labels=12, size.y.labels=12,
                                  is.negLog=TRUE, color=c("#fc8d59", "#b30000"))

ggplot2::ggsave(plot=p_manhattan, height=4, width=10,
                file=paste(directory,"/bgenie_", name, "_lm_manhattanplot.pdf",
                           sep=""))

p_qq <- plots$qqplot(genomewide[,index_logp], size.text=14, size.title=14)
ggplot2::ggsave(plot=p_qq, height=7, width=7,
                file=paste(directory,"/bgenie_", name, "_lm_qqplot.pdf", sep=""))

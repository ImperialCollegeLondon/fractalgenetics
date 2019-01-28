#################
## libraries ####
#################

options(bitmapType = 'cairo', device = 'pdf')
modules::import_package('ggplot2', attach=TRUE)
optparse <- modules::import_package('optparse')

############
## data ####
############

## read command line args ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
                        dest="directory", type="character",
                        help="Path to directory with imputation QC results
                        [default: %default].", default=NULL),
    optparse$make_option(c("--batch1"), action="store",
                        dest="batch1", type="character",
                        help="Name of first dataset [default: %default].",
                        default=NULL)),
    optparse$make_option(c("--batch2"), action="store",
                        dest="batch2", type="character",
                        help="Name of second dataset [default: %default].",
                        default=NULL)),
    optparse$make_option(c("--batch3"), action="store",
                        dest="batch3", type="character",
                        help="Name of third dataset [default: %default].",
                        default=NULL)),
    optparse$make_option(c("-n", "--name"), action="store",
                        dest="alg", type="character",
                        help="Name of genotype caller [default: %default].",
                        default=NULL)),
    optparse$make_option(c("--suffix"), action="store",
                        dest="suffix", type="character",
                        help="Common suffix across datasets [default: %default].",
                        default="clean.related")),
    optparse$make_option(c("--bgen"), action="store",
                        dest="bgen", type="character",
                        help="Path to bgen sample file[default: %default].",
                        default=NULL)),
    optparse$make_option(c("--HVOL"), action="store",
                        dest="HVOL", type="character",
                        help="Path to file with HVOL samples [default: %default].",
                        default=NULL),
    optparse$make_option(c("--debug"), action="store_true",
                        dest="debug", default=FALSE, type="logical",
                        help="If set, predefined arguments are used to test the
                        script [default: %default]."))

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    directory <- "~/data/genotype"
    alg <- "gencall"
    batch1 <- "sanger12"
    batch2 <- "singapore12"
    batch3 <- "singapore3"
    suffix <- "clean.related"
    bgen <- "~/data/genotype/imputation/combined/genotypes/gencall.combined.clean.related.chr1.sample"
    HVOL <- "~/data/genotype/QC/combined/European.HVOL.gencall.combined.txt"
}

directory <- args$directory
alg <- args$alg
batch1 <- args$batch1
batch2 <- args$batch2
batch3 <- args$batch3
ukb <- args$ukb
HVOL <- args$HVOL
bgen <- args$bgen
suffix <- args$suffix

data1 <- data.table::fread(paste(directory, "/QC/", batch1, "/", alg, ".",
                                 batch1, ".fam", sep=""),
                           header=FALSE, stringsAsFactors=FALSE,
                           data.table=FALSE)

data2 <- data.table::fread(paste(directory, "/QC/", batch2, "/", alg, ".",
                                 batch2,".fam", sep=""),
                           header=FALSE, stringsAsFactors=FALSE,
                           data.table=FALSE)

data3 <- data.table::fread(paste(directory, "/QC/", batch3, "/", alg, ".",
                                 batch3, ".fam", sep=""),
                           header=FALSE, stringsAsFactors=FALSE,
                           data.table=FALSE)

combined <- data.table::fread(paste(directory, "/QC/combined/", alg,
                                    ".combined.", suffix, ".fam", sep=""),
                           header=FALSE, stringsAsFactors=FALSE,
                           data.table=FALSE)

allsamples <- data.table::fread(bgen, stringsAsFactors=FALSE, data.table=FALSE)
allsamples <- allsamples[-1,]

HVOLsamples <- data.table::fread(HVOL, stringsAsFactors=FALSE, data.table=FALSE)
name <- paste("HVOL", ".", alg, ".combined.", suffix, sep="")

################
## analysis ####
################

## create and format data for BGENIE GWAS against batch ####
colnames(combined) <- c("FID", "IID", "PID","MID","SEX", "PHENO")
combined$PHENO <- 1
combined$PHENO[combined[,2] %in% data2[,2]] <- 2
combined$PHENO[combined[,2] %in% data3[,2]] <- 3

write.table(combined,
            paste(directory, "/control/ALL.", alg, ".combined.", suffix,
                  ".batchPhenotypes.fam", sep=""),
            sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)

combined$PHENO[!(combined[,1] %in% HVOLsamples[,2])] <- "-9"
write.table(combined,
            paste(directory, "/control/HVOL.", alg, ".combined.", suffix,
                  ".batchPhenotypes.fam", sep=""),
            sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)

bgenie <- merge(allsamples, combined, by=1, all.x=TRUE, sort=FALSE)
bgenie <- bgenie[match(allsamples[,1], bgenie[,1]),]
bgenie$PHENO[!(bgenie[,1] %in% HVOLsamples[,2])] <- "-999"

write.table(bgenie$PHENO,
            paste(directory, "/control/Batch_phenotypes_", name, "_bgenie.txt",
                  sep=""),
            sep=" ", row.names=FALSE, col.names="batch", quote=FALSE)

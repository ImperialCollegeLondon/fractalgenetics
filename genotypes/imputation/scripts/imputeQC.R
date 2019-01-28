#################
### libraries ###
#################

options(bitmapType = 'cairo', device = 'pdf')
modules::import_package('ggplot2', attach=TRUE)
optparse <- modules::import_package('optparse')

################
## analysis ####
################
## read command line args ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
                        dest="directory", type="character",
                        help="Path to directory with imputation QC results
                        [default: %default].", default=NULL),
    optparse$make_option(c("-n", "--name"), action="store",
                        dest="name", type="character",
                        help="Name of dataset [default: %default].",
                        default=NULL))
args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))
directory <- args$directory
name <- args$name
if (!is.null(name)) name <- paste(name, ".", sep="")

#directory <- "/Users/hannah/data/genotype/imputation/combined/counts"

## Display maximum posterior probability of imputation per chunk per chr ####
concordance_files = dir(directory, pattern=".chunkConcordance",
                        full.names=TRUE, ignore.case=TRUE)
stats_files = dir(directory, pattern=".chunkStats", full.names=TRUE,
                  ignore.case=TRUE)

pdf(paste(directory, "/", name, "imputeQC.perChrMPP.pdf", sep=""), width = 12,
    height = 6)
concordance <- lapply(seq_along(stats_files), function(chr) {
    concordance <- read.table(concordance_files[chr], header=TRUE, sep="\t",
                              stringsAsFactors=FALSE)
    stats <- read.table(stats_files[chr], header=TRUE, sep="\t",
                        stringsAsFactors=FALSE)
    concordance$MPP <- as.factor(paste("[",concordance$Start, "-",
                                       concordance$End, "]", sep=""))
    concordance_list <- split(concordance, f=concordance$Chunk)
    concordance_list <- lapply(concordance_list, function(x) {
        x$total <- sum(x$Genotypes)
        x$percent <- x$Genotypes/x$total*100
        return(x)
    })
    concordance <- do.call(rbind, concordance_list)
    upper <- dplyr::filter_(concordance, ~MPP %in%
                                c("[0.7-0.8]", "[0.8-0.9]", "[0.9-1]"))
    r  <- ggplot()
    r <- r + geom_bar(data=upper, aes_string(x="Chunk", y="percent",
                                             fill="MPP"),
                      stat="identity", position="dodge") +
        scale_fill_manual(values=c('#edf8b1','#7fcdbb','#2c7fb8')) +
        geom_point(data=stats, aes_string(x="Chunk", y="concordance_overall"),
                   color="red" ) +
        annotate("text", label="Cumulative concordance [%]", x=0, y=104,
                 color="red", hjust=0) +
        xlab(paste("Chromosome", chr, "[chunks]", sep=" ")) +
        ylab("Number of SNPs [%]") +
        labs(title=paste("Percentage of SNPs with maximum posterior",
                         "probability (MPP) greater than 70%")) +
        theme_bw()
    print(r)
    return(concordance)
})
dev.off()

## overview of SNP counts after each genotyping and imputation step ####
overview <- read.table(paste(directory, "/", name, "SNPsPerChr.txt", sep=""),
                       sep="\t", header=TRUE,stringsAsFactors=FALSE)
chromosomes <- unique(overview$Chr)
overview$Chr[overview$Chr == "X"] <- 23
overview$Chr <- as.numeric(overview$Chr)
overview <- reshape2::melt(overview, id.vars="Chr", variable.name="type",
                           value.name="SNP")


pdf(paste(directory, "/", name, "SNPsPerChr.pdf", sep=""), onefile = TRUE,
    paper = "a4r", width = 14, height = 6)
p <- ggplot(data=overview, aes_string(x="Chr", y="SNP", fill="type"))
p + geom_bar(stat="identity", position="dodge") +
    scale_x_continuous(breaks=unique(overview$Chr),
                       labels=paste("chr", chromosomes, sep="")) +
    scale_fill_manual(values=c('#41b6c4','#2c7fb8','#253494')) +
    xlab("Chromosomes") +
    ylab("Number of SNPs") +
    labs(title ="SNP numbers after genotyping, imputation and imputation QC") +
    theme_bw() +
    theme(legend.position="bottom",
          axis.text.x =element_text(angle=90, hjust=1))
dev.off()

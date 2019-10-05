#################
## libraries ####
#################

modules::import_package('ggplot2', attach=TRUE)
modules::import_package('tidyr', attach=TRUE)
modules::import_package('optparse', attach=TRUE)


#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-d", "--dir"), action="store", dest="dir",
                type="character", help="Path to directory with expression atlas
                data [default: %default].", default=NULL),
    make_option(c("-p", "--prefix"), action="store", dest="prefix",
                type="character", help="Common prefix of expression atlas files
                [default: %default].", default=NULL),
    make_option(c("--debug"), action="store_true",
                dest="debug", default=FALSE, type="logical",
                help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- optparse$parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$dir <- "~/data/ukbb/GTEX"
    args$prefix <- "expression_atlas-homo_sapiens"
}

#############
## data  ####
#############

## expression values ####
genefiles <- list.files(path=args$dir, pattern=args$prefix, full.names=TRUE)
gene_order <- c("chr1:PKP1",  "chr3:PDZRN3", "chr5:SAP30L", "chr6:C4A",
                "chr6:SSXP10",   "chr8:TDH", "chr8:MTSS1", "chr17:GOSR2",
                "chr19:ZNF358")
names(gene_order) <- gsub("chr.*:(.*)", "\\1", gene_order)

tpm <- lapply(genefiles, function(s) {
    g <- gsub(paste(".*", args$prefix, "_(.*)\\.tsv", sep=""), "\\1", s)
    g <- gene_order[names(gene_order) == g]
    dat <- data.table::fread(s, header=TRUE, data.table=FALSE,
                             stringsAsFactors=FALSE, na.strings=c("NA", "NaN"))
    consortium <- dat$V1
    dat <- data.frame(t(dat[,-1]))
    colnames(dat) <- consortium
    dat$tissue <- rownames(dat)
    dat <- dat %>% pivot_longer(-tissue, names_to="Consortium",
                                values_to = "count")
    dat$gene <- rep(g, nrow(dat))
    return(dat)
})

tpm_all <- do.call(rbind, tpm)
tpm_all$gene <- factor(tpm_all$gene, levels=rev(gene_order))

################
## analysis ####
################

p <- ggplot(data=filter(tpm_all, Consortium=="GTEx"))
p <- p + geom_tile(aes(x=tissue, y=gene, fill=log10(count+1)), width=0.9,
                   height=0.99, color="white") +
    scale_fill_distiller(type='seq', palette = 2, direction=1,
                         name="log10(TPM)",
                         na.value ="lightgrey") +
    labs(x="GTEx tissues", y="GTEx target genes") +
    theme_bw() +
    scale_y_discrete(expand=c(0,0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5,
                                     size = 8)) +
    theme(axis.text.y = element_text(hjust = 1, size = 10)) +
    theme(axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank())

ggsave(plot=p, file=paste(args$dir, "/gTEX_geneexpression.pdf", sep=""),
       height=4, width=10, units="in")

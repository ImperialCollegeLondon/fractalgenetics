#################
## libraries ####
#################
options(import.path=c("~/analysis/fractalgenetics",
                      "~/projects"))
options(bitmapType = 'cairo', device = 'pdf')


modules::import_package('ggplot2', attach=TRUE)
optparse <- modules::import_package('optparse')
autofd <- modules::import('fractal-analysis-processing')

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-d", "--dir"), action="store", dest="dir",
                type="character", help="Path to output directory
                data [default: %default].", default=NULL),
    make_option(c("--ukb"), action="store", dest="ukb",
                type="character", help="Path to file with UKB per slice FD
                [default: %default].", default=NULL),
    make_option(c("-i", "--interpolate"), action="store", dest="interpolate",
                type="integer", help="Number of slices to interpolate to
               [default: %default].", default=9),
    make_option(c("--debug"), action="store_true",
                dest="debug", default=FALSE, type="logical",
                help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- optparse$parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$dir <- "~/data/ukbb/ukb-hrt/phenotypes"
    args$ukb <- "~/data/ukbb/ukb-hrt/phenotypes/180628_fractal_dimension/FD_slices_EUnorel.csv"
    args$interpolate <- 9
}

#############
## data  ####
#############

## ukb FD measurments ####
ukb <-  data.table::fread(args$ukb, data.table=FALSE)[,-1]
colnames(ukb) <- gsub("Slice_", "", colnames(ukb))

################
## analysis ####
################

pca_UKB <- prcomp(ukb)
pb_12 <- ggbiplot(pca_UKB, labels=NULL, choices=1:2, var.axes=TRUE, alpha=0.1,
                  color.axes = "blue", color.axes.text = "blue",
                  color.points="grey", repel=TRUE) +
    theme_bw()

pb_34 <- ggbiplot(pca_UKB, labels=NULL, choices=3:4, var.axes=TRUE, alpha=0.1,
                  color.axes = "blue", color.axes.text = "blue",
                  color.points="grey", repel=TRUE) +
    theme_bw()

ps <- ggscreeplot(pca_UKB) + theme_bw() +
    theme(panel.grid.minor =element_blank())

all <- egg::ggarrange(ps, pb_12, pb_34, nrow=1, labels=c("A", "B", "C"),
                       label.args = list(gp = grid::gpar(font = 2, cex =
                                                             1.2)))
ggsave(plot=all,
       paste(args$dir, "/UKB_slices", args$interpolate, "pca.pdf", sep=""),
       width=18, height=6)

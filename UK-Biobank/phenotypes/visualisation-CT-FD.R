#################
## libraries ####
#################
options(bitmapType = 'cairo', device = 'pdf')


modules::import_package('tidyverse', attach=TRUE)
modules::import_package('optparse', attach=TRUE)

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-d", "--dir"), action="store", dest="dir",
                type="character", help="Path to output directory
                data [default: %default].", default=NULL),
    make_option(c("--ct"), action="store", dest="ct",
                type="character", help="Path to file with CT per slice FD
                [default: %default].", default=NULL),
    make_option(c("--debug"), action="store_true",
                dest="debug", default=FALSE, type="logical",
                help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- optparse$parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$dir <- "~/data/ukbb/ukb-hrt/phenotypes"
    args$ct <- "~/data/ukbb/ukb-hrt/rawdata/ct_fd.csv"
}
#############
## data  ####
#############

## CT FD measurements ####
FDalongHeart <-  read_csv(file.path(args$ct))
colnames(FDalongHeart) <- c("ID", "Slice", "FD")
FDalongHeart <- FDalongHeart %>% drop_na()

################
## analysis ####
################

# plot distribution of FD along heart
FDalongHeart$Slice <- as.factor(FDalongHeart$Slice)
FDalongHeart$Location <- "Apical section"
FDalongHeart$Location[as.numeric(FDalongHeart$Slice) <= 3] <- "Basal section"
FDalongHeart$Location[as.numeric(FDalongHeart$Slice) <= 6 &
                          as.numeric(FDalongHeart$Slice) > 3] <- "Mid section"
FDalongHeart$Location <- factor(FDalongHeart$Location,
                               levels=c("Basal section", "Mid section",
                                         "Apical section"))


p_fd <- ggplot(data=FDalongHeart, aes(x=Slice, y=FD))
p_fd <- p_fd +
    geom_boxplot(aes(color=Location), outlier.size = 0.2) +
    scale_color_manual(values=c('#67a9cf','#1c9099','#016c59')) +
    labs(x="Slice", y="FD") +
    theme_bw()

ggsave(plot=p_fd,
       file=paste(args$dir, "/FDAlongHeart_CT_all_slices",
                  args$interpolate, ".pdf", sep=""),
       height=2, width=5, units="in")



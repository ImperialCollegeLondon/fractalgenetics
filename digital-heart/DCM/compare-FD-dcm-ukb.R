#################
## libraries ####
#################
options(import.path=c("~/analysis/fractalgenetics",
                      "~/projects"))
options(bitmapType = 'cairo', device = 'pdf')


modules::import_package('tidyverse', attach=TRUE)
modules::import_package('dplyr', attach=TRUE)
modules::import_package('ggplot2', attach=TRUE)
modules::import_package('optparse', attach=TRUE)
autofd <- modules::import('fractal-analysis-processing')

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-d", "--dir"), action="store", dest="dir",
                type="character", help="Path to output directory
                data [default: %default].", default=NULL),
    make_option(c("--dcm"), action="store", dest="dcm",
                type="character", help="Path to file with DCM per slice FD
                [default: %default].", default=NULL),
    make_option(c("-i", "--interpolate"), action="store", dest="interpolate",
                type="integer", help="Number of slices to interpolate to
               [default: %default].", default=9),
    make_option(c("--ukb"), action="store", dest="ukb",
                type="character", help="Path to file with UKB per slice FD
                [default: %default].", default=NULL),
    make_option(c("--debug"), action="store_true",
                dest="debug", default=FALSE, type="logical",
                help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- optparse$parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$dir <- "~/data/ukbb/ukb-hrt/DCM"
    args$dcm <- "~/data/digital-heart/phenotype/FD/FD_slices_EUnorel.csv"
    args$ukb <- "~/data/ukbb/ukb-hrt/phenotypes/180628_fractal_dimension/FD_slices_EUnorel.csv"
    args$interpolate <- 9
}

#############
## data  ####
#############

## DCM FD measurements ####
DCM_slices <-  data.table::fread(args$dcm, data.table=FALSE)[,-1]

## ukb FD measurments ####
UKB_slices <-  data.table::fread(args$ukb, data.table=FALSE)[,-1]

################
## analysis ####
################

## combine datasets ####
DCM_slices$study <- "DCM"
UKB_slices$study <- "UKB"
combined <- rbind(DCM_slices, UKB_slices)
saveRDS(combined, file.path(args$outdir, "ukb_DCM_FD.rds"))


# plot distribution of FD along heart
FDalongHeart <- tidyr::pivot_longer(combined,  cols=starts_with("Slice") ,
                                    values_to="FD", names_to="Slice")

FDalongHeart$Slice <- as.factor(as.numeric(gsub("Slice_", "",
                                                FDalongHeart$Slice)))
FDalongHeart$Location <- "Apical section"
FDalongHeart$Location[as.numeric(FDalongHeart$Slice) <= 3] <- "Basal section"
FDalongHeart$Location[as.numeric(FDalongHeart$Slice) <= 6 &
                          as.numeric(FDalongHeart$Slice) > 3] <- "Mid section"
FDalongHeart$Location <- factor(FDalongHeart$Location,
                               levels=c("Basal section", "Mid section",
                                         "Apical section"))
FDalongHeart$study <- as.factor(FDalongHeart$study)

p_fd <- ggplot(data=FDalongHeart, aes(x=Slice, y=FD, fill=study,))
p_fd <- p_fd +
    geom_boxplot(aes(color=Location), outlier.size = 0.2) +
    scale_color_manual(values=c('#67a9cf','#1c9099','#016c59')) +
    scale_fill_manual(values=c('#fc8d62', '#8da0cb'), name="Cohort") +
    labs(x="Slice", y="FD") +
    theme_bw()

ggsave(plot=p_fd,
       file=paste(args$dir, "/FDAlongHeart_DCM_UKB_all_slices",
                  args$interpolate, ".pdf", sep=""),
       height=2, width=5, units="in")

lmm_study <- nlme::lme(FD ~ study, random=~1|Slice, data=FDalongHeart)
p.value  <- 2*pt(-abs(-53.48570), df=165617, log=TRUE)

tTable <- summary(lmm_study)$tTable %>%
    as_tibble %>%
    select(-`p-value`)

tTable$log_pvalues <- 2*pt(-abs(tTable$`t-value`), df=tTable$DF, log=TRUE)
write.table(tTable,
            file.path(args$dir, "LMM_fixed_study_random_slices_tTable.csv"),
            sep=",", col.names=TRUE, row.names=FALSE)


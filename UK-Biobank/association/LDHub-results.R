#################
## libraries ####
#################

library(ggplot2)
library(dplyr)

getHeritability <- function(fn, pheno) {
    if (file.exists(fn)) {
        tmp <- readLines(fn)
        h2_string <- tmp[grepl("Total Observed scale h2:", tmp)]
        h2 <- as.numeric(gsub(".* ([0-9\\.]*) \\(([0-9\\.]*)\\)",
                              "\\1", h2_string))
        h2_se <- as.numeric(gsub(".* ([0-9\\.]*) \\(([0-9\\.]*)\\)",
                                 "\\2", h2_string))
        lambda_string <- tmp[grepl("Lambda GC:",tmp)]
        lambda_gc <- as.numeric(gsub(".* ([0-9\\.]*)", "\\1", lambda_string))
        chi2_string <- tmp[grepl("Mean Chi",tmp)]
        chi2 <- as.numeric(gsub(".* ([0-9\\.]*)", "\\1", chi2_string))
        intercept_string <- tmp[grepl("Intercept",tmp)]
        intercept <- as.numeric(gsub(".* ([0-9\\.]*) \\(([0-9\\.]*)\\)",
                                     "\\1", intercept_string))
        intercept_se <- as.numeric(gsub(".* ([0-9\\.]*) \\(([0-9\\.]*)\\)",
                                        "\\2", intercept_string))
        ratio_string <- tmp[grepl("Ratio",tmp)]
        if (ratio_string == "Ratio < 0 (usually indicates GC correction).") {
            ratio <- NA
            ratio_se <- NA
            ratio_character <- "Ratio < 0"
            ratio_se_character <- "-"
        } else {
            ratio <- as.numeric(gsub(".* ([0-9\\.]*) \\(([0-9\\.]*)\\)",
                                     "\\1", ratio_string))
            ratio_se <- as.numeric(gsub(".* ([0-9\\.]*) \\(([0-9\\.]*)\\)",
                                        "\\2", ratio_string))
            ratio_character <- as.character(ratio)
            ratio_se_character <- as.character(ratio_se)
        }
        h2=data.frame(pheno=pheno, h2=h2, h2_se=h2_se, chi2=chi2,
                      lambda=lambda_gc, ratio=ratio, ratio_se=ratio_se,
                      intercept=intercept, intercept_se=intercept_se,
                      ratio_character=ratio_character,
                      ratio_se_character=ratio_se_character)
        return(h2)
    }
}

getGeneticCorrelation <- function(fn, pheno) {
    if (file.exists(fn)) {
        tmp <- read.table(fn, stringsAsFactors=FALSE,
                          sep=",", quote="", header=TRUE, fill=TRUE)
        # remove NA results
        tmp <- tmp[!is.na(tmp$p),]
        # use only European summary stats
        tmp <- tmp[tmp$ethnicity == "European",]
        # order by category
        tmp <- tmp[order(tmp$Category),]
        tmp$Category <- as.factor(tmp$Category)
        # Numeric traits position
        tmp$pos <- 1:nrow(tmp)
        return(tmp)
    }
}



############
## data ####
############
# LDscore regression was conducted using LDhub at http://ldsc.broadinstitute.org
# LDhub collected summary statistics for 732 traits, including 516 ukbb traits.
# Formated, summary stats data (via association.smk and
# association/association-results.R were uploaded and all LDhub traits chosen
# for genetic correlation analyses. The application returns a log file,
# containing heritability estimates and .rg.txt files containing the genetic
# correlation estimates.
# These files are read below.

directory <-"~/data/ukbb/ukb-hrt/ldhub"
interpolate <- 9
regions <- c("MeanApicalFD", "MeanMidFD", "MeanBasalFD")

# from heritability.log
heritability_slices <- lapply(1:interpolate, function(pheno) {
    fn <- paste(directory, "/sumstats_slices_Slice_", pheno, "-h2.log", sep="")
    getHeritability(fn, pheno)
})
h2_slices <- do.call(rbind, heritability_slices)
write.table(h2_slices, file.path(directory,
                                 "LDscore_heritability_estimates_slices.csv"),
            sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

heritability_summary <- lapply(regions, function(pheno) {
    fn <- paste(directory, "/sumstats_summary_", pheno, "-h2.log", sep="")
    getHeritability(fn, pheno)
})
h2_summary <- do.call(rbind, heritability_summary)
write.table(h2_summary, file.path(directory,
                                  "LDscore_heritability_estimates_summary.csv"),
            sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

# from rg.csv
rg_ld_slices <- lapply(1:interpolate, function(slice) {
    fn <- paste(directory, "/sumstats_slices_Slice_", slice, "_rg.results.csv",
                sep="")
    getGeneticCorrelation(fn)
})
rg_ld_slices <- do.call(rbind, rg_ld_slices)

rg_ld_summary <- lapply(regions, function(pheno) {
    fn <- paste(directory, "/sumstats_summary_", pheno, "_rg.results.csv",
                sep="")
    getGeneticCorrelation(fn)
})
rg_ld_summary <- do.call(rbind, rg_ld_summary)
rg_ld_summary$trait1 <- factor(gsub("Mean(.*)FD", "\\1", rg_ld_summary$trait1),
                               levels=c("Basal", "Mid", "Apical"))

# Manual selection of any heart/cardiovascular phenotypes from ldhub
ukbb_names <- c(
    "Diastolic blood pressure_ automated reading",
    "Systolic blood pressure_ automated reading",
    "Pulse rate",
    "Pulse wave reflection index",
    "Pulse wave peak to peak time",
    "Target heart rate achieved",
    "Pulse wave Arterial Stiffness index",
    "Non-cancer illness code_ self-reported: hypertension",
    "Non-cancer illness code_ self-reported: angina",
    "Non-cancer illness code_ self-reported: heart attack/myocardial infarction",
    "Non-cancer illness code_ self-reported: hypertrophic cardiomyopathy (hcm / hocm)",
    "Illnesses of father: Heart disease",
    "Illnesses of father: High blood pressure",
    "Illnesses of father: Diabetes",
    "Illnesses of mother: Heart disease",
    "Illnesses of mother: High blood pressure",
    "Illnesses of siblings: Heart disease",
    "Illnesses of siblings: High blood pressure",
    "Vascular/heart problems diagnosed by doctor: Heart attack",
    "Vascular/heart problems diagnosed by doctor: None of the above",
    "Vascular/heart problems diagnosed by doctor: Angina",
    "Vascular/heart problems diagnosed by doctor: High blood pressure",
    "Medication for cholesterol_ blood pressure_ diabetes_ or take exogenous hormones: Blood pressure medication",
    "Medication for cholesterol_ blood pressure or diabetes: Blood pressure medication",
    "Diagnoses - main ICD10: I20 Angina pectoris",
    "Diagnoses - main ICD10: I21 Acute myocardial infarction",
    "Diagnoses - main ICD10: I25 Chronic ischaemic heart disease",
    "Diagnoses - main ICD10: I30 Acute pericarditis",
    "Diagnoses - main ICD10: I48 Atrial fibrillation and flutter"
)
short_names <- c("Diastolic blood pressure",
                "Systolic blood pressure",
                "Pulse rate",
                "Pulse wave reflection index",
                "Pulse wave peak to peak time",
                "Target heart rate achieved",
                "Pulse wave Arterial Stiffness Index",
                "self-reported: hypertension",
                "self-reported: angina",
                "self-reported: heart attack/myocardial infarction",
                "self-reported: hypertrophic cardiomyopathy",
                "Illnesses of father: Heart disease",
                "Illnesses of father: High blood pressure",
                "Illnesses of father: Diabetes",
                "Illnesses of mother: Heart disease",
                "Illnesses of mother: High blood pressure",
                "Illnesses of siblings: Heart disease",
                "Illnesses of siblings: High blood pressure",
                "diagnosed by doctor: Heart attack",
                "diagnosed by doctor: Vascular/heart problems",
                "diagnosed by doctor: Angina",
                "diagnosed by doctor: High blood pressure",
                "Medication (1) : Blood pressure medication",
                "Medication (2): Blood pressure medication",
                "ICD10: I20 Angina pectoris",
                "ICD10: I21 Acute myocardial infarction",
                "ICD10: I25 Chronic ischaemic heart disease",
                "ICD10: I30 Acute pericarditis",
                "ICD10: I48 Atrial fibrillation and flutter")

heart_pheno <- data.frame(ukbb_names=ukbb_names, names=short_names,
                          stringsAsFactors=FALSE)

write.table(heart_pheno,
            file=file.path(directory, "ldhub_heart_phenotypes.csv"),
            sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

################
## analysis ####
################

## Depict heritability estimates per slice ####
h2_slices$region <- "Mid"
h2_slices$region[h2_slices$pheno <=3] <- "Basal"
h2_slices$region[h2_slices$pheno >=7] <- "Apical"
h2_slices$region <- factor(h2_slices$region, levels=c("Basal", "Mid", "Apical"))

h_slices <- ggplot(data=h2_slices, aes(x=pheno, y=h2, color=region))
h_slices <- h_slices + geom_point() +
    geom_pointrange(aes(ymin=h2-h2_se, ymax=h2+h2_se)) +
    scale_x_continuous(labels=1:interpolate, breaks=1:interpolate) +
    scale_color_manual(values=c('#67a9cf','#1c9099','#016c59'),
                       name="Region") +
    xlab("Slice") +
    ylab("Additive heritability") +
    theme_bw()
ggsave(plot=h_slices, filename=file.path(directory, "heritability_slices.pdf"),
       width=5, height=3)

## Depict genetic correlation estimates in 'manhattan plot';
## any heart-related traits highlighted in magenta
rg_ld_slices$heart <- 0
rg_ld_slices$heart[rg_ld_slices$trait2 %in% heart_pheno$ukbb_name] <- 1
heart_slices <- filter(rg_ld_slices, heart == 1) %>%
    droplevels
heart_slices$name <- heart_pheno$name

rg_slices <- ggplot()
rg_slices <- rg_slices + geom_point(data=rg_ld_slices,
                                    aes(x=pos, y=-log10(p), color=Category),
                    size=0.5) +
    geom_point(data=filter(rg_ld_slices, heart == 1), aes(x=pos, y=-log10(p)),
               color='#dd1c77', size=0.5) +
    scale_x_discrete(labels=NULL) +
    facet_wrap(~trait1, nrow=interpolate) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                '#fb9a99', '#e31a1c','#fdbf6f','#ff7f00',
                                '#cab2d6','#6a3d9a', '#ffff99','#b15928',
                                '#8dd3c7','#ffffb3','#bebada', '#fb8072',
                                '#80b1d3','#fdb462','#b3de69','#fccde5',
                                '#d9d9d9','#bc80bd','#ccebc5','#ffed6f')) +
    xlab("LD Hub traits") +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"))

ggsave(plot=rg_slices, filename=file.path(directory, "Rg_slices_manhattan.pdf"),
       width=8, height=10)

rgh_slices <- ggplot()
rgh_slices <- rgh_slices + geom_point(data=heart_slices,
                                      aes(x=-log10(p), y=name, color=trait1)) +
    scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3','#e7298a',
                                '#66a61e','#e6ab02','#a6761d', "#e41a1c",
                                "#377eb8"),
                       name="Slice") +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"),
          axis.title.y = element_blank())
          #axis.text.x = element_text(angle=45, hjust=1, vjust=1))
ggsave(plot=rgh_slices, filename=file.path(directory, "Rg_slices_heart.pdf"),
       width=8, height=10)


sig_slices <- filter(rg_ld_slices,p < 5*10^-4) %>% droplevels
sig_slices$rg <- as.numeric(sig_slices$rg)

rgsig_slices <- ggplot()
rgsig_slices <- rgsig_slices + geom_point(data=sig_slices, aes(x=rg, y=trait2,
                                          color=Category, size=-log10(p))) +
    geom_point(data=filter(sig_slices, heart == 1), aes(x=rg, y=trait2,
                                                        size=-log10(p)),
               color='#dd1c77') +
    facet_wrap(~trait1, ncol=interpolate) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                '#fb9a99', '#e31a1c','#fdbf6f','#ff7f00',
                                '#cab2d6','#6a3d9a','#ffff99','#b15928',
                                '#8dd3c7','#ffffb3','#bebada', '#fb8072',
                                '#80b1d3','#fdb462','#b3de69','#fccde5',
                                '#d9d9d9','#bc80bd','#ccebc5','#ffed6f')) +
    scale_size(range=c(0.2, 3)) +
    xlab("Genetic correlation") +
    ylab("") +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"),
          legend.position = 'bottom',
          legend.box = "vertical",
          plot.margin = unit(c(5.5,12,5.5,5.5), "pt"))
ggsave(plot=rgsig_slices, filename=file.path(directory,
                                             "Rg_slices_sig_heart.pdf"),
       width=12, height=10)

## Depict heritability estimates per region ####
h2_summary$region <- gsub("Mean(.*)FD", "\\1", h2_summary$pheno)
h2_summary$region <- factor(h2_summary$region,
                            levels=c("Basal", "Mid", "Apical"))

h_summary <- ggplot(data=h2_summary, aes(x=region, y=h2, color=region))
h_summary <- h_summary + geom_point() +
    geom_pointrange(aes(ymin=h2-h2_se, ymax=h2+h2_se)) +
    scale_color_manual(values=c('#67a9cf','#1c9099','#016c59'),
                       guide=FALSE) +
    xlab("Cardiac region") +
    ylim(c(0,0.3)) +
    ylab("Additive heritability") +
    theme_bw()
ggsave(plot=h_summary, filename=file.path(directory, "heritability_summary.pdf"),
       width=2, height=2)

## Depict genetic correlation estimates in 'manhattan plot';
## any heart-related traits highlighted in magenta
rg_ld_summary$heart <- 0
rg_ld_summary$heart[rg_ld_summary$trait2 %in% heart_pheno$ukbb_name] <- 1
heart_summary <- filter(rg_ld_summary, heart == 1) %>%
    droplevels
heart_summary$name <- heart_pheno$name

rg_summary <- ggplot()
rg_summary <- rg_summary + geom_point(data=rg_ld_summary,
                                      aes(x=pos, y=-log10(p), color=Category),
                      size=0.5) +
    geom_point(data=filter(rg_ld_summary, heart == 1), aes(x=pos, y=-log10(p)),
               color='#dd1c77', size=0.5) +
    scale_x_discrete(labels=NULL, expand=c(0.05,0)) +
    facet_wrap(~trait1, ncol=length(regions)) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                '#fb9a99', '#e31a1c','#fdbf6f','#ff7f00',
                                '#cab2d6','#6a3d9a', '#ffff99','#b15928',
                                '#8dd3c7','#ffffb3','#bebada', '#fb8072',
                                '#80b1d3','#fdb462','#b3de69','#fccde5',
                                '#d9d9d9','#bc80bd','#ccebc5','#ffed6f')) +
    xlab("LD Hub traits") +
    guides(color=guide_legend(nrow=4,byrow=TRUE)) +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"),
          legend.position='bottom',
          legend.box = "vertical",
          plot.margin = unit(c(5.5,12,5.5,5.5), "pt"))

ggsave(plot=rg_summary, filename=file.path(directory, "Rg_summary_manhattan.pdf"),
       width=10, height=6)

rgh_summary <- ggplot()
rgh_summary <- rgh_summary +
    geom_point(data=heart_summary,
               aes(x=-log10(p), y=name, color=trait1, size=as.numeric(rg))) +
    scale_color_manual(values=c('#67a9cf','#1c9099','#016c59'),
                       name="Region") +
    scale_size_continuous(name="Genetic correlation", limits=c(-0.45, 0.3)) +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"),
          axis.title.y = element_blank(),
          legend.position="bottom",
          legend.box="vertical")
ggsave(plot=rgh_summary, filename=file.path(directory, "Rg_summary_heart.pdf"),
       width=8, height=6)

sig_heart_summary <- filter(heart_summary,p < 5*10^-2)# %>% droplevels
rgh_sig_summary <- ggplot()
rgh_sig_summary <- rgh_sig_summary +
    geom_point(data=sig_heart_summary,
               aes(x=-log10(p), y=name, color=trait1, size=as.numeric(rg))) +
    scale_color_manual(values=c('#67a9cf','#1c9099','#016c59'),
                       name="Region") +
    scale_size_continuous(name="Genetic correlation", limits=c(-0.3, 0.3)) +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"),
          axis.title.y = element_blank(),
          legend.position="bottom",
          legend.box="vertical")
ggsave(plot=rgh_sig_summary,
       filename=file.path(directory, "Rg_sig_summary_heart.pdf"),
       width=6, height=5)


sig_summary <- filter(rg_ld_summary,p < 5*10^-2) %>% droplevels
sig_summary$rg <- as.numeric(sig_summary$rg)

rgsig_summary <- ggplot()
rgsig_summary <- rgsig_summary + geom_point(data=sig_summary,
                                            aes(x=rg, y=trait2,
                                                color=Category,
                                                size=-log10(p))) +
    geom_point(data=filter(sig_summary, heart == 1), aes(x=rg, y=trait2,
                                                        size=-log10(p)),
               color='#dd1c77') +
    facet_wrap(~trait1, ncol=interpolate) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                '#fb9a99', '#e31a1c','#fdbf6f','#ff7f00',
                                '#cab2d6','#6a3d9a','#ffff99','#b15928',
                                '#8dd3c7','#ffffb3','#bebada', '#fb8072',
                                '#80b1d3','#fdb462','#b3de69','#fccde5',
                                '#d9d9d9','#bc80bd','#ccebc5','#ffed6f')) +
    scale_size(range=c(0.2, 3)) +
    xlab("Genetic correlation") +
    ylab("") +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"),
          legend.position = 'bottom',
          legend.box = "vertical",
          plot.margin = unit(c(5.5,12,5.5,5.5), "pt"))
ggsave(plot=rgsig_summary,
       filename=file.path(directory, "Rg_summary_sig_heart.pdf"),
       width=12, height=10)

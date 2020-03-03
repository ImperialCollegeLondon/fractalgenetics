###############################
### Libraries and functions ###
###############################
options(bitmapType = 'cairo', device = 'pdf')

modules::import_package('tidyverse', attach=TRUE)
modules::import_package('data.table', attach=TRUE)
modules::import_package('ggplot2', attach=TRUE)
modules::import_package('TwoSampleMR', attach=TRUE)
optparse <- modules::import_package('optparse')

## Pierce 2013 American Journal of Epidemiology
Fstat <- function(r2, nsample, ninstruments) {
    if (length(unique(c(length(r2), length(nsample),
                        length(ninstruments)))) != 1) {
        stop("Length of all input variables has to be the same")
    }
    rep1 <- rep(1, length(nsample))
    r2*(nsample-rep1-ninstruments)/((rep1-r2)*ninstruments)
}

## Burgess 2016 Genetic Epidemiology
findlowerlimitF <- function(f, ninstruments, nsample) {
    lambda <- f*ninstruments*(nsample-2)/nsample-ninstruments
    lower <- f - 1
    while (pf(lower, df1=ninstruments, df2=nsample, ncp=lambda) > 0.05) {
        lower <- lower-1
    }
    upper <- lower + 1
    while (abs(pf((lower+upper)/2, df1=ninstruments, df2=nsample,
                  ncp=lambda)-0.05) > 0.0001) {
        if (pf((lower+upper)/2, df1=ninstruments, df2=nsample,
               ncp=lambda) > 0.05) {
            upper = (lower+upper)/2
        }
        if (pf((lower+upper)/2, df1=ninstruments, df2=nsample,
               ncp=lambda) <0.05) {
            lower = (lower+upper)/2
        }
    }
    return((lower+upper)/2)
}

estimateI2 <- function(y,s) {
    k <- length(y)
    w <- 1/s^2
    sum.w <- sum(w)
    mu.hat <- sum(y*w)/sum.w
    Q <- sum(w*(y-mu.hat)^2)
    Isq <- (Q - (k-1))/Q
    Isq <-  max(0,Isq)
    return(Isq)
}

getMRdata <- function(filename_exposure, data_exposure=NULL, outcomes=NULL,
                      access, name='FD', data_outcome=NULL,
                      filename_outcome=NULL,
                      sep = " ", phenotype_col = "Phenotype", snp_col = "SNP",
                      beta_col = "beta", se_col = "se", eaf_col = "eaf",
                      effect_allele_col = "effect_allele", pval_col='pval',
                      other_allele_col = "other_allele") {
    if (is.null(data_exposure)) {
        data_exposure <- read_exposure_data(filename = filename_exposure)
    } else {
        data_exposure <- format_data(data_exposure, type='exposure',
                                    phenotype_col=phenotype_col,
                                    snp_col=snp_col, beta_col=beta_col,
                                    se_col=se_col, eaf_col=eaf_col,
                                    pval_col=pval_col,
                                    effect_allele_col=effect_allele_col,
                                    other_allele_col=other_allele_col)
    }
    if (!is.null(outcomes)) {
        data_outcome <- extract_outcome_data(data_exposure$SNP, outcomes,
                                        proxies=TRUE, rsq=0.8, align_alleles=1,
                                        palindromes=1, maf_threshold=0.3,
                                        access_token=access)
    } else if (!is.null(filename_outcome)) {
        data_outcome <- read_outcome_data(filename=filename_outcome, sep=sep,
                                         phenotype_col=phenotype_col,
                                         snp_col=snp_col, beta_col=beta_col,
                                         se_col=se_col, eaf_col=eaf_col,
                                         effect_allele_col=effect_allele_col,
                                         other_allele_col=other_allele_col)
    }  else if (!is.null(data_outcome)) {
        data_outcome <- format_data(data_outcome, type='outcome',
                                   phenotype_col=phenotype_col,
                                   snp_col=snp_col, beta_col=beta_col,
                                   se_col=se_col, eaf_col=eaf_col,
                                   pval_col=pval_col,
                                   effect_allele_col=effect_allele_col,
                                   other_allele_col=other_allele_col)
    }

    return(list(exposure=data_exposure, outcome=data_outcome))
}

MRanalysis <- function(exposure_dat, outcome_dat,
                       methods=c("mr_egger_regression", "mr_ivw",
                                 "mr_weighted_median",
                                 "mr_weighted_mode"), verbose=TRUE) {
    if (verbose) message("Harmonise exposure and outcome")
    dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

    if (verbose) message("MR regression analysis")
    mr_results <- mr(dat, method_list=methods)

    if (verbose) message("MR heterogeneity analysis")
    het_results <- mr_heterogeneity(dat, method_list=c("mr_egger_regression",
                                                       "mr_ivw"))
    if (verbose) message("MR pleiotropy analysis")
    plei_results <- mr_pleiotropy_test(dat)

    if (verbose) message("MR leave-one-out analysis")
    loo_results <- mr_leaveoneout(dat)

    if (verbose) message("MR directionality analysis")
    directionality_results <- directionality_test(dat)

    if (verbose) message("MR F statistic")
    per_study_F <-
        data.frame(study=mr_results$outcome[!duplicated(mr_results$outcome)])
    per_study_F$ninstruments <- mr_results$nsnp[!duplicated(mr_results$outcome)]
    per_study_F$samplesize.exposure <- unique(dat$samplesize.exposure)
    per_study_F$samplesize.outcome <-
        dat$samplesize.outcome[!duplicated(dat$outcome)]
    per_study_F$r2 <- directionality_results$snp_r2.exposure
    per_study_F$Fstat <- Fstat(per_study_F$r2, per_study_F$samplesize.exposure,
                               per_study_F$ninstruments)
    per_study_F$lowerBound <- apply(as.matrix(per_study_F[,-1]), 1, function(x){
        findlowerlimitF(x[5], x[1], x[2])
    })

    if (verbose) message("MR I^2 analysis")
    weighted_beta <- dat$beta.exposure/dat$se.outcome
    weighted_se <- dat$se.exposure/dat$se.outcome
    I2 <- estimateI2(weighted_beta,weighted_se)
    return(list(dat=dat, mr_results=mr_results, het_results=het_results,
           plei_results=plei_results, loo_results=loo_results,
           directionality_results=directionality_results,
           I2=I2,
           Fstat=per_study_F))
}

MRtables <- function(MRbase_results, pheno) {
    lapply(seq_along(MRbase_results), function(x, pheno) {
        region <- MRbase_results[[x]]
        name_region <- names(MRbase_results)[x]
        write.table(region$mr_results,
                    file.path(mrdir, str_c(name_region, pheno, "MR_results.csv",
                                           sep="_")),
                    col.names=TRUE, row.names=FALSE, sep=',', quote=FALSE)
        write.table(region$plei_results,
                    file.path(mrdir, str_c(name_region, pheno,
                                           "MR_pleiotropy_results.csv",
                                           sep="_")),
                    col.names=TRUE, row.names=FALSE, sep=',', quote=FALSE)
        write.table(region$directionality_results,
                    file.path(mrdir, str_c(name_region, pheno,
                                           "MR_directionality_results.csv",
                                           sep="_")),
                    col.names=TRUE, row.names=FALSE, sep=',', quote=FALSE)

        write.table(region$I2,
                    file.path(mrdir, str_c(name_region, pheno,
                                           "MR_I2_results.csv", sep="_")),
                    col.names=FALSE, row.names=FALSE, sep=',', quote=FALSE)
        write.table(region$Fstat,
                    file.path(mrdir, str_c(name_region, pheno,
                                           "MR_Fstat_results.csv", sep="_")),
                    col.names=FALSE, row.names=FALSE, sep=',', quote=FALSE)
    }, pheno=pheno)
}

MRplot <- function(dat, mr_results) {
    p_mr <- mr_scatter_plot(mr_results, dat)
    res_single <- mr_singlesnp(dat,
                               all_method=c("mr_egger_regression", "mr_ivw",
                                            "mr_weighted_median",
                                            "mr_weighted_mode"))
    p_forest <- mr_forest_plot(res_single)

    res_loo <- mr_leaveoneout(dat)
    p_loo <- mr_leaveoneout_plot(res_loo)

    p_funnel <- mr_funnel_plot(res_single)
    legendMR <- cowplot::get_legend(p_funnel[[1]] + theme_bw() +
                                        theme(legend.position='left') )

    all_plots <- lapply(seq_along(p_forest), function(x) {
        mr <- p_mr[[x]] + theme_bw()  + theme(legend.position='none')
        forest <- p_forest[[x]] + theme_bw() + theme(legend.position='none')
        loo <- p_loo[[x]] + theme_bw() + theme(legend.position='none')
        funnel <-  p_funnel[[x]] + theme_bw() + theme(legend.position='none')
        plots_1 <- cowplot::plot_grid(mr, funnel, nrow=2, align='v',
                                      axis="lr")
        plots_2 <- cowplot::plot_grid(forest,  loo, nrow=2, align='v',
                                      axis="lr")
        plots <- cowplot::plot_grid(plots_1, plots_2, ncol=2, align='h',
                                    axis="tb")
        p_all <- cowplot::plot_grid(legendMR, plots, ncol=2, rel_widths=c(1,4))
    })
    mr_panel_plots <- lapply(seq_along(p_forest), function(x) {
        mr <- p_mr[[x]] + theme_bw()  + theme(legend.position='none')
        forest <- p_forest[[x]] + theme_bw() + theme(legend.position='none')
        cowplot::plot_grid(mr, forest, ncol=2, align='h', axis="tb")
    })
    return(list(all_plots=all_plots, forest=p_forest, mr_panel=mr_panel_plots))
}


MR <- function(exposure_list=NULL, outcome_list=NULL, exposure=NULL,
               outcome=NULL, name, mrdir) {
    if (!is.null(exposure_list) & !is.null(exposure)) {
        stop("Only one type of exposure can be provided: specify either",
             "exposure_list or exposure")
    }
    if (!is.null(outcome_list) & !is.null(outcome)) {
        stop("Only one type of outcome can be provided: specify either",
             "outcome_list or outcome")
    }
    if (!is.null(exposure_list)) {
        if (!is.null(outcome_list)) {
            MRbase <- list(basal=getMRdata(data_exposure = exposure_list[[1]],
                                           data_outcome = outcome_list[[1]]),
                           mid=getMRdata(data_exposure = exposure_list[[2]],
                                         data_outcome = outcome_list[[2]]),
                           apical=getMRdata(data_exposure = exposure_list[[3]],
                                            data_outcome = outcome_list[[3]]))
        } else {
            MRbase <- list(basal=getMRdata(data_exposure = exposure_list[[1]],
                                           data_outcome = outcome),
                           mid=getMRdata(data_exposure = exposure_list[[2]],
                                         data_outcome = outcome),
                           apical=getMRdata(data_exposure = exposure_list[[3]],
                                            data_outcome = outcome))
        }
    } else {
        if (!is.null(outcome_list)) {
            MRbase <- list(basal=getMRdata(data_exposure = exposure,
                                           data_outcome = outcome_list[[1]]),
                           mid=getMRdata(data_exposure = exposure,
                                         data_outcome = outcome_list[[2]]),
                           apical=getMRdata(data_exposure = exposure,
                                            data_outcome = outcome_list[[3]]))
        } else {
            MRbase <- list(mr=getMRdata(data_exposure = exposure,
                                           data_outcome = outcome))
        }
    }
    saveRDS(MRbase, file.path(mrdir,
                              str_c("MRbase", name, ".rds", sep="_")))

    MRresults <- lapply(MRbase, function(x) {
        MRanalysis(x$exposure, x$outcome)
    })

    mr_tables <-  MRtables(MRresults, name)

    mr_plots <- lapply(MRresults, function(region) {
        dat <- region$dat
        res <- region$mr_results
        tmp <- MRplot(dat, res)
    })

    panel_plots <- lapply(mr_plots, function(x) x$all_plots)
    panels <- lapply(panel_plots, function(x) x[[1]])
    p_panels <- cowplot::plot_grid(plotlist=panels, nrow=3)
    ggsave(plot=p_panels,
           file.path(mrdir, str_c("MR_panels", name, ".pdf", sep="_")),
           height=20, width=12)

    forests <- lapply(mr_plots, function(x) x$forest[[1]] + theme_bw() +
                          theme(legend.position='none',
                                axis.title.x = element_blank()))
    p_forests <- cowplot::plot_grid(plotlist=forests, ncol=3)
    ggsave(plot=p_forests,
           file.path(mrdir, str_c("MR_forest", name, ".pdf", sep="_")),
           height=4, width=12)

    mr_forest <- lapply(mr_plots, function(x) x$mr_panel)
    return(list(plot=mr_plots, forest=p_forests, mr_forest=mr_forest))
}

############
## data  ###
############

## command line arguments ####
option_list <- list(
    optparse$make_option(c("--gwasdir"), action="store",
               dest="gwasdir",
               type="character", help="Path the ukbb gwas directory
                [default: %default].", default=NULL),
    optparse$make_option(c("--mrdir"), action="store",
               dest="mrdir",
               type="character", help="Path to output MR directory
                [default: %default].", default=NULL),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out ",
               "[default: %default]."),
    optparse$make_option(c("--debug"), action="store_true",
                dest="debug", default=FALSE, type="logical",
                help="If set, predefined arguments are used to test the script",
                "[default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$gwasdir <- "~/data/ukbb/ukb-hrt/gwas/180628_fractal_dimension"
    args$mrdir <- "~/data/ukbb/ukb-hrt/MR/180628_fractal_dimension"
    args$verbose <- TRUE
}
directory <- args$directory
gwasdir <- args$gwasdir
mrdir <- args$mrdir
verbose <- args$verbose

## ld-filtered, significant genome-wide association results ukb ####
slices_sig <- read.table(paste(gwasdir,
                               "/Pseudomultitrait_Slices_sig5e08_ldFiltered.txt",
                               sep=""),
                         sep=",", stringsAsFactors=FALSE, header=TRUE)
# LD filter misses these two SNPs, manually remove
slices_sig <- slices_sig[!slices_sig$rsid %in% c("rs12214483", "rs117953218"),]

slices_sig$SNPID <- paste(slices_sig$chr, ":", slices_sig$pos, "_",
                          slices_sig$a_0, "_", slices_sig$a_1, sep="")

## genotypes of ld-filtered, significant genome-wide association results ####
geno_sig <- fread(file.path(gwasdir,
                            "Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz"),
                  stringsAsFactors=FALSE, data.table=FALSE)
geno_sig$SNPID <- paste(geno_sig$chromosome, ":", geno_sig$position, "_",
                        geno_sig$alleleA, "_", geno_sig$alleleB, sep="")
geno_sig <- geno_sig[geno_sig$SNPID %in% slices_sig$SNPID,]
geno_sig <- merge(slices_sig[, c(2,6)], geno_sig, by='rsid')

###############
## analysis ###
###############

## format slice betas and p-values per slice ####
slices_beta <- cbind(rsid=slices_sig$rsid,
                     a_0=slices_sig$a_0,
                     a_1=slices_sig$a_1, af=slices_sig$af, samplesize=18097,
                     slices_sig[,grepl('beta', colnames(slices_sig))])
colnames(slices_beta) <- gsub("_beta", "", colnames(slices_beta))
beta <- reshape2::melt(slices_beta, id.var=c('rsid',
                                             'a_0', 'a_1', 'af',
                                             'samplesize'),
                                             value.name='beta',
                       variable.name='slice')

slices_se <- cbind(rsid=slices_sig$rsid, a_0=slices_sig$a_0,
                   a_1=slices_sig$a_1,
                   slices_sig[,grepl('se', colnames(slices_sig))])
colnames(slices_se) <- gsub("_se", "", colnames(slices_se))
se <- reshape2::melt(slices_se, id.var=c('rsid', 'a_0', 'a_1'),
                     value.name='se',
                     variable.name='slice')

slices_logp <- cbind(rsid=slices_sig$rsid, a_0=slices_sig$a_0,
                     a_1=slices_sig$a_1,
                     slices_sig[, grepl('log10p', colnames(slices_sig))])
colnames(slices_logp) <- gsub("\\.log10p", "", colnames(slices_logp))
logp <- reshape2::melt(slices_logp, id.var=c('rsid', 'a_0', 'a_1'),
                       value.name='logp',
                       variable.name='slice')

slices <- cbind(beta, p=10^(-logp$logp), se=se$se)
slices$rsid <- as.character(slices$rsid)
slices$a_0 <- as.character(slices$a_0)
slices$a_1 <- as.character(slices$a_1)

## slices significant slice pvalues and betas ####
sig_per_slice <- dplyr::filter(slices, p  < 5e-8)
write.table(sig_per_slice,
            file.path(gwasdir, 'Significant_per_slice.csv'),
            sep=',', col.names=TRUE, row.names=FALSE)
# sig_per_slice <- read.table(file.path(gwasdir, 'Significant_per_slice.csv'),
# sep=',', header=TRUE, stringsAsFactors = FALSE)

# analyse betas slices ####
# from https://jmarchini.org/bgenie-usage/: beta coefficient refers to the
# effect of having an extra copy of a_1
colnames(sig_per_slice)[1:9] <- c('rsid', 'a_0', 'a_1', 'af', 'samplesize',
                             'slices', 'slices_beta', 'slices_p', 'slices_se')

base <- as.character(sig_per_slice$slices) %in%
    c('Slice_1', 'Slice_2', 'Slice_3')
mid <- as.character(sig_per_slice$slices) %in%
    c('Slice_4', 'Slice_5', 'Slice_6')
sig_per_slice$area <- "apical"
sig_per_slice$area[base] <- 'basal'
sig_per_slice$area[mid] <- 'mid'
sig_per_slice$area <- factor(sig_per_slice $area,
                              levels=c('basal', 'mid', 'apical'))

sig_per_area <- split(sig_per_slice , f=sig_per_slice $area)
names(sig_per_area) <- levels(sig_per_slice$area)

unique_per_area <- lapply(seq_along(sig_per_area), function(test_area) {
    area <- sig_per_area[[test_area]]
    colnames(area)[1:9] <- c('SNP', 'other_allele', 'effect_allele', 'eaf',
                             'samplesize', 'slices', 'beta', 'pval', 'se')
    if (any(duplicated(area$SNP))) {
        dup <- area[area$SNP %in% area$SNP[duplicated(area$SNP)],]
        remove <- sapply(unique(dup$SNP), function(x) {
            minbeta <- sort(abs(area$beta[area$SNP %in% x]),decreasing=TRUE)[-1]
            which(abs(area$beta) %in% minbeta)
        })
        area <- area[-unlist(remove),]
    }
    area$id <- 'FD'
    area$Phenotype <- 'FD'
    return(area)
})
names(unique_per_area) <- names(sig_per_area)

exposure_per_area <- lapply(seq_along(unique_per_area), function(test_area) {
    area <- unique_per_area[[test_area]]
    area_exposure <- format_data(area, type='exposure')
    area_exposure$id.exposure <- 'FD'
    area_exposure$exposure <- 'FD'
    write.table(area_exposure, paste(mrdir, '/',
                                  names(unique_per_area)[test_area],
                            "_association_results.txt", sep=""),
                sep=' ', col.names=TRUE, row.names=FALSE, quote=FALSE)
    return(area_exposure)
})
names(exposure_per_area) <- names(unique_per_area)

## MRbase MR using BOLT LMM output files provided by Gib Hermani, 2/19/20 ####
# A1: effect allele
# UKB-b:6025 - n = 4525, UKB-b:2240 - n = 7796 (Gib Hermani, 2/24/20)

## Stroke volume
sv <- read_delim("~/data/ukbb/ukb-hrt/MR/UKB-b6025.txt", col_names = TRUE,
                     delim="\t")
sv <- sv %>%
    select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, P_LINREG, SE) %>%
    rename(other_allele=ALLELE0, effect_allele=ALLELE1, eaf=A1FREQ, beta=BETA,
           pval=P_LINREG, se=SE) %>%
    distinct() %>%
    mutate(Phenotype="SV") %>%
    mutate(id="UKB-b:6025") %>%
    mutate(samplesize=4525) %>%
    filter(SNP %in% sig_per_slice$rsid)

write_csv(sv, "~/data/ukbb/ukb-hrt/MR/UKB-b6025_SV_outcome.csv")

mr_sv <- MR(exposure_list=unique_per_area, outcome=sv, name="SV",
            mrdir=mrdir)

## QRS duration
qrs <- read_delim("~/data/ukbb/ukb-hrt/MR/UKB-b2240.txt", col_names = TRUE,
                 delim="\t")
qrs <- qrs %>%
    select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, P_LINREG, SE) %>%
    rename(other_allele=ALLELE0, effect_allele=ALLELE1, eaf=A1FREQ, beta=BETA,
           pval=P_LINREG, se=SE) %>%
    distinct() %>%
    mutate(Phenotype="QRS") %>%
    mutate(id="UKB-b:2240")  %>%
    mutate(samplesize=7796) %>%
    filter(SNP %in% sig_per_slice$rsid)

write_csv(qrs, "~/data/ukbb/ukb-hrt/MR/UKB-b2240_QRS_outcome.csv")

mr_qrs<- MR(exposure_list=unique_per_area, outcome=qrs, name="QRS",
            mrdir=mrdir)

# HERMES MR
# A1: effect allele, beta:log odds ratio of heart failure per extra copy of A1
hermes <- read_delim("~/data/ukbb/HERMES/lookup_results.tsv", col_names = TRUE,
                     delim="\t")

hermes <- hermes %>%
    select(SNP, A2, A1, freq, N, b, p, se) %>%
    rename(other_allele=A2, effect_allele=A1, eaf=freq, samplesize=N,
           beta=b, pval=p) %>%
    distinct() %>%
    mutate(Phenotype="HF") %>%
    mutate(id="HERMES")

mr_hermes <- MR(exposure_list=unique_per_area, outcome=hermes,
                name="HF", mrdir=mrdir)

all_formated <- do.call(rbind, unique_per_area)
dup <- all_formated[all_formated$SNP %in% all_formated$SNP[duplicated(all_formated$SNP)],]
remove <- sapply(unique(dup$SNP), function(x) {
    minbeta <- sort(abs(all_formated$beta[all_formated$SNP %in% x]),decreasing=TRUE)[-1]
    which(abs(all_formated$beta) %in% minbeta)
})
all_formated <- all_formated[-unlist(remove),]

mr_hermes_all <- MR(exposure=all_formated, outcome=hermes,
                name="HF_all", mrdir=mrdir)

## HERMES <-> SV,QRS
exposure_sv <- list(basal=sv[sv$SNP %in% unique_per_area$basal$SNP,],
                    mid=sv[sv$SNP %in% unique_per_area$mid$SNP,],
                    apical=sv[sv$SNP %in% unique_per_area$apical$SNP,])
exposure_qrs <- list(basal=qrs[qrs$SNP %in% unique_per_area$basal$SNP,],
                    mid=qrs[qrs$SNP %in% unique_per_area$mid$SNP,],
                    apical=qrs[qrs$SNP %in% unique_per_area$apical$SNP,])
outcome_hermes <- list(basal=hermes[hermes$SNP %in% unique_per_area$basal$SNP,],
                    mid=hermes[hermes$SNP %in% unique_per_area$mid$SNP,],
                    apical=hermes[hermes$SNP %in% unique_per_area$apical$SNP,])

mr_hermes_sv <- MR(exposure_list=exposure_sv, outcome_list=outcome_hermes,
                   name="SV_HF", mrdir=mrdir)

mr_hermes_sv_all <- MR(exposure=sv, outcome=hermes, name="SV_HF_all",
                       mrdir=mrdir)

mr_hermes_qrs <- MR(exposure_list=exposure_qrs, outcome_list=outcome_hermes,
                   name="QRS_HF", mrdir=mrdir)

mr_hermes_qrs_all <- MR(exposure=qrs, outcome=hermes, name="QRS_HF_all",
                       mrdir=mrdir)

## DCM MR
# A1: effect allele, beta:log odds ratio of heart failure per extra copy of A1
dcm <- read_delim("~/data/digital-heart/gwas/DCM/HVOL.DCM.FD.Pseudomultitrait_Slices_sig5e08.sigSNPs.clean.assoc.logistic.all", col_names=TRUE, delim = " ")

dcm <- dcm %>%
    select(SNP, A1, A2, MAF, NMISS, BETA, SE, ` EMP1`) %>%
    rename(effect_allele=A1, other_allele=A2, eaf=MAF, samplesize=NMISS,
           beta=BETA, pval=` EMP1`, se=SE) %>%
    mutate(Phenotype="DCM") %>%
    mutate(id="DCM")

dcm <- sig_per_slice %>%
    select(rsid, a_0, a_1) %>%
    inner_join(dcm, by=c("rsid"="SNP")) %>%
    mutate(other_allele=case_when(effect_allele == a_0 ~ a_1,
                                  TRUE ~ a_0)) %>%
    select(rsid, other_allele, effect_allele, eaf, beta, se, pval, samplesize,
           Phenotype, id) %>%
    distinct() %>%
    rename(SNP="rsid")


mr_dcm <- MR(exposure_list=unique_per_area, outcome=dcm, name="DCM",
             mrdir=mrdir)

mr_dcm_all <- MR(exposure=all_formated, outcome=dcm,
                    name="DCM_all", mrdir=mrdir)

## Combine mr_forest plots for HERMES and DCM

hermes_dcm <- cowplot::plot_grid(mr_hermes$mr_forest[[1]][[1]],
                                 mr_dcm$mr_forest[[1]][[1]],
                   ncol=1)
ggsave(plot=hermes_dcm,
       file.path(mrdir, "MR_hermes_dcm.pdf"),
       height=8, width=12)

hermes_dcm_all <- cowplot::plot_grid(mr_hermes_all$mr_forest[[1]][[1]],
                                 mr_dcm_all$mr_forest[[1]][[1]],
                                 ncol=1)
ggsave(plot=hermes_dcm_all,
       file.path(mrdir, "MR_hermes_dcm_all.pdf"),
       height=8, width=12)





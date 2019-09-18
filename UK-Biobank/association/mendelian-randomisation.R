###############################
### Libraries and functions ###
###############################
options(bitmapType = 'cairo', device = 'pdf')

modules::import_package('dplyr', attach=TRUE)
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
    }
    if (!is.null(outcomes)) {
        ao <- available_outcomes(access)
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
    per_study_F$samplesize.outcome <- dat$samplesize.outcome[!duplicated(dat$outcome)]
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
           directionality_results=directionality_results, I2=I2,
           Fstat=per_study_F))
}

plotMR <- function(dat, mr_results) {
    p_mr <- mr_scatter_plot(mr_results, dat)
    res_single <- mr_singlesnp(dat,
                               all_method=c("mr_egger_regression", "mr_ivw",
                                            "mr_weighted_median",
                                            "mr_weighted_mode"))
    p_forrest <- mr_forest_plot(res_single)
    
    res_loo <- mr_leaveoneout(dat)
    p_loo <- mr_leaveoneout_plot(res_loo)
    
    p_funnel <- mr_funnel_plot(res_single)
    legendMR <- cowplot::get_legend(p_funnel[[1]] + theme_bw() +
                                        theme(legend.position='left') )
    
    all_plots <- lapply(seq_along(p_forrest), function(x) {
        mr <- p_mr[[x]] + theme_bw()  + theme(legend.position='none')
        forrest <- p_forrest[[x]] + theme_bw() + theme(legend.position='none')
        loo <- p_loo[[x]] + theme_bw() + theme(legend.position='none')
        funnel <-  p_funnel[[x]] + theme_bw() + theme(legend.position='none')
        plots_1 <- cowplot::plot_grid(mr, funnel, nrow=2, align='v')
        plots_2 <- cowplot::plot_grid(forrest,  loo, nrow=2, align='v')
        plots <- cowplot::plot_grid(plots_1, plots_2, ncol=2, align='h')
        p_all <- cowplot::plot_grid(legendMR, plots, ncol=2, rel_widths=c(1,4))
    })
    return(list(all_plots=all_plots, forrest=p_forrest))
}

############
## data  ###
############

## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path ukbb root data rootdirectory
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
    args$directory <- "~/data/ukbb/ukb-hrt"
    args$verbose <- TRUE
    args$Teff <- 6.6
}
directory <- args$directory
Teff <- args$Teff
verbose <- args$verbose

## ld-filtered, significant genome-wide association results ukb ####
slices_sig <- read.table(paste(directory,
                               "/gwas/Pseudomultitrait_Slices_sig5e08_ldFiltered.txt",
                               sep=""),
                         sep=",", stringsAsFactors=FALSE, header=TRUE)
# LD filter misses these two SNPs, manually remove
slices_sig <- slices_sig[!slices_sig$rsid %in% c("rs12214483", "rs117953218"),]

slices_sig$SNPID <- paste(slices_sig$chr, ":", slices_sig$pos, "_",
                          slices_sig$a_0, "_", slices_sig$a_1, sep="")

## genotypes of ld-filtered, significant genome-wide association results ####
cmd=paste("zcat ", directory,
          "/gwas/Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz", sep="")
geno_sig <- data.table::fread(cmd=cmd, stringsAsFactors=FALSE, data.table=FALSE)
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
            paste(directory, '/gwas/Significant_per_slice.csv', sep=""),
            sep=',', col.names=TRUE, row.names=FALSE)

# analyse betas slices ####
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
    colnames(area)[1:9] <- c('SNP', 'effect_allele', 'other_allele', 'eaf',
                             'samplesize', 'slices', 'beta', 'pval', 'se')
    if (any(duplicated(area$SNP))) {
        dup <- area[area$SNP %in% area$SNP[duplicated(area$SNP)],]
        remove <- sapply(unique(dup$SNP), function(x) {
            minbeta <- sort(abs(area$beta[area$SNP %in% x]),decreasing=TRUE)[-1]
            which(abs(area$beta) %in% minbeta)
        })
        area <- area[-unlist(remove),]
    }
    return(area)
})
names(unique_per_area) <- names(sig_per_area)

exposure_per_area <- lapply(seq_along(unique_per_area), function(test_area) {
    area <- unique_per_area[[test_area]]
    area_exposure <- format_data(area, type='exposure')
    area_exposure$id.exposure <- 'FD'
    area_exposure$exposure <- 'FD'
    write.table(area_exposure, paste(directory, '/',
                                  names(unique_per_area)[test_area],
                            "_association_results.txt", sep=""),
                sep=' ', col.names=TRUE, row.names=FALSE, quote=FALSE)
    return(area_exposure)
})
names(exposure_per_area) <- names(unique_per_area)


## MR base: `two-sample` MR with Biobank data ####
token <- googleAuthR::gar_auth(paste(directory, "/MR/mrbase.oauth", sep=''))
access <- token$credentials$access_token

#SV, QRS duration, HR
outcomes <- c('UKB-b:6025', 'UKB-b:2240', '1056')

basalMRbase <- getMRdata(data_exposure=exposure_per_area[[1]],
                       outcomes=outcomes, access=access)
midMRbase <- getMRdata(data_exposure=exposure_per_area[[2]],
                       outcomes=outcomes, access=access)
apicalMRbase <- getMRdata(data_exposure=exposure_per_area[[3]],
                          outcomes=outcomes, access=access)

MRbase <- list(basal=basalMRbase, mid=midMRbase, apical=apicalMRbase)
saveRDS(MRbase, paste(directory, "/MR/MRbase.rds", sep=""))

MRbase_results <- lapply(MRbase, function(x) MRanalysis(x$exposure, x$outcome))

mr_tables <- lapply(seq_along(MRbase_results), function(x) {
    region <- MRbase_results[[x]]
    name_region <- names(MRbase_results)[x]
    write.table(region$mr_results,  paste(directory, "/MR/", name_region,
                                          "_MR_results.csv", sep=""),
                col.names=TRUE, row.names=FALSE, sep=',', quote=FALSE)
    write.table(region$plei_results,  paste(directory, "/MR/", name_region,
                                          "_MR_pleiotropy_results.csv", sep=""),
                col.names=TRUE, row.names=FALSE, sep=',', quote=FALSE)
    write.table(region$directionality_results,
                paste(directory, "/MR/", name_region,
                      "_MR_directionality_results.csv", sep=""),
                col.names=TRUE, row.names=FALSE, sep=',', quote=FALSE)

    write.table(region$I2, paste(directory, "/MR/", name_region,
                                            "_MR_I2_results.csv", sep=""),
                col.names=FALSE, row.names=FALSE, sep=',', quote=FALSE)
    write.table(region$Fstat, paste(directory, "/MR/", name_region,
                                 "_MR_Fstat_results.csv", sep=""),
                col.names=FALSE, row.names=FALSE, sep=',', quote=FALSE)
})

mr_plots <- lapply(MRbase_results, function(region) {
    dat <- region$dat
    res <- region$mr_results
    tmp <- plotMR(dat, res)
})

panels_hr <- lapply(panel_plots, function(x) x[[1]])
p_panels_hr <- cowplot::plot_grid(plotlist=panels_hr, nrow=3)
ggsave(plot=p_panels_hr, paste(directory, "/MR/MR_panels_HR.pdf", sep=""),
       height=20, width=12)

panels_qrs <- lapply(panel_plots, function(x) x[[4]])
p_panels_qrs <- cowplot::plot_grid(plotlist=panels_qrs, nrow=3)
ggsave(plot=p_panels_qrs, paste(directory, "/MR/MR_panels_QRS.pdf", sep=""),
       height=20, width=12)

panel_plots <- lapply(mr_plots, function(x) x$all_plots)
panels_sv <- lapply(panel_plots, function(x) x[[5]])
p_panels_sv <- cowplot::plot_grid(plotlist=panels_sv, nrow=3)
ggsave(plot=p_panels_sv, paste(directory, "/MR/MR_panels_SV.pdf", sep=""),
       height=20, width=12)


forrests_sv <- lapply(mr_plots, function(x) x$forrest[[3]] + theme_bw() +
                          scale_x_continuous(limits=c(-1,4)) +
                          theme(legend.position='none',
                              axis.title.x = element_blank()))
p_forrests_sv <- cowplot::plot_grid(plotlist=forrests_sv, ncol=3)
ggsave(plot=p_forrests_sv, paste(directory, "/MR/MR_forrest_SV.pdf", sep=""),
       height=4, width=12)

forrests_qrs <- lapply(mr_plots, function(x) x$forrest[[2]] + theme_bw() +
                          scale_x_continuous(limits=c(-1,4)) +
                          theme(legend.position='none',
                                axis.title.x = element_blank()))
p_forrests_qrs <- cowplot::plot_grid(plotlist=forrests_qrs, ncol=3)
ggsave(plot=p_forrests_qrs, paste(directory, "/MR/MR_forrest_QRS.pdf", sep=""),
       height=4, width=12)


###############################
### Libraries and functions ###
###############################
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('dplyr', attach=TRUE)
optparse <- modules::import_package('optparse')
rmr <- modules::import_package('RadialMR')

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

alleleScore <- function(bx, by, bxse, byse) {
    # for equal weights, otherwise specify weights
    wts <- rep(1, length(bx))
    # equivalent to original IVW method when wts = bx
    # standard error from delta method (first-order approximation)
    se_SSw_first <- sqrt(sum(wts^2/byse^2)/sum(wts*bx/byse^2))
    # standard error from delta method (second-order approximation)
    # theta is the correlation between the numerator and denominator of the estimate
    # if the correlation is not known, it can be taken as the observational
    # correlation between the risk factor and outcome; # a sensitivity analysis can also be performed for its value
    se_SSw_second = sqrt(sum(wts^2/byse^2)/sum(wts*bx/byse^2)^2 + sum(wts*by/byse^2)^2/sum(wts*bx/byse^2)^4*sum(wts^2/byse^2) - 2*theta*sum(wts*by/byse^2)/sum(wts*bx/byse^2)^3) 
    beta_SSw = sum(wts*by/byse^2)/sum(wts*bx/byse^2) 
If the genetic variants are correlated

getMRdata <- function(filename, outcomes, access) {
    ao <- available_outcomes(access)
    exposure_dat <- read_exposure_data(
        filename = filename
    )
    outcome_dat <- extract_outcome_data(exposure_dat$SNP, outcomes,
                                        proxies=TRUE, rsq=0.8, align_alleles=1,
                                        palindromes=1, maf_threshold=0.3,
                                        access_token=access)
    return(list(exposure=exposure_dat, outcome=outcome_dat))
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

    if (verbose) message("MR I^2 analysis")
    weighted_beta <- dat$beta.exposure/dat$se.outcome
    weighted_se <- dat$se.exposure/dat$se.outcome
    I2 <- estimateI2(weighted_beta,weighted_se)
    return(list(dat=dat, mr_results=mr_results, het_results=het_results,
           plei_results=plei_results, loo_results=loo_results,
           directionality_results=directionality_results, I2=I2))
}
############
## data  ###
############

## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path to directory with digital-heart
               bgenie association results [default: %default].", default=NULL),
    optparse$make_option(c("-ukb", "--ukbdir"), action="store",
               dest="ukbdir",
               type="character", help="Path to directory with ukbb significant
               association results [default: %default].", default=NULL),
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
    args$directory <- "~/data/ukbb/ukb-hrt/gwas"
    args$name <- 'volumes'
    args$verbose <- TRUE
}
directory <- args$directory
name <- args$name
verbose <- args$verbose

## genome-wide association results volumes ####
lvv <- data.table::fread(paste(directory, "/bgenie_", name, "_lm_st_genomewide.csv",
                               sep=""), data.table=FALSE, stringsAsFactors=FALSE)

## ld-filtered, significant genome-wide association results ukb ####
slices_sig <- read.table(paste(directory, "/Slices_sig5e08_ldFiltered.txt", sep=""),
                         sep=",", stringsAsFactors=FALSE, header=TRUE)
slices_sig$SNPID <- paste(slices_sig$chr, ":", slices_sig$pos, "_",
                          slices_sig$a_0, "_", slices_sig$a_1, sep="")

## genotypes of ld-filtered, significant genome-wide association results ####
cmd=paste("zcat ", directory,
          "/Pseudomultitrait_slices_sig5e08_genotypes.dosage.gz", sep="")
geno_sig <- data.table::fread(cmd=cmd, stringsAsFactors=FALSE, data.table=FALSE)
geno_sig$SNPID <- paste(geno_sig$chromosome, ":", geno_sig$position, "_",
                        geno_sig$alleleA, "_", geno_sig$alleleB, sep="")
geno_sig <- geno_sig[geno_sig$SNPID %in% slices_sig$SNPID,]
geno_sig <- merge(slices_sig[, c(2,6)], geno_sig, by='rsid')
geno_recode <- apply(geno_sig[, c(2,8:ncol(geno_sig))], 1, function(x) {
    if(x[1] > 0.5) return(abs(x[-1]-2))
    return(x[-1])
})
geno_recode <- cbind(geno_sig[,1:7], t(geno_recode))
colnames(geno_recode) <- colnames(geno_sig)


## phenotyes and covariates ####
slices <- data.table::fread('~/data/ukbb/ukb-hrt/phenotypes/FD_slices_EUnorel.csv',
                           stringsAsFactors=FALSE, data.table=FALSE)
summary <- data.table::fread('~/data/ukbb/ukb-hrt/phenotypes/FD_summary_EUnorel.csv',
                           stringsAsFactors=FALSE, data.table=FALSE)
covs <- data.table::fread('~/data/ukbb/ukb-hrt/phenotypes/FD_covariates_EUnorel.csv',
                          stringsAsFactors=FALSE, data.table=FALSE)

###############
## analysis ###
###############



## format slice betas and p-values per slice ####
slices_beta <- cbind(rsid=slices_sig$rsid, a_0=slices_sig$a_0,
                     a_1=slices_sig$a_1, af=slices_sig$af, samplesize=18097,
                     slices_sig[,grepl('beta', colnames(slices_sig))])
colnames(slices_beta) <- gsub("_beta", "", colnames(slices_beta))
beta <- reshape2::melt(slices_beta, id.var=c('rsid', 'a_0', 'a_1', 'af',
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
write.table(sig_per_slice, paste(directory, '/Significant_per_slice.csv',
                                 sep=''),
            sep=",", col.names=TRUE, row.names=FALSE)

sig_per_slice <- sig_per_slice[sig_per_slice$rsid %in% lvv$rsid, ]

## betas and pvalues of significant snps and slice in ukb ####
lvv_beta <- lvv[,c(2, which(grepl('beta', colnames(lvv))))]
colnames(lvv_beta) <- gsub("_beta", "", colnames(lvv_beta))
lvv_logp <- lvv[,c(2, which(grepl('log10p', colnames(lvv))))]
colnames(lvv_logp) <- gsub("-log10p", "", colnames(lvv_logp))
lvv_se <- lvv[,c(2, which(grepl('se', colnames(lvv))))]
colnames(lvv_se) <- gsub("_se", "", colnames(lvv_se))

lvv_sig <- do.call(rbind, apply(sig_per_slice, 1, function(x) {
    pos <- which(lvv_beta$rsid %in%  x[[1]])
    beta_tmp <- lvv_beta[pos, -1]
    names(beta_tmp) <- paste(names(beta_tmp), "_beta", sep="")
    p_tmp <- 10^(-lvv_logp[pos, -1])
    colnames(p_tmp) <- paste(colnames(p_tmp), "_p", sep="")
    se_tmp <- lvv_se[pos, -1]
    names(se_tmp) <- paste(names(se_tmp), "_se", sep="")
    tmp <- data.frame(beta_tmp, p_tmp, se_tmp, stringsAsFactors=FALSE)
    return(tmp)
}))

write.table(lvv_sig, paste(directory, '/lvv_sig.csv',
                                 sep=''),
            sep=",", col.names=TRUE, row.names=FALSE)

## analyse betas slices and volumes ####
combined <- cbind(sig_per_slice, lvv_sig)
colnames(combined)[1:9] <- c('rsid', 'a_0', 'a_1', 'af', 'samplesize',
                             'slices', 'slices_beta', 'slices_p', 'slices_se')

base <- as.character(combined$slices) %in% c('Slice_1', 'Slice_2', 'Slice_3')
mid <- as.character(combined$slices) %in% c('Slice_4', 'Slice_5', 'Slice_6')
combined$area <- "apical"
combined$area[base] <- 'basal'
combined$area[mid] <- 'mid'
combined$area <- factor(combined$area, levels=c('basal', 'mid', 'apical'))

combined_per_area <- split(combined, f=combined$area)
names(combined_per_area) <- levels(combined$area)

unique_per_area <- lapply(seq_along(combined_per_area), function(test_area) {
    area <- combined_per_area[[test_area]]
    colnames(area)[1:9] <- c('SNP', 'effect_allele', 'other_allele', 'eaf',
                             'samplesize', 'slices', 'beta', 'pval', 'se')
    if (any(duplicated(area$SNP))) {
        dup <- area[area$SNP %in% area$SNP[duplicated(area$SNP)],]
        remove <- sapply(unique(dup$SNP), function(x) {
                   minbeta <- sort(abs(area$beta[area$SNP %in% x]),
                                   decreasing=TRUE)[-1]
                   which(abs(area$beta) %in% minbeta)
        })
        area <- area[-unlist(remove),]
    }
    write.table(area[,1:9], paste(directory, '/',
                                  names(combined_per_area)[test_area],
                            "_association_results.txt", sep=""),
                sep=' ', col.names=TRUE, row.names=FALSE, quote=FALSE)
    return(area)
})
names(unique_per_area) <- names(combined_per_area)

## estimate % of variance of sig loci ####
genotypes <- t(geno_sig[, -c(1:7)])
colnames(genotypes) <- geno_sig$rsid

ve_sig <- sapply(seq_along(unique_per_area), function(test_area) {
    area <- unique_per_area[[test_area]]
    pheno <- summary[, c(1, 1+test_area)]
    colnames(pheno) <- c("IID", "pheno")
    #gt <- rowSums(genotypes[, colnames(genotypes) %in% area$SNP])
    gt <- genotypes[, colnames(genotypes) %in% area$SNP]

    df_null <- merge(pheno, covs, by=1)
    df_alt <- merge(df_null, gt, by.x=1, by.y=0)

    lm_null <- lm(pheno ~ ., df_null[, -1])
    lm_alt <- lm(pheno ~ ., df_alt[, -1])
    r2 <- summary(lm_alt)$adj.r.squared - summary(lm_null)$adj.r.squared
    N <- dim(geno_sig)[2]
    k <- dim(gt)[2]
    F <- (N - k -1)/k * r2/(1-r2)
    low <- findlowerflimit(F, k, N)
    return(c(F, low))
})

findlowerflimit <- function(f, nu1, nu2) {
    lambda <- f*nu1*(nu2-2)/nu2-nu1
    lower <- f - 1
    while (pf(lower, df1=nu1, df2=nu2, ncp=lambda) > 0.05) {
        lower <- lower-1
    }
    upper <- lower + 1
    while (abs(pf((lower+upper)/2, df1=nu1, df2=nu2, ncp=lambda)-0.05) > 0.0001 ) {
        if (pf((lower+upper)/2, df1=nu1, df2=nu2, ncp=lambda) > 0.05) {
            upper = (lower+upper)/2
        }
        if (pf((lower+upper)/2, df1=nu1, df2=nu2, ncp=lambda) <0.05) {
        lower = (lower+upper)/2
        }
    }
    return((lower+upper)/2)
}

## MR base ####
library(TwoSampleMR)
token <- googleAuthR::gar_auth( paste(directory, "/mrbase.oauth", sep=''))
access <- token$credentials$access_token


filename <- paste(directory, '/mid_association_results.txt', sep="")
outcomes <- c('UKB-b:6025','1056')

basalMRdata <- getMRdata(paste(directory, '/basal_association_results.txt', sep=""),
                       outcomes, access)
midMRdata <- getMRdata(paste(directory, '/mid_association_results.txt', sep=""),
                       outcomes, access)
apicalMRdata <- getMRdata(paste(directory, '/apical_association_results.txt', sep=""),
                       outcomes, access)

MRdata <- list(basal=basalMRdata, mid=midMRdata, qpical=apicalMRdata)
saveRDS(MRdata, paste(directory, "/MR/MRdata.rds", sep=""))
results <- lapply(MRdata, function(x) MRanalysis(x$exposure, x$outcome))
saveRDS(results, paste(directory, "/MR/MRresults.rds", sep=""))


p1 <- mr_scatter_plot(res, dat)

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]




mr_area <- lapply(seq_along(combined_per_area), function(test_area) {
    area <- combined_per_area[[test_area]]
    tmp <- lapply(c('CO', 'SV', 'HR'), function(x){
        if (x == 'CO') {
            BYG=area$CO_beta
            seBYG=area$CO_se
        }
        if (x == 'SV') {
            BYG=area$SV_beta
            seBYG=area$SV_se
        }
        if (x == 'HR') {
            BYG=area$HR_beta
            seBYG=area$HR_se
        }
        formated <- rmr$format_radial(BXG=area$slices_beta, BYG=BYG,
                                  seBXG=area$slices_se,
                                  seBYG=seBYG, RSID=area$rsid)
        ivw <- rmr$ivw_radial(r_input=formated, alpha=0.05, weights=1)
        egger <- rmr$egger_radial(r_input=formated, alpha=0.05, weights=1)
        if (any(!ivw$outliers %in% "No significant outliers")) {
            formated_wo_outliers <- formated[!as.character(formated$SNP) %in%
                                                as.character(ivw$outliers$SNP),]
            ivw <- rmr$ivw_radial(r_input=formated_wo_outliers, alpha=0.05,
                                  weights=1)
            egger <- rmr$egger_radial(r_input=formated_wo_outliers, alpha=0.05,
                                  weights=1)
        }
        ivw$data$area <- names(combined_per_area)[test_area]
        ivw$data$type <- x
        egger$data$area <- names(combined_per_area)[test_area]
        egger$data$type <- x
        return(list(ivw=ivw, egger=egger))
    })
    names(tmp) <- c('CO', 'SV', 'HR')
    return(tmp)
})
names(mr_area) <- names(combined_per_area)

summaryMR <- lapply(seq_along(mr_area), function(test_area) {
    area <- mr_area[[test_area]]
    perX <- lapply(seq_along(area), function(test_x) {
        x <- area[[test_x]]
        ivw_tmp <- data.frame(x$ivw$coef)
        ivw_tmp$df <- x$ivw$df
        ivw_tmp$Q <- x$ivw$qstatistic
        ivw_tmp$p_q <- pchisq(ivw_tmp$Q, ivw_tmp$df, lower.tail=FALSE)
        rownames(ivw_tmp) <- c("beta_ivw")
        ivw_tmp$test <- c("beta_ivw")
        egger_tmp <- data.frame(x$egger$coef)
        egger_tmp$df <- x$egger$df
        egger_tmp$Q <- x$egger$qstatistic
        egger_tmp$p_q <- pchisq(egger_tmp$Q, egger_tmp$df, lower.tail=FALSE)
        rownames(egger_tmp) <- c("beta_1E", "beta1E")
        egger_tmp$test <- c("beta_1E", "beta1E")
        tmp <- rbind(ivw_tmp, egger_tmp)
        colnames(tmp)[1:4] <- c("Estimate", "se", "t.value", 'p.value')
        tmp$type <- names(area)[test_x]
        return(tmp)
    })
    tmp <- do.call(rbind, perX)
    tmp$area <- names(mr_area)[test_area]
    return(tmp)
})

summaryMR <- do.call(rbind, summaryMR)
summaryMR$p.adjust <- p.adjust(summaryMR$p.value)
MRsig <- summaryMR[summaryMR$p.adjust < 0.05,]

p_CO_mid <- rmr$plot_radial(r_object=mr_area$mid$CO$ivw,
                            show_outliers=FALSE, radial_scale = TRUE) +
    ggtitle("CO in mid region") + xlim(c(0, 12))

p_CO_apical <- rmr$plot_radial(r_object=mr_area$apical$CO$ivw,
                            show_outliers=TRUE, radial_scale = TRUE) +
    ggtitle("CO in apical region") + xlim(c(0, 12)) + ylim(c(-1,5))

p_SV_basal <- rmr$plot_radial(r_object=mr_area$basal$SV$ivw,
                               show_outliers=FALSE, radial_scale = TRUE) +
    ggtitle("SV in basal region") + xlim(c(0, 12)) 

p_HR_apical <- rmr$plot_radial(r_object=mr_area$apical$HR$ivw,
                              show_outliers=FALSE, radial_scale = TRUE) +
    ggtitle("HR in apical region") + xlim(c(0, 12)) + ylim(c(-0.5,5))



pvl <- select(combined, rsid, slices, slices_beta, slices_se,
                     HR_beta, HR_se, SV_beta, SV_se, CO_beta, CO_se, area)
pvl <- reshape2::melt(pvl, id.vars=c('rsid', 'slices', 'slices_se','slices_beta',
                                     'area', 'HR_beta', 'SV_beta', 'CO_beta'),
                      value.name='se', variable.name='type_se')
pvl <- reshape2::melt(pvl, id.vars=c('rsid', 'slices', 'slices_se',
                                     'slices_beta', 'area',
                                     'type_se', 'se'),
                      value.name='beta', variable.name='type_beta')
pvl$type_beta <- factor(gsub("_beta", "",as.character(pvl$type_beta)),
                        levels=c("SV", "CO", "HR"))

#reverse <- pvl$slices_beta < 0
#pvl$slices_beta[reverse] <- -pvl$slices_beta[reverse] 
#pvl$beta[reverse] <- -pvl$beta[reverse] 

p_pvl <- ggplot(data=pvl)
p_pvl <- p_pvl + geom_point(aes(x=slices_beta, y=beta, color=slices)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_color_manual(name="Slices",
                       values=c('#8c96c6','#810f7c', '#7bccc4','#43a2ca',
                                '#0868ac', '#fdbb84','#fc8d59')) +
    xlab(expression(hat(beta)[slices])) +
    ylab(expression(hat(beta)[pheno])) +
    geom_errorbar(aes(ymin=beta - se, ymax=beta + se, x=slices_beta,
                      color=slices)) +
    geom_errorbarh(aes(xmin=slices_beta-slices_se, xmax=slices_beta + slices_se,
                       y=beta, color=slices)) +
    facet_grid(type_beta~area) +
    theme_bw() +
    theme(strip.background = element_rect(fill='white', colour = 'white'))
print(p_pvl)

ggsave(plot=p_pvl, paste(directory, "/MR_PVL.pdf", sep=""),
       height=6.5, width=8)

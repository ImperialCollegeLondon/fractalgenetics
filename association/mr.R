###############################
### Libraries and functions ###
###############################
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('dplyr', attach=TRUE)
optparse <- modules::import_package('optparse')
rmr <- modules::import_package('RadialMR')

###############
## analysis ###
###############

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

## format slice betas and p-values per slice ####
slices_beta <- cbind(rsid=slices_sig$rsid,
                          slices_sig[,grepl('beta', colnames(slices_sig))])
colnames(slices_beta) <- gsub("_beta", "", colnames(slices_beta))
beta <- reshape2::melt(slices_beta, id.var='rsid', value.name='beta',
                       variable.name='slice')

slices_se <- cbind(rsid=slices_sig$rsid,
                          slices_sig[,grepl('se', colnames(slices_sig))])
colnames(slices_se) <- gsub("_se", "", colnames(slices_sig))
se <- reshape2::melt(slices_se, id.var='rsid', value.name='se',
                       variable.name='slice')

slices_logp <- cbind(rsid=slices_sig$rsid,
                          slices_sig[, grepl('log10p', colnames(slices_sig))])
colnames(slices_logp) <- gsub("\\.log10p", "", colnames(slices_logp))
logp <- reshape2::melt(slices_logp, id.var='rsid', value.name='logp',
                       variable.name='slice')

slices <- cbind(beta, logp=logp$logp, se=se$se)
slices$rsid <- as.character(slices$rsid)

## slices significant slice pvalues and betas ####
slices <- slices[slices$rsid %in% lvv$rsid, ]
sig_per_slice <- dplyr::filter(slices, logp  > -log10(5e-8))

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
    logp_tmp <- lvv_logp[pos, -1]
    names(logp_tmp) <- paste(names(logp_tmp), "_logp", sep="")
    se_tmp <- lvv_se[pos, -1]
    names(se_tmp) <- paste(names(se_tmp), "_se", sep="")
    tmp <- data.frame(beta_tmp, logp_tmp, se_tmp, stringsAsFactors=FALSE)
    return(tmp)
}))

## analyse betas slices and volumes ####
combined <- cbind(sig_per_slice, lvv_sig)
colnames(combined)[1:5] <- c('rsid', 'slices', 'slices_beta', 'slices_logp',
                           'slices_se')

base <- as.character(combined$slices) %in% c('Slice_1', 'Slice_2', 'Slice_3')
mid <- as.character(combined$slices) %in% c('Slice_4', 'Slice_5', 'Slice_6')
combined$area <- "apical"
combined$area[base] <- 'basal' 
combined$area[mid] <- 'mid' 
combined$area <- factor(combined$area, levels=c('basal', 'mid', 'apical'))

combined_per_area <- split(combined, f=combined$area)
names(combined_per_area) <- levels(combined$area)

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

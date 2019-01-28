#################
## libraries ####
#################
options(import.path=c("/homes/hannah/analysis/fd",
                      "/homes/hannah/projects"))
options(bitmapType = 'cairo', device = 'pdf')

modules::import_package('ggplot2', attach=TRUE)
autofd <- modules::import('AutoFD_interpolation')
smooth <- modules::import('utils/smoothAddR2')

exploreInterpolate <- function(interpolate, data, raw, sections) {
    FDi <- autofd$interpolate$fracDecimate(data=data, interpNoSlices=interpolate,
                               id.col.name='rownames')

    summary_interpolate <- data.frame(t(apply(as.matrix(FDi), 1,
                                              autofd$stats$summaryStatistics,
                           discard=FALSE, sections=sections)))
    summary_interpolate <- reshape2::melt(summary_interpolate[,-1])

    all <- data.frame(type=raw$variable, raw=raw$value,
                              interpolate=summary_interpolate$value,
                              slices=raw$SlicesUsed)
    all$interpolateSlices <- interpolate
    return(all)
}

#################################
## parameters and input data ####
#################################

outdir <- "~/data/digital-heart/phenotype/FD"
pheno <- "~/data/digital-heart/phenotype/FD/20181025_HVOLSegmentations_FD.csv"

NaN_values <- c("Sparse myocardium", "Meagre blood pool","FD measure failed")
NrSlices <- 7:12

################
## analysis ####
################

## FD measurements ####
dataFD <- data.table::fread(pheno, data.table=FALSE,
                            stringsAsFactors=FALSE, na.strings=c("NA", "NaN"))
rownames(dataFD) <- dataFD[, 1]
colnames(dataFD)[colnames(dataFD) == 'FD - Slice 1'] <- 'Slice 1'
dataFD <- dataFD[,grepl("Slice \\d{1,2}", colnames(dataFD))]
colnames(dataFD) <- gsub(" ", "", colnames(dataFD))

# Exclude individuals where less than 6 slices were measured
fd_notNA <- apply(dataFD, 1,  function(x) {
                length(which(!(is.na(x) | x %in% NaN_values))) > 5
                            })
dataFD <- dataFD[fd_notNA, ]


# manually look at non-numerics in FD slices
nn <- sort(unique(unlist(dataFD)), decreasing =TRUE)
dataFD <- as.data.frame(apply(dataFD, 2, function(x) {
    x[x %in% nn[2]] <- NA
    x[x %in% nn[c(1,3,4)]] <- NaN
    return(as.numeric(x))
}))

slices <- paste("Slice", 1:20, sep="")
fd_slices <- dplyr::select(dataFD, slices)

# plot distribution of nas
all_nas <- apply(fd_slices, 2, function(x) length(which(is.na(x))))
nans <- apply(fd_slices, 2, function(x) length(which(is.nan(x))))
nas <- all_nas - nans
complete <- nrow(fd_slices) - all_nas
data_na <- rbind(data.frame(Slice=1:ncol(fd_slices), samples=nas, type="NA"),
                 data.frame(Slice=1:ncol(fd_slices), samples=nans, type="NaN"),
                 data.frame(Slice=1:ncol(fd_slices), samples=complete,
                            type="Complete"))

p_na <- ggplot(data_na, aes(x=Slice, y=samples, color=type))
p_na <- p_na + geom_point(size=1) +
    facet_wrap(~type, nrow=3) +
    scale_color_brewer(type="qual", palette=6) +
    theme_bw()
ggsave(plot=p_na, file=paste(outdir, "/NAdist_FD.pdf", sep=""),
       height=4, width=4, units="in")

## Divide slices into Basal and Apical FD ####
summary_raw_BA <- data.frame(t(apply(as.matrix(dataFD), 1,
                                     autofd$stats$summaryStatistics,
                           discard=FALSE, sections="BA")))
rownames(summary_raw_BA) <- rownames(dataFD)

basalNA_BA <- summary_raw_BA[apply(summary_raw_BA[,3:4], 1, function(x)
                                    any(is.na(x))),]
apicalNA_BA <- summary_raw_BA[apply(summary_raw_BA[,5:6], 1, function(x)
                                     any(is.na(x))),]

summary_raw_BA <- reshape2::melt(summary_raw_BA, id.var="SlicesUsed")
summary_BA <- lapply(c(7,8,9,10,11,12), exploreInterpolate, data=dataFD,
                      raw=summary_raw_BA, sections="BA")
summary_BA <- do.call(rbind, summary_BA)
summary_BA <- summary_BA[!apply(summary_BA, 1, function(x) any(is.na(x))),]

p_BA <- ggplot(summary_BA, aes(x=raw, y=interpolate, color=type))
p_BA <- p_BA +
    smooth$stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE,
                            xpos=0.9, ypos=1.4, vjust=0, color="black") +
    geom_smooth(method="lm", se=FALSE) +
    geom_point(size=1) +
    facet_grid(interpolateSlices ~ type) +
    scale_color_brewer(type="qual", palette=6) +
    theme_bw() +
    theme(legend.position = "bottom",
            strip.background = element_rect(fill="white", color="white"),
            axis.text  =element_text(size=12),
            legend.title  =element_text(size=12),
            legend.text  =element_text(size=12),
            strip.text  =element_text(size=12))
ggsave(plot=p_BA,
       file=paste(outdir, "/Raw_vs_interpolated_compareBA.pdf", sep=""),
       height=15, width=15, units="in")

## Divide slices into Basal, Mid and Apical FD ####
summary_raw_BMA <- data.frame(t(apply(as.matrix(dataFD), 1,
                            autofd$stats$summaryStatistics,
                            discard=FALSE, sections="BMA")))
rownames(summary_raw_BMA) <- rownames(dataFD)

basalNA_BMA <- summary_raw_BMA[apply(summary_raw_BMA[,3:4], 1, function(x)
                                    any(is.na(x))),]
midNA_BMA <- summary_raw_BMA[apply(summary_raw_BMA[,5:6], 1, function(x)
                                     any(is.na(x))),]
apicalNA_BMA <- summary_raw_BMA[apply(summary_raw_BMA[,7:8], 1, function(x)
                                     any(is.na(x))),]
summary_raw_BMA <- reshape2::melt(summary_raw_BMA, id.var="SlicesUsed")
summary_BMA <- lapply(c(7,8,9,10,11,12), exploreInterpolate, data=dataFD,
                      raw=summary_raw_BMA, sections="BMA")
summary_BMA <- do.call(rbind, summary_BMA)
summary_BMA <- summary_BMA[!apply(summary_BMA, 1, function(x) any(is.na(x))),]

p_BMA <- ggplot(summary_BMA, aes(x=raw, y=interpolate, color=type))
p_BMA <- p_BMA +
    smooth$stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE,
                            xpos=0.9, ypos=1.4, vjust=0, color="black") +
    geom_smooth(method="lm", se=FALSE) +
    geom_point(size=1) +
    facet_grid(interpolateSlices ~ type) +
    scale_color_brewer(type="qual", palette=6) +
    theme_bw() +
    theme(legend.position = "bottom",
            strip.background = element_rect(fill="white", color="white"),
            axis.text  =element_text(size=12),
            legend.title  =element_text(size=12),
            legend.text  =element_text(size=12),
            strip.text  =element_text(size=12))
ggsave(plot=p_BMA,
       file=paste(outdir, "/Raw_vs_interpolated_compareBMA.pdf", sep=""),
       height=15, width=15, units="in")

p_BMA <- ggplot(dplyr::filter(summary_BMA, grepl("Mean", type)),
                              aes(x=raw, y=interpolate, color=type))
p_BMA <- p_BMA +
    smooth$stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE,
                            xpos=0.9, ypos=1.35, vjust=0, color="black") +
    geom_smooth(method="lm", se=FALSE) +
    geom_point(size=1) +
    facet_grid(interpolateSlices ~ type) +
    scale_color_brewer(type="qual", palette=6) +
    theme_bw() +
    theme(legend.position = "bottom",
            strip.background = element_rect(fill="white", color="white"),
            axis.text  =element_text(size=12),
            legend.title  =element_text(size=12),
            legend.text  =element_text(size=12),
            strip.text  =element_text(size=12))
ggsave(plot=p_BMA,
       file=paste(outdir, "/Raw_vs_interpolated_compareBMA_meanOnly.png",
                  sep=""),
       height=15, width=12, units="in")

raw_NA <- data.frame(apical_BMA=nrow(apicalNA_BMA), mid_BMA=nrow(midNA_BMA),
                     basal_BMA=nrow(apicalNA_BMA), apical_BA=nrow(apicalNA_BA),
                     basal_BA=nrow(basalNA_BMA)
write.table(raw_NA, file=paste(outdir, "/summary_NA.pdf", sep=""),
            col.names=TRUE, row.names=FALSE, quote=FALSE)

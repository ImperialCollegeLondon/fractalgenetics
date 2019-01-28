###############################
### Libraries and functions ###
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

modules::import_package('ggplot2', attach=TRUE)
modules::import_package('grid', attach=TRUE)
modules::import_package('dplyr', attach=TRUE)

optparse <- modules::import_package('optparse')
garfield <- modules::import_package('garfield')
zip <- modules::import_package('zip')
prepGarfield <- modules::import("prepGarfield")

## functions ####
garfieldR <- function(trait_name, directory,
                             garfielddir=paste("/nfs/research1/birney/",
                                               "resources/human_reference/",
                                               "GARFIELD/garfield-data", sep=""),
                             nperm=100000, chrs=1:22, is.NegLog10=TRUE,
                             run.option='complete'){

    resdir <- paste(directory, "/", trait_name, sep="")
    prep.file <- paste(resdir, "/garfield_results.prep", sep="")
    message("Run Garfield analyses for trait: ", trait_name)
    if (run.option %in% c('complete', 'perm', 'prep')) {
        garfieldRun <- garfield::garfield.run(out.file=paste(resdir,
                                                            "/garfield_results",
                                                            sep=""),
                             chrs=chrs, data.dir=garfielddir, trait=trait_name,
                             nperm=nperm, run.option=run.option, prep.file=prep.file)
    }
    if (run.option %in% c('complete', 'plot')) {
        message("Plot Garfield results for trait: ", trait_name)
        garfieldPlot <- garfield::garfield.plot(input_file=paste(resdir,
                                                            "/garfield_results.perm",
                                                            sep=""),
                                           num_perm=nperm,
                                           output_prefix=paste(resdir, "/",
                                                               trait_name, "_plot",
                                                               sep="")
                                           )
    }
}


prepAnalyses <- function(traitindex, bgenie, directory,
                             garfielddir=paste("/nfs/research1/birney/",
                                               "resources/human_reference/",
                                               "GARFIELD/garfield-data", sep=""),
                             is.NegLog10=TRUE,
                             annotations="1-1005"){
    trait_name <- gsub("-log10p", "", colnames(bgenie)[traitindex])
    gwas <- bgenie[, c(1:3, traitindex)]
    colnames(gwas)[4] <- "P"
    if (is.NegLog10) gwas$P <- 10^(-gwas$P)
    resdir <- paste(directory, "/", trait_name, sep="")

    message("Prepare Garfield files for trait: ", trait_name)
    garfieldPrep <- prepGarfield$prepGarfield(gwas=gwas,
            trait_name=trait_name, directory=directory,
            chr_name="chr", bp_name="pos",
            garfielddir=paste(garfielddir, "/pval", sep=""))
}

garfieldRun <- function(trait_name, garfielddir, annotations="1-1005",
                        penrichment=0.05) {
    message("Run Garfield analyses for trait: ", trait_name)
    garfieldRun <- system(paste("garfield", trait_name, garfielddir,
                                annotations, penrichment))
}

prepData <- function(input, name, link, filterTh=10, num_perm=100000,
                     category='Peaks') {
    input <- dplyr::filter(input,  NThresh >= filterTh,  Category == category)

    thresholdsP <- sort(unique(input$PThresh[which(!is.na(input$Pvalue))]))
    annotations <- unique(as.character(input$ID))
    tissues <- as.character(input$Tissue[match(annotations, input$ID)])

    enrichment <- matrix(NA, nrow=length(thresholdsP) + 1,
                         ncol=length(annotations))
    empPvalues <- matrix(NA, nrow=length(thresholdsP) + 1,
                         ncol=length(annotations))
    for (j in 1:length(annotations)) {
        for (i in 1:length(thresholdsP)) {
            enrichment[i,j] <- input$OR[which(input$ID == annotations[j] &
                                               input$PThresh == thresholdsP[i])]
            empPvalues[i,j] <- input$Pvalue[which(input$ID == annotations[j] &
                                                   input$PThresh == thresholdsP[i])]
        }
    }
    enrichment[length(thresholdsP) + 1, ] <- 1
    empPvalues[length(thresholdsP) + 1, ] <- 1
    empPvalues <- -log10(empPvalues)
    rownames(enrichment) <- c(thresholdsP, 1)
    rownames(empPvalues) <- c(thresholdsP, 1)
    colnames(enrichment) <- annotations
    colnames(empPvalues) <- annotations

    empPvalues_long <- reshape2::melt(empPvalues, value.name='empP')
    colnames(empPvalues_long)[1:2] <- c('Threshold', 'Index')

    enrichment_long <- reshape2::melt(enrichment, value.name='OR')
    colnames(enrichment_long)[1:2] <- c('Threshold', 'Index')
    enrichment_long$empP <- empPvalues_long$empP

    enrichment_long <- merge(enrichment_long, link, by='Index')
    enrichment_long$Tissue <- as.character(enrichment_long$Tissue)
    enrichment_long$Name <- name

    return(enrichment_long)
}

###########
## data ###
###########
directory <- '/homes/hannah/data/ukbb/ukb-hrt/gwas'
garfielddir <- paste("/nfs/research1/birney/resources/human_reference/",
                  "GARFIELD/garfield-data", sep="")

penrichment <- 1e-3

annotation_link <- data.table::fread(paste(garfielddir,
                                           '/annotation/link_file.txt', sep=""),
                                 data.table=FALSE, stringsAsFactors=FALSE)
peaks <- annotation_link$Index[annotation_link$Category == "Peaks"]
peaks_ranges <- "166-290,590-888"

gwas <- data.table::fread(paste(directory,
                                '/bgenie_summary_lm_st_genomewide.csv', sep=""),
                          data.table=FALSE, stringsAsFactors=FALSE)
gwas$SNPID <- paste(gwas$chr, ":", gwas$pos, sep="")

index_logp <- which(grepl("log10p", colnames(gwas)))
traits <- gsub("-log10p", "", colnames(gwas)[index_logp])

###############
## analysis ###
###############

## prepare garfield data ####
perSlicePrep <- sapply(index_logp, prepAnalyses, bgenie=gwas,
                        annotations=peaks_ranges,
                        directory=directory)

## submit garfield jobs ####
summaryGarfieldprep <- clustermq::Q(garfieldR,
                                 trait_name=traits[c(2,4,6)][c(2,4,6)],
                                 const=list(garfielddir=garfielddir,
                                            directory=directory,
                                            run.option='prep'),
                                 n_jobs=3, memory=20000)

summaryGarfieldperm <- clustermq::Q(garfieldR,
                                 trait_name=traits[c(2,4,6)],
                                 const=list(garfielddir=garfielddir,
                                            directory=directory,
                                            run.option='perm'),
                                 n_jobs=3, memory=20000)

summaryGarfieldplot <- clustermq::Q(garfieldR,
                                 trait_name=traits[c(2,4,6)],
                                 const=list(garfielddir=garfielddir,
                                            directory=directory,
                                            run.option='plot'),
                                 n_jobs=3, memory=1000)
## submit garfield jobs ####
perSummaryGarfield <- clustermq::Q(garfieldRun,
                                 trait_name=traits[c(2,4,6)],
                                 const=list(garfielddir=garfielddir,
                                            annotations=peaks_ranges,
                                            penrichment=penrichment),
                                 n_jobs=3, memory=50000)

## read garfield results ####
basalFD <- data.table::fread(paste(garfielddir, '/output/MeanBasalFD/',
                                 'garfield.test.MeanBasalFD.out', sep=""),
                           data.table=FALSE, stringsAsFactors=FALSE)

midFD <- data.table::fread(paste(garfielddir, '/output/MeanMidFD/',
                                 'garfield.test.MeanMidFD.out', sep=""),
                           data.table=FALSE, stringsAsFactors=FALSE)

apicalFD <- data.table::fread(paste(garfielddir, '/output/MeanApicalFD/',
                                 'garfield.test.MeanApicalFD.out', sep=""),
                           data.table=FALSE, stringsAsFactors=FALSE)

apical_variants <- data.table::fread(paste(garfielddir, '/output/MeanApicalFD/',
                                 'garfield.test.MeanApicalFD.out.significant.annotations.1e-5.0.001.variants', sep=""),
                           data.table=FALSE, stringsAsFactors=FALSE)

apical_variants$SNPID <- gsub("(\\d{1,2}:\\d*)\\(.*", "\\1", apical_variants$VAR_INFO)
apical_variants <- merge(apical_variants, annotation_link, by.x='ID', by.y="Index")
apical_fh <- dplyr::filter(apical_variants, Tissue %in% "fetal_heart")

gwas_apical <- gwas[gwas$SNPID %in% apical_fh$SNPID,]
apical_results <- dplyr::select(gwas_apical, SNPID, rsid, chr, pos, af,
                                'MeanApicalFD_beta', 'MeanApicalFD-log10p')

mid_variants <- data.table::fread(paste(garfielddir, '/output/MeanMidFD/',
                                 'garfield.test.MeanMidFD.out.significant.annotations.1e-5.0.001.variants', sep=""),
                           data.table=FALSE, stringsAsFactors=FALSE)

mid_variants$SNPID <- gsub("(\\d{1,2}:\\d*)\\(.*", "\\1", mid_variants$VAR_INFO)
mid_variants <- merge(mid_variants, annotation_link, by.x='ID', by.y="Index")
mid_fh <- dplyr::filter(mid_variants, Tissue %in% "fetal_heart")

gwas_mid <- gwas[gwas$SNPID %in% mid_fh$SNPID,]
mid_results <- dplyr::select(gwas_mid, SNPID, rsid, chr, pos, af,
                                'MeanMidFD_beta', 'MeanMidFD-log10p')

## format garfield results ####
mid <- prepData(input=midFD, link=annotation_link, name='Mid')
apical <- prepData(input=apicalFD, link=annotation_link, name='Apical')
basal <- prepData(input=basalFD, link=annotation_link, name='Basal')

combined <- rbind(basal, mid, apical)
combined$Name <- factor(combined$Name, levels=c("Basal", "Mid", "Apical"))


## select tissues of interest and represeentative colors
tissues_color <- c("tomato", "skyblue3", "yellow", "brown2", "lightgreen",
                   "lightgoldenrod3", "purple", "pink", "darkblue", "gray",
                   "darkgreen")
toi <- c("fetal_heart", "heart", 'fetal_muscle', 'muscle',
         "blood", "blood_vessel", 'epithelium')
toi_color <- c( '#542788', '#8073ac', '#4575b4', '#74add1',
                '#e6f598' ,'#abdda4', '#66c2a5', '#666666')
all_color <- colorRampPalette(toi_color)(length(unique(combined$Tissue)))
section_color <- c('#fdcc8a','#fc8d59','#e34a33')

## depict region-wise annotation enrichments ####
# a) only tissues of interest
selected <- combined
selected$Tissue[!selected$Tissue %in% toi] <- 'other'
selected$Tissue <- factor(selected$Tissue, levels=c(toi, 'other'),
                          labels=c(gsub("_", " ", as.character(toi)),
                                          'other tissues'))
selected_red <- dplyr::filter(selected, Threshold == 1e-5, empP > -log10(5e-3))

p_selected <- ggplot(selected_red, aes(x=Tissue, y=OR, fill=Tissue))
p_selected <- p_selected +
    facet_grid(~Name, scales = "free_x", space="free_x", switch = "x") +
    geom_boxplot() +
    scale_fill_manual(values=toi_color, name="Tissue") +
    theme_bw() +
    theme(legend.position='bottom',
          strip.background=element_rect(fill='white'),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
    )
g_selected <- ggplot_gtable(ggplot_build(p_selected))
strip <- which(grepl('strip-b', g_selected$layout$name))
k <- 1
for (i in strip) {
    j <- which(grepl('rect', g_selected$grobs[[i]]$grobs[[1]]$childrenOrder))
    g_selected$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- section_color[k]
    k <- k+1
}
grid.draw(g_selected)
ggsave(plot=g_selected,
       filename=paste(directory, '/annotation/Funtional_enrichment.pdf', sep=""),
       height=5, width=8)


# b) all GARFIELD tissues
all_red <- dplyr::filter(combined, Threshold == 1e-5, empP > -log10(5e-2))

p_all <- ggplot(all_red, aes(x=Tissue, y=OR, fill=Tissue))
p_all <- p_all +
    facet_wrap(~Name, scales = "free_x", strip.position = 'top',
               ncol=1) +
    geom_boxplot() +
    scale_fill_manual(values=all_color, name="Tissue") +
    theme_bw() +
    theme(legend.position='bottom',
          strip.background=element_rect(fill='white'),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
    )


g_all <- ggplot_gtable(ggplot_build(p_all))
strip <- which(grepl('strip-b', g_all$layout$name))
k <- 1
for (i in strip) {
    j <- which(grepl('rect', g_all$grobs[[i]]$grobs[[1]]$childrenOrder))
    g_all$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- section_color[k]
    k <- k+1
}
grid.draw(g_all)
ggsave(plot=g_all,
       filename=paste(directory, '/annotation/Funtional_enrichment_all.pdf',
                      sep=""),
       height=12, width=9)

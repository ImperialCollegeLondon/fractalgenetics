###############################
### Libraries and functions ###
###############################

library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(grid)

prepData <- function(input, name, link, filterTh=10, num_perm=100000,
                     category='Peaks') {
    input <- dplyr::filter(input,  NThresh >= filterTh,  Category == category)
    input$FE[which(input$FE == -1)] <-  0
    input$EmpPval[which(input$EmpPval == 0)] <- 1/num_perm
    input$EmpPval[which(input$EmpPval == -1)] <- NA
    
    thresholdsP <- sort(unique(input$PThresh[which(!is.na(input$EmpPval))]))
    annotations <- unique(as.character(input$Index))
    tissues <- as.character(input$Tissue[match(annotations, input$Index)])
    
    enrichment <- matrix(NA, nrow=length(thresholdsP) + 1,
                         ncol=length(annotations))
    empPvalues <- matrix(NA, nrow=length(thresholdsP) + 1,
                         ncol=length(annotations))
    for (j in 1:length(annotations)) {
        for (i in 1:length(thresholdsP)) {
            enrichment[i,j] <- input$FE[which(input$Index == annotations[j] & 
                                               input$PThresh == thresholdsP[i])]
            empPvalues[i,j] <- input$EmpPval[which(input$Index == annotations[j] & 
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
    
    enrichment_long <- reshape2::melt(enrichment, value.name='FE')
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
directory <- '~/data/ukbb/ukb-hrt/gwas'

annotation_link <- data.table::fread('~/data/GARFIELD/garfield-data/annotation/link_file.txt',
                                 data.table=FALSE, stringsAsFactors=FALSE)
midFD <- data.table::fread(paste(directory, '/MeanMidFD/garfield_results.perm', 
                                 sep=""), data.table=FALSE,
                           stringsAsFactors=FALSE)
apicalFD <- data.table::fread(paste(directory,
                                    '/MeanApicalFD/garfield_results.perm',
                                    sep=""), data.table=FALSE,
                              stringsAsFactors=FALSE)
basalFD <- data.table::fread(paste(directory,
                                   '/MeanBasalFD/garfield_results.perm',
                                   sep=''), data.table=FALSE,
                             stringsAsFactors=FALSE)

tissues_color <- c("tomato", "skyblue3", "yellow", "brown2", "lightgreen",
                   "lightgoldenrod3", "purple", "pink", "darkblue", "gray",
                   "darkgreen")
toi <- c("fetal_heart", "heart", 'fetal_muscle', 'muscle',
         "blood", "blood_vessel", 'epithelium')
toi_color <- c( '#e31a1c', '#fb9a99', '#08519c', '#1f78b4',
                   '#984ea3', '#ff7f00','#ffff33', '#666666')
section_color <- c('#fdcc8a','#fc8d59','#e34a33')

###############
## analysis ###
###############

mid <- prepData(input=midFD, link=annotation_link, name='Mid')
apical <- prepData(input=apicalFD, link=annotation_link, name='Apical')               
basal <- prepData(input=basalFD, link=annotation_link, name='Basal')   

combined <- rbind(basal, mid, apical)
combined$Name <- factor(combined$Name, levels=c("Basal", "Mid", "Apical"))

#toi <- c("fetal_heart", "heart", 'fetal_muscle', 'muscle', 'myometrium',
#         "blood", "blood_vessel", 'epithelium',   'blastula')
#tissue_color <- c( '#e31a1c', '#fb9a99', '#08519c', '#1f78b4', '#a6cee3',
 #                  '#984ea3', '#ff7f00','#ffff33','#a65628', '#666666')


selected <- combined
selected$Tissue[!selected$Tissue %in% toi] <- 'other'
selected$Tissue <- factor(selected$Tissue, levels=c(toi, 'other'),
                          labels=c(gsub("_", " ", as.character(toi)),
                                          'other tissues'))
selected_red <- dplyr::filter(selected, Threshold == 1e-6, empP > -log10(1e-3))

p_selected <- ggplot(selected_red, aes(x=Tissue, y=FE, fill=Tissue))
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
    
all_red <- dplyr::filter(combined, Threshold == 1e-6, empP > -log10(1e-3))

p_all <- ggplot(all_red, aes(x=Tissue, y=FE, fill=Tissue))
p_all <- p_all +
    facet_wrap(~Name, scales = "free_x", strip.position = 'bottom',
               ncol=1) +
    geom_boxplot() +
    scale_fill_manual(values=colorRampPalette(toi_color)(length(unique(combined$Tissue))),
                      name="Tissue") +
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

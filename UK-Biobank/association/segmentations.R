options(import.path="~/analysis/fractalgenetics/fractal-analysis-processing/")
rr <- modules::import('radial-registration')
modules::import_package('ggplot2', attach=TRUE)
modules::import_package('plotrix', attach=TRUE)
makeCall <- function(genotype) {
    x <- rep(1, length(genotype))
    x[genotype < 0.5] <- 0
    x[genotype >= 1.5] <- 2
    return(x)
}
############
## data ####
############
directory <- "~/data/ukbb/ukb-hrt"
nonseg <- c("Sparse myocardium",  "Meagre blood pool", "FD measure failed")

coregistered <- read.table(paste(directory,
                                 "/segmentations/BiobankCoregisteredIDs.txt",
                                 sep=""), stringsAsFactors=FALSE)[,1]

## Loci associated with FD ####
sig_loci <- read.table(paste(directory, '/gwas/Significant_per_slice.csv',
                             sep=''),
                         header=TRUE, sep=',')
sig_loci <- sig_loci[!sig_loci$rsid %in% c('rs71394376', 'rs71105784'),]

## genotypes of associated loci ####
genotypes <- data.table::fread(paste(directory, "/gwas/",
                                     "Pseudomultitrait_slices_sig5e08_genotypes.dosage",
                                     sep=""),
                               stringsAsFactors=FALSE, data.table=FALSE)
geno <- genotypes[,-c(1,2,4,5,6)]
geno <- dplyr::filter(geno, rsid %in% sig_loci$rsid)
geno <- geno[!duplicated(geno$rsid),]
rownames(geno) <- geno$rsid
geno <- as.data.frame(t(geno[,-1]))
geno$IID <- rownames(geno)

ukb_ids <- dplyr::select(geno, IID)
ukb_geno <- dplyr::select(geno, as.character(unique(sig_loci$rsid)))
ukb_called <- as.data.frame(apply(ukb_geno, 2, makeCall))
ukb <- cbind(ukb_ids, ukb_called)
ukb <- dplyr::filter(ukb, IID %in% coregistered)

## FD measurements ####
slices <- data.table::fread(paste(directory, '/rawdata/FD.csv', sep=''),
                            data.table=FALSE, stringsAsFactors=FALSE)
slices <- slices[, c(1, 10:29)]
colnames(slices) <- c("IID", paste('Slice_', 1:20, sep=""))
slices <- data.frame(t(apply(slices, 1, function(x) {
    x[x %in% c("NaN", "NA")] <- NA
    return(x)
})), stringsAsFactors = FALSE)
slices <- dplyr::filter(slices, IID %in% ukb$IID)

## available image list per individual ####
folders <- list.dirs(paste(directory, "/segmentations/aligned", sep=""))
folders <- folders[-1]
filelist <- list.files(folders, pattern="Edge-Image-Slice-.*-ED.png",
                       full.names=TRUE)
idlist <- gsub(".*/", "", folders)

################
## analysis ####
################

## individuals with same number of measured slices ####
same_slices <- slices[apply(slices, 1, function(x) !any(x[3:10] %in% nonseg)),
                      c(1,3:10)]

## measurements and genotypes per slice ####
sig_per_slice <- split(dplyr::select(sig_loci, rsid, slice), f=sig_loci$slice)
ukb_per_slice <- lapply(sig_per_slice, function(x) {
    pos <- c(1,which(colnames(ukb) %in% x$rsid))
    ukb[,pos]
})

geno_per_slice <- lapply(seq_along(ukb_per_slice), function(x){
    s <- names(ukb_per_slice)[x]
    iid <- ukb_per_slice[[x]]$IID
    g <- ukb_per_slice[[x]][, -1]
    fd <- dplyr::select(slices, IID, s)
    tmp <- lapply(g, function(gs) {
        new <- fd
        new$geno <- 0
        new$geno[new$IID %in% iid[gs == 1]] <- 1
        new$geno[new$IID %in% iid[gs == 2]] <- 2
        return(new)
    })
    return(tmp)
})
names(geno_per_slice) <- names(ukb_per_slice)
saveRDS(geno_per_slice,
        paste(directory, "/genotypes_per_slice_association.rds", sep=""))

## median FD per slice and genotype
median_per_slice <- sapply(geno_per_slice, function(s) {
    lapply(s, function(g) {
        g[,2] <- as.numeric(g[,2])
        gg <- split(g, f=as.factor(g$geno))
        median_id <- data.frame(t(sapply(gg, function(x) {
            pos <- which.min(abs(x[,2] - median(x[,2], na.rm = TRUE)))
            cc <- c(x$IID[pos], x[pos,2])
            return(cc)
        })), stringsAsFactors=FALSE)
        if (length(gg[[1]]) < length(gg[[2]])) {
            median_id <- median_id[nrow(median_id):1,]
        }
        names(median_id) <- c("IID", "FD")
        median_id$FD <- as.numeric(median_id$FD)
        return(median_id)
    })
})

## read, process and register all measured slices per individual ####
processed <- parallel::mclapply(idlist, rr$processSlices,
#processed <- lapply(idlist, rr$processSlices,
                                directory=paste(directory,
                                                "/segmentations/aligned",
                                                sep=""),
                                sliceInfo=slices, n=2,
                    rotation=FALSE, mc.cores=8)
names( processed) <- lapply(processed, function(x) {
    if (length(x) == 1) return(NULL)
    if (length(x) != 1) return(unique(x[[1]]$registered$perimeter$IID))
    })
#saveRDS(processed, paste(directory, "/registered_segmentations.rds", sep=""))

## extract measurments of relevant individuals ####
ioi <- processed[names(processed) %in% same_slices$IID]

edges_all <- do.call(rbind, lapply(processed, function(ind) {
    do.call(rbind, lapply(ind, function(x) x$registered$edges))
}))

edges <- dplyr::filter(edges_all, IID %in% same_slices$IID)

tmp <- lapply(idlist, function(id) {
    dplyr::filter(edges, slice %in% "Slice_4", IID %in% id)
})
r.transformed <- do.call(cbind, lapply(tmp, function(x) x$r.transformed))
colnames(r.transformed) <- sapply(tmp, function(x) unique(x$IID))
theta <- do.call(cbind, lapply(tmp, function(x) x$theta))
colnames(theta) <- sapply(tmp, function(x) unique(x$IID))

edges_slice5 <- list(r=r.transformed, theta=theta)
median_slice5 <- median_per_slice[[3]][[1]]

#pdf(paste(directory, "/segmentations/Slice5_rs35006907.pdf", sep=""),
pdf(paste(directory, "/segmentations/Slice4_rs17608766.pdf", sep=""),
    height=6, width=2)
layout(matrix(c(1, 2, 3), 3))
radial.plot(edges_slice5$r[,1],
            radial.pos=edges_slice5$theta[,1],
            labels="", main=median_slice5$FD[1],
            rp.type = "s", point.symbols = 20, cex=0.5,
            radial.lim=c(0, 1), point.col='#386cb0',
            show.grid.labels = FALSE)
radial.plot(edges_slice5$r[,2],
            radial.pos=edges_slice5$theta[,2],
            labels="", rp.type = "s", main=median_slice5$FD[2],
            radial.lim=c(0, 1), point.symbols = 20, cex=0.5,
            point.col='#e7298a',add=FALSE,
            show.grid.labels = FALSE)
radial.plot(edges_slice5$r[,3],
            radial.pos=edges_slice5$theta[,3],
            labels="", main=median_slice5$FD[3],
            rp.type = "s", point.symbols = 20, cex=0.5,
            radial.lim=c(0, 1), add=FALSE, point.col='#66a61e',
            show.grid.labels = FALSE)
dev.off()



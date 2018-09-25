if (FALSE) {
    args <- list()
    args$directory <-"/homes/hannah/data/ukbb/ukb-hrt/gwas"
    args$name <-"slices"
    args$valueMeff <- NULL
    args$phenofile <- "/homes/hannah/data/ukbb/ukb-hrt/phenotypes/FD_slices_EUnorel.csv"
    args$verbose <- TRUE
    args$interpolate <- 9
}

directory <- args$directory
name <- args$name
phenofile <- args$phenofile
verbose <- args$verbose

genomewide <- data.table::fread(file=paste(directory, "/bgenie_", name,
                                           "_lm_st_genomewide.csv", sep=""),
                                data.table=FALSE, stringsAsFactors=FALSE)

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))

d$pos <- NA
ticks <- NULL
lastbase <- 0
numchroms <- length(unique(d$CHR))

if (numchroms == 1) {
	d$pos <- d$BP
} else {
	for (i in unique(d$CHR)) {
		if (i == 1) {
			d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP
		} else {
			lastbase <- lastbase + max(subset(d, CHR==i-1)$BP)
			d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP + lastbase
		}
		ticks <- c(ticks,
				d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
	}
	ticklim <- c(min(d$pos),max(d$pos))
}



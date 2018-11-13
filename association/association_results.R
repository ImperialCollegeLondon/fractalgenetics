###############################
### Libraries and functions ###
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('optparse')
garfield <- modules::import_package('garfield')
zip <- modules::import_package('zip')
prepGarfield <- modules::import("prepGarfield")
bgenie <- modules::import("bgenieResults")
plots <- modules::import("plots")
meta <- modules::import("metaanalysis")
ldfilter <- modules::import("LDfilter")

stPlots <- function(trait_index, bgenie_result, directory, is.negLog=FALSE,
                         Meff=1, ymin=1, ymax="max", name=NULL) {
    if (!is.null(name)) name <- paste("_", name, sep="")
    if (!is.negLog) {
        pvalues <- sapply(bgenie_result[,trait_index],
                          function(x) min(x * valueMeff, 1))
        which.sig <- pvalues < ymin
    } else {
        pvalues <- sapply(bgenie_result[,trait_index],
                          function(x) max(x - log10(valueMeff), 0))
        which.sig <- pvalues > -log10(ymin)
    }
    p_manhattan <- plots$manhattan(data.frame(bgenie_result[which.sig, 1:3],
                                              P=pvalues[which.sig]),
                                   title=colnames(bgenie_result)[trait_index],
                                   size.x.labels=12, size.y.labels=12,
                                   is.negLog=is.negLog,
                                   color=c("#fc8d59", "#b30000"), chr='chr',
                                   snp="rsid", bp="pos", max.y=ymax)
    ggplot2::ggsave(plot=p_manhattan, height=4, width=10,
                    file=paste(directory,"/bgenie", name, "_lm_",
                               colnames(bgenie_result)[trait_index],
                               "_manhattanplot.pdf", sep=""))
    p_qq <- plots$qqplot(bgenie_result[,trait_index], size.text=14,
                      size.title=14, is.negLog=is.negLog, raster=FALSE)
    ggplot2::ggsave(plot=p_qq, height=7, width=7,
                  file=paste(directory,"/bgenie", name, "_lm_",
                         colnames(bgenie_result)[trait_index],"_qqplot.png",
                    sep=""))
}

garfieldAnalyses <- function(traitindex, bgenie, directory,
                             garfielddir=paste("/nfs/research1/birney/",
                                               "resources/human_reference/",
                                               "GARFIELD/garfield-data", sep=""),
                             nperm=100000, chrs=1:22, is.NegLog10=TRUE){
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

    message("Run Garfield analyses for trait: ", trait_name)
    garfieldRun <- garfield$garfield.run(out.file=paste(resdir,
                                                        "/garfield_results",
                                                        sep=""),
                         chrs=chrs, data.dir=garfielddir, trait=trait_name,
                         nperm=nperm)

    message("Plot Garfield results for trait: ", trait_name)
    garfieldPlot <- garfield$garfield.plot(input_file=paste(resdir,
                                                        "/garfield_results.perm",
                                                        sep=""),
                                       num_perm=nperm,
                                       output_prefix=paste(resdir, "/",
                                                           trait_name, "_plot",
                                                           sep="")
                                       )
}


###############
## analysis ###
###############


## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="directory",
               type="character", help="Path to directory with bgenie association
               results [default: %default].", default=NULL),
    optparse$make_option(c("-n", "--name"), action="store", dest="name",
               type="character", help="Name of analysis; has to be the same as
               in naming bgenie files [default: %default].", default=NULL),
    optparse$make_option(c("-p", "--pheno"), action="store", dest="phenofile",
               type="character", help="Path to fd phenotype file; if mEffective
               is not provided, --pheno can be used to estimate --mEffective. If
               neither are set, default mEffective = NrTraits_tested
               [default:%default].", default=NULL),
    optparse$make_option(c("-m", "--mEffective"), action="store",
               dest="valueMeff",
               type="double", help="Estimated number of effective test; based
               on Galwey (2009) Genetic Epidemiology [default:
               NrTraits_tested].", default=NULL),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out [default: %default]."),
    optparse$make_option(c("-i", "--interpolate"), action="store",
               dest="interpolate",
               type="integer", help="Number of slices to interpolate to
               [default: %default].", default=9)
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

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
tags_prefix <- "European_ukb_imp_chr"
tags_suffix <- "_v3_maf0.001_250kb_r0.6.tags.list"
tags_dir <- "~/data/ukbb/ukb-hrt/tags"

## Read results ####
if (verbose) message("Read files with genome-wide association results")
genomewide <- lapply(1:22, bgenie$readBgenieOutput, directory=directory,
                     name=paste("bgenie_", name, "_lm_st_chr", sep=""),
                     maf=0.001, biallelicOnly=FALSE)
genomewide <- do.call(rbind, genomewide)

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))
index_beta <- which(grepl("beta", colnames(genomewide)))
index_t <- which(grepl("_t", colnames(genomewide)))

## Write combined results ####
if (verbose) message("Write file with genome-wide association results")
write.table(genomewide,
            file=paste(directory, "/bgenie_", name, "_lm_st_genomewide.csv",
                       sep=""),
            sep=",",quote=FALSE, col.names=TRUE, row.names=FALSE)

## Meta-analysis single-trait summary statistics ####
if (name %in% c('summary', 'slices')) {
    if (verbose) message("Estimate pseudo-multitrait association results")
    multitrait_res <- meta$pseudoMultitrait(genomewide[, index_t])
    multitrait <- cbind(genomewide[,index_snp], multitrait_res$multitrait_p)
    colnames(multitrait) <- c("CHR", "SNP","BP", "P")

    if (verbose) message("Write pseudo-multitrait association results")
    write.table(multitrait,
                file=paste(directory, "/bgenie_", name, "_lm_pseudomt_all.csv",
                           sep=""),
                sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
    if (name == 'slices') {
        title="Pseudomulti-trait all FD slices"
    } else {
        title="Pseudomulti-trait mean basal, mid and apical FD"
    }

    if (verbose) message("Plot pseudo-multitrait association results")
    ymin <- 10^(-1)
    ymax <- max(c(-log10(min(multitrait$P)), max(genomewide[,index_logp])))
    sig <- multitrait[multitrait$P < ymin, ]
    p_manhattan_mt <- plots$manhattan(d=sig,
                        title=title,
                        size.x.labels=12, size.y.labels=12, is.negLog=FALSE,
                        color=c("#fc8d59", "#b30000"), max.y=ymax)

    ggplot2::ggsave(plot=p_manhattan_mt, height=4, width=10,
           file=paste(directory,"/bgenie_", name,
                      "_lm_pseudomt_manhattanplot.pdf", sep=""))

    p_qq_mt <- plots$qqplot(multitrait$P, size.text=14, size.title=14)
    ggplot2::ggsave(plot=p_qq_mt, height=7, width=7,
           file=paste(directory,"/bgenie_", name,
                      "_lm_pseudomt_qqplot.pdf", sep=""))

    if (verbose) message("Filter pseudo-multitrait association results for LD")
    multitrait_ld <-
        ldfilter$filterSigLoci4LD(gwas=cbind(multitrait[,1:3], P=multitrait$P,
                                             genomewide[,c(4:7, index_beta)]),
                                  tags_dir=tags_dir,
                                  tags_prefix=tags_prefix,
                                  tags_suffix=tags_suffix,
                                  matchBy="SNP", tags_ID="SNP",
                                  tagsSplit="|")
    ldfilter$writeSig(multitrait_ld, threshold=5*10^(-8), directory=directory,
                   name="Pseudomultitrait_Slices")
    sig <- genomewide[genomewide$rsid %in% multitrait_ld$sig_no_ld,]
    write.table(sig, paste(directory, "/", name, "_sig5e08_ldFiltered.txt",
                                  sep=""),
                col.names=TRUE, row.names=FALSE, quote=FALSE)

    if (name == 'summary') {
        ## GARFIELD analysis ####
        perTraitGarfield <- sapply(index_logp, garfieldAnalyses,
                                   bgenie=genomewide, directory=directory)
    }
}
## per trait qq and manhattan plots ####
# effective number of tests to adjust single-trait assocation p-values by
if (!is.null(phenofile)) {
    pheno <- data.table::fread(phenofile, data.table=FALSE,
                               stringsAsFactors=FALSE)
    valueMeff <- meta$Teff(as.matrix(pheno[,-1]))
} else if (is.null(args$valueMeff)) {
    valueMeff <- length(index_logp)
} else {
    valueMeff <- args$valueMeff
}

# plot single-trait GWAS
if (verbose) message("Plot single-trait association results")
plots_perTrait <- sapply(index_logp, stPlots, bgenie_result=genomewide, ymin=0.1,
                         is.negLog=TRUE, directory=directory, Meff=valueMeff,
                         name=name,  ymax=ymax)

## Format single-trait association statistics for LD score regression ####
if (verbose) message("Format results for LD score regression")
ldshub_file <- "~/data/ukbb/ukb-hrt/gwas/ldc_hub/w_hm4.noMHC.snplist"
ldshub_snps <- data.table::fread(ldshub_file, data.table=FALSE,
                                 stringsAsFactors=FALSE)

ldsc_single <- sapply(index_beta, function(x) {
                   columns <- c(1:7, x:(x+3))
                   nn <- paste(name, "_",
                               gsub("_beta", "", colnames(genomewide)[x]),
                               sep="")
                   tmp <- bgenie$bgenie2ldsc(genomewide[,columns],
                                             sumstat="_beta",
                                             ldshub_snps=ldshub_snps)
                   message("Write LDSC output for ", nn)
                   write.table(tmp, file=paste(directory, "/sumstats_", nn,
                                               ".txt", sep=""),
                               sep="\t", quote=FALSE, col.names=TRUE,
                               row.names=FALSE)
                   zip$zip(zipfile=paste(directory, "/sumstats_", nn,".gz",
                                         sep=""),
                       paste(directory, "/sumstats_", nn, ".txt", sep=""))
                }
)



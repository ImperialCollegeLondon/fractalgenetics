###############################
### Libraries and functions ###
###############################
options(import.path="/homes/hannah/projects/GWAS")
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('optparse')
garfield <- modules::import_package('garfield')
prepGarfield <- modules::import("prepGarfield")
bgenie <- modules::import("bgenieResults")
plots <- modules::import("plots")
meta <- modules::import("metaanalysis")

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

garfieldAnalyses <- function(traitindex, bgenie, directory, garfielddir,
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
    args$name <-"summary"
    args$valueMeff <- NULL
    args$phenofile <- "/homes/hannah/data/ukbb/ukb-hrt/phenotypes/FD_slices_EUnorel.csv"
    args$verbose <- TRUE
    args$interpolate <- 9
}

directory <- args$directory
name <- args$name
phenofile <- args$phenofile
verbose <- args$verbose

genomewide <- lapply(1:22, bgenie$readBgenieOutput, directory=directory,
                     name=paste("bgenie_", name, "_lm_st_chr", sep=""),
                     maf=0.001, biallelicOnly=FALSE)
genomewide <- do.call(rbind, genomewide)

index_snp <- 1:3
index_logp <- which(grepl("log10p", colnames(genomewide)))
index_beta <- which(grepl("beta", colnames(genomewide)))
index_t <- which(grepl("_t", colnames(genomewide)))

## write results ####
if (verbose) message("Write file with genome-wide association results")
write.table(genomewide,
            file=paste(directory, "/bgenie_", name, "_lm_st_genomewide.csv",
                       sep=""),
            sep=",",quote=FALSE, col.names=TRUE, row.names=FALSE)

## Meta-analysis single-trait summary statistics ####
if (verbose) message("Estimate pseudo-multitrait association results")
multitrait_res <- meta$pseudoMultitrait(genomewide[, index_t])
multitrait <- cbind(genomewide[,index_snp], multitrait_res$multitrait_p)
colnames(multitrait) <- c("CHR", "SNP","BP", "P")

if (verbose) message("Write pseudo-multitrait association results")
write.table(multitrait,
            file=paste(directory, "/bgenie_", name, "_lm_pseudomt_all.csv",
                       sep=""),
            sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

if (verbose) message("Plot pseudo-multitrait association results")
ymin <- 10^(-1)
ymax <- max(c(-log10(min(multitrait$P)), max(genomewide[,index_logp])))
sig <- multitrait[multitrait$P < ymin, ]
p_manhattan_mt <- plots$manhattan(d=sig,
                    title="Pseudomulti-trait all",
                    size.x.labels=12, size.y.labels=12, is.negLog=FALSE,
                    color=c("#fc8d59", "#b30000"), max.y=ymax)

ggplot2::ggsave(plot=p_manhattan_mt, height=4, width=10,
       file=paste(directory,"/bgenie_", name,
                  "_lm_pseudomt_manhattanplot.pdf", sep=""))

p_qq_mt <- plots$qqplot(multitrait$P, size.text=14, size.title=14)
ggplot2::ggsave(plot=p_qq_mt, height=7, width=7,
       file=paste(directory,"/bgenie_", name,
                  "_lm_pseudomt_qqplot.pdf", sep=""))

if (grepl("summary", name)) {
    if (verbose) message("Estimate pseudo-multitrait association results for
                          mean basal, mean mid and mean apical FD")
    multitrait_meanBasalMidApical <-
        meta$pseudoMultitrait(genomewide[, index_t][,c(2,4,6)])
    multitrait_meanBasalMidApical <- cbind(genomewide[,index_snp],
        multitrait_meanBasalMidApical$multitrait_p)
    colnames(multitrait_meanBasalMidApical) <- c("CHR", "SNP","BP", "P")
    if (verbose) message("Write pseudo-multitrait association results for
                          max basal and max apical FD")
    write.table(multitrait_meanBasalMidApical,
                file=paste(directory, "/bgenie_", name,
                           "_lm_pseudomt_meanBasalMidApical.csv", sep=""),
                sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

    if (verbose) message("Plot pseudo-multitrait association results for
                          mean basal, mean mid and mean apical FD")
     sig <- multitrait_meanBasalMidApical[multitrait_meanBasalMidApical$P < ymin, ]
     p_manhattan_mt <- plots$manhattan(d=sig,
                    title="Pseudomulti-trait mean basal, mid and apical FD",
                    size.x.labels=12, size.y.labels=12, is.negLog=FALSE,
                    color=c("#fc8d59", "#b30000"),  max.y=ymax)
     ggplot2::ggsave(plot=p_manhattan_mt, height=4, width=10,
            file=paste(directory,"/bgenie_", name,
                 "_lm_pseudomt_manhattanplot_meanBasalMidApical.pdf", sep=""))

     p_qq_mt <- plots$qqplot(multitrait_meanBasalMidApical$P, size.text=14,
                             size.title=14, is.negLog=FALSE, raster=FALSE)
     ggplot2::ggsave(plot=p_qq_mt, height=7, width=7,
       file=paste(directory,"/bgenie_", name,
                  "_lm_pseudomt_qqplot_meanBasalMidApical.png",
                  sep=""))
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
ldsc_single <- sapply(index_beta, function(x) {
                   columns <- c(1:7, x:(x+3))
                   nn <- paste(name, "_",
                               gsub("_beta", "", colnames(genomewide)[x]),
                               sep="")
                   tmp <- bgenie$bgenie2ldsc(genomewide[,columns],
                                             sumstat="_beta")
                   message("Write LDSC output for ", nn)
                   write.table(tmp, file=paste(directory, "/sumstats_", nn,
                                               ".txt", sep=""),
                               sep="\t", quote=FALSE, col.names=TRUE,
                               row.names=FALSE)
                         }
)

multitrait_ldsc <-
    bgenie$bgenie2ldsc(cbind(genomewide[,1:7], chi2=multitrait_res$chi2stats,
                             "-log10p"=-log10(multitrait_res$multitrait_p)),
                       sumstat="chi2")
write.table(multitrait_ldsc, file=paste(directory, "/sumstats_", name,
                                        "_pseudomt.txt", sep=""),
                               sep="\t", quote=FALSE, col.names=TRUE,
                               row.names=FALSE)

## GARFIELD analysis ####
trait_name <- "meanBasalMidApicalFD"
gdir <- "/homes/hannah/software/garfield-data"
resdir <- paste(directory, "/", trait_name, sep="")
nperm <- 100000

perTraitGarfield <- sapply(index_logp[c(2,4,6)], garfieldAnalyses,
                           bgenie=genomewide, directory=directory,
                           garfielddir=gdir)




    garfieldPrep <- prepGarfield$prepGarfield(gwas=multitrait_meanBasalMidApical,
            trait_name="meanBasalMidApicalFD", directory=directory,
            garfielddir=paste(gdir, "/pval", sep=""))




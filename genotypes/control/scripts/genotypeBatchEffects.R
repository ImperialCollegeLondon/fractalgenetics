#################
## libraries ####
#################

options(bitmapType = 'cairo', device = 'pdf')
modules::import_package('ggplot2', attach=TRUE)
optparse <- modules::import_package('optparse')

############
## data ####
############

## read command line args ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
                        dest="directory", type="character",
                        help="Path to directory with imputation QC results
                        [default: %default].", default=NULL),
    optparse$make_option(c("-n", "--name"), action="store",
                        dest="name", type="character",
                        help="Name of dataset [default: %default].",
                        default=NULL))
args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))
directory <- args$directory
alg <- args$alg
batch1 <- args$batch1
batch2 <- args$batch2
batch3 <- args$batch3
ukb <- args$ukb


directory <- "~/data/genotype"
alg <- "gencall"
batch1 <- "sanger12"
batch2 <- "singapore12"
batch3 <- "singapore3"
suffix <- "clean.related"
bgen <- "~/data/genotype/imputation/combined/genotypes/gencall.combined.clean.related.chr1.sample"
HVOL <- "~/data/genotype/QC/combined/European.HVOL.gencall.combined.txt"

data1 <- data.table::fread(paste(directory, "/QC/", batch1, "/", alg, ".",
                                 batch1, ".fam", sep=""),
                           header=FALSE, stringsAsFactors=FALSE,
                           data.table=FALSE)

data2 <- data.table::fread(paste(directory, "/QC/", batch2, "/", alg, ".",
                                 batch2,".fam", sep=""),
                           header=FALSE, stringsAsFactors=FALSE,
                           data.table=FALSE)

data3 <- data.table::fread(paste(directory, "/QC/", batch3, "/", alg, ".",
                                 batch3, ".fam", sep=""),
                           header=FALSE, stringsAsFactors=FALSE,
                           data.table=FALSE)

combined <- data.table::fread(paste(directory, "/QC/combined/", alg,
                                    ".combined.", suffix, ".fam", sep=""),
                           header=FALSE, stringsAsFactors=FALSE,
                           data.table=FALSE)

allsamples <- data.table::fread(bgen, stringsAsFactors=FALSE, data.table=FALSE)
allsamples <- allsamples[-1,]

HVOLsamples <- data.table::fread(HVOL, stringsAsFactors=FALSE, data.table=FALSE)
name <- paste("HVOL", ".", alg, ".combined.", suffix, sep="")

################
## analysis ####
################

## create and format data for BGENIE GWAS against batch ####
colnames(combined) <- c("FID", "IID", "PID","MID","SEX", "PHENO")
combined$PHENO <- 1
combined$PHENO[combined[,2] %in% data2[,2]] <- 2
combined$PHENO[combined[,2] %in% data3[,2]] <- 3

write.table(combined,
            paste(directory, "/control/ALL.", alg, ".combined.", suffix,
                  ".batchPhenotypes.fam", sep=""),
            sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)

combined$PHENO[!(combined[,1] %in% HVOLsamples[,2])] <- "-9"
write.table(combined,
            paste(directory, "/control/HVOL.", alg, ".combined.", suffix,
                  ".batchPhenotypes.fam", sep=""),
            sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE)

bgenie <- merge(allsamples, combined, by=1, all.x=TRUE, sort=FALSE)
bgenie <- bgenie[match(allsamples[,1], bgenie[,1]),]
bgenie$PHENO[!(bgenie[,1] %in% HVOLsamples[,2])] <- "-999"


write.table(bgenie$PHENO,
            paste(directory, "/control/Batch_phenotypes_", name, "_bgenie.txt",
                  sep=""),
            sep=" ", row.names=FALSE, col.names="batch", quote=FALSE)


# 1. Run PCA on directly genotyped SNPs from entire HVOL cohort
system(paste("plink --bfile ",  directory, "/gencall_all", " --pca --out ", dir, "/gencall_all.batch.pruned.pca", sep=""))

# 2. Read PCA results (eigenvectors)
pca_data <- read.table(paste(directory,"/gencall_all.batch.pruned.pca.eigenvec", sep=""), header=FALSE)
colnames(pca_data) <- c("FAMID", "IID", paste("PC", seq(1, (ncol(pca_data)-2), 1), sep=""))

# 3. Assign samples to corresponding sample subset
pca_data$type <- "subset"
pca_data$type[which(pca_data[,1] %in% samples_HVOL$omnix)] <- "allHVOL"
pca_data$type[which(pca_data[,1] %in% samples_original$omnix)] <- "original"
pca_data$type[which(pca_data[,1] %in% samples_subset$omnix)] <- "subset"
pca_data$type <- as.factor(pca_data$type)

#£ 4. Depict first three PCs as scatterplots
pdf(file=paste("~/GWAS/data/genotype/MRI_genotype/QC/",  strftime(Sys.time(), "%Y%m%d"),"_HVOL_original_batcheffects.pdf", sep=""), onefile=TRUE, height=12,width=16, paper = "a4r")
p <- ggplot(data=pca_data, aes(x=PC1, y=PC2))
p + geom_point(aes(color=type))
p <- ggplot(data=pca_data, aes(x=PC1, y=PC2))
p + geom_point(aes(color=type)) +
    ylim(c(-0.1,0.25)) +
    xlim(c(-0.1,0.25))

p <- ggplot(data=pca_data, aes(x=PC1, y=PC3))
p + geom_point(aes(color=type))
p <- ggplot(data=pca_data, aes(x=PC1, y=PC3))
p + geom_point(aes(color=type)) +
    ylim(c(-0.1,0.25)) +
    xlim(c(-0.1,0.25))

p <- ggplot(data=pca_data, aes(x=PC3, y=PC2))
p + geom_point(aes(color=type))
p <- ggplot(data=pca_data, aes(x=PC3, y=PC2))
p + geom_point(aes(color=type)) +
    ylim(c(-0.1,0.25)) +
    xlim(c(-0.1,0.25))
dev.off()


# outlier samples 
outlier_samples <- filter(pca_data, PC1 >0.25 | PC2 >0.25| PC3 >0.25)
write.table(outlier_samples, paste(directory, "/HVOL_pca_outlier_samples.txt", sep=""), sep="\t", col.names=TRUE, row.names=FALSE)




#######################################################################################
### overview of any ID matching/manipulation done with data received from London    ###
###                                                                                 ###
###         * fix typos in sanger manifest file                                     ###
###         * generate overview of genotype IDs                                     ###
###             (raw, preQC -> manual removal of IDs)                               ###
###         * generate overview of all BRU IDs                                      ###
###             (HVOL, non-HVOL, family membership)                                 ###
###         * write keep_lists (for plink) and exclude lists (for impute/snptest)   ###
###             (per diganosis across all batches)                                  ###
###                                                                                 ###
###         by Hannah Meyer                                                         ###
###                                                                                 ###
#######################################################################################



#################
### functions ###
#################

dob2age <- function(dob) {
             year = as.numeric(gsub(".*/", "", dob))
             month = gsub("\\d{1,2}/(\\d{1,2})/\\d{4}", "\\1", dob)
             month = as.numeric(gsub("^0","", month))
             if (month < 4 ) {
                return(2016-year)
             } else {
                 return(2016-year-1)
             }}

IDconvert <-  function(id, convert, to=2, from=1) {
            if (any(convert[,from] == id))  {
                return(convert[which(convert[,from] == id), to])
            } else {
                return (id)
            }
}

#####################################
### 1. Sanger sequencing id file: ###
#####################################

phenodir <- "/homes/hannah/data/digital-heart/phenotype/2Dphenotype"
qcdir <- "/homes/hannah/data/digital-heart/genotype/QC"
sangerdir <- "/homes/hannah/data/digital-heart/genotype/omnix_hhtmri_20141121"
singaporedir <- "/homes/hannah/data/digital-heart/genotype/omnix_NHCS_20160217"

# maps BRU numbers to omnix identifiers
### a)  overview table contains samples that were excluded from genotyping
### -> mapping to BRU IDs will clash with later batches as samples were sent for
### sequencing to singapore with BRU ID
# i) save copy of original file
system(paste("cp ", sangerdir, "/gencall_qc/pipeline_summary.csv ", sangerdir,
             "/gencall_qc/pipeline_summary.originalSanger.BRUtypos.csv",
             sep=""))
# ii) remove samples from overview table that were not genoytyped
system(paste("head -n 1 ", sangerdir,
             "gencall_qc/pipeline_summary.originalSanger.BRUtypos.csv > ",
             sangerdir, "gencall_qc/pipeline_summary_genotyped.csv;",
             "grep -f <(cut -d ' ' -f 1 ", sangerdir,
             "/omnix_hhtmri_20141121.gencall.smajor.fam) ", sangerdir,
             "/gencall_qc/pipeline_summary.originalSanger.BRUtypos.csv >> ",
             sangerdir, "/gencall_qc/pipeline_summary_genotyped.csv", sep=""))

# iii) two of the original BRU numbers contained typos (rt Rachel's email April
# 13, 2016), so will be manually changed in the manifest file:
# were before:  14KU01740, 14JG01942 -> should be: 14KV01740, 14JS01942
manifest <- read.table(paste(sangerdir,
                             "/gencall_qc/pipeline_summary_genotyped.csv",
                             sep=""),
                       sep=",", stringsAsFactors=FALSE, header=TRUE)

manifest[manifest$supplier_name == "14KU01740", "supplier_name"] <- "14KV01740"
manifest[manifest$supplier_name == "14JG01942", "supplier_name"] <- "14JS01942"

write.table(manifest, paste(sangerdir, "/gencall_qc/pipeline_summary.csv",
                            sep=""),
            sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)

########################
### 2. BRU database: ###
########################


### a) read conversion overview file (generated in genotypeQC.R)
id_convert <- read.table(paste(qcdir,
                               "/gencall.sanger12/gencall_sampleID.txt",
                               sep=""),
                         sep="\t", header=TRUE,stringsAsFactors=FALSE)

### b) split into HVOL (Email Paz: 09.02.2016) and non-HVOL samples
# (Email Paz: 09.02.2016)
# HVOL data
bru_HVOL <- read.table(paste(phenodir, "/20160412_HVOL_BRU.csv", sep=""),
                       sep=",", header=TRUE, stringsAsFactors=FALSE)[,c(1:4,7:9)]
colnames(bru_HVOL) <- c("Bru.Number", "Ethnicity", "Sex", "ConfirmedDiagnosis",
                        "Age", "Height", "Weight")
# non-HVOL data
bru_nonHVOL <- read.table(paste(phenodir, "/All_BRU_Hannah_9.2.16.csv", sep=""),
                          sep=",", header=TRUE, stringsAsFactors=FALSE)[,1:8]
bru_nonHVOL$Age <- sapply(bru_nonHVOL$Dob.Str, dob2age)
bru_nonHVOL <- bru_nonHVOL[,c(1:4,9,7:8)]
colnames(bru_nonHVOL) <- c("Bru.Number", "Ethnicity", "Sex",
                           "ConfirmedDiagnosis", "Age", "Height", "Weight")

# i) combine both databases
bru <- rbind(bru_nonHVOL, bru_HVOL)

# remove duplicated samples (check before that measurements are the same
# (in 3 cases weight/height of by +1)
bru <- bru[!duplicated(bru$Bru.Number),]
write.table(bru, paste(phenodir, "/20160412_All_BRU_format.txt", sep=""),
            sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

# ii) format family-id file
# family information: match fam info (all BRU ids) to geno ids from first
# batch (omnix/sanger ids)
fam_matching <- read.table(paste(phenodir,
                                 "/20160209_All_BRU_Hannah_FamilyNo.txt",
                                 sep=""), sep="\t", header=TRUE,
                           stringsAsFactors=FALSE)
fam_matching$Geno.Number <- sapply(fam_matching$Bru.Number, IDconvert,
                                   convert=id_convert, from=2, to=1)
fam_matching <- fam_matching[, c(6,1,5,2,3,4)]
colnames(fam_matching) <- c("Geno.Number", "Bru.Number", "Fam.Number",
                            "Ethnicity", "Sex", "ConfirmedDiagnosis")

# Email Katie April 18: add family id to one individual
fam_matching[which(fam_matching$Geno.number == "hhtmri5420810"),
             "Fam.Number"] <- 272
write.table(fam_matching, paste(phenodir, "/20160209_All_BRU_family_format.txt",
                                sep=""), sep="\t",col.names=TRUE,
            row.names=FALSE, quote=FALSE)

### c)  Create overview of all BRU samples based on diagnosis
# (independent of genotype status)
samples_by_diagnosis <- split(fam_matching,
                              f=as.factor(fam_matching$ConfirmedDiagnosis))
samples_by_diagnosis_overview <-
    data.frame(diagnosis=names(samples_by_diagnosis),
               nr_samples = sapply(samples_by_diagnosis, nrow))
write.table(samples_by_diagnosis_overview,
            paste(phenodir, "/All_BRU_by_diagnosis.csv", sep=""),
            quote=FALSE, row.names=FALSE, col.names=TRUE)

################################################################
### 3. generate overview of genotypes from different batches ###
################################################################

### a) read .fam files of preQC genotype files
### (IDs processed -renamed/regex matched- via formatRawGenotypes.sh)
sanger1.raw <- read.table(paste(sangerdir, "/gencall.tmp.fam", sep=""),
                          stringsAsFactors=FALSE, header=FALSE)
colnames(sanger1.raw) <- c("FamID", "IndID", "FatherID", "MotherID", "Sex",
                           "Pheno")
sanger1.raw$FamID <- sapply(sanger1.raw$FamID, IDconvert, convert=id_convert)

sanger12.raw <- read.table(paste(sangerdir, "/gencall.sanger12.raw.fam", sep=""),
                               stringsAsFactors=FALSE, header=FALSE)
colnames(sanger12.raw) <- c("FamID", "IndID", "FatherID", "MotherID", "Sex",
                                "Pheno")
sanger12.raw$FamID <- sapply(sanger12.raw$FamID, IDconvert, convert=id_convert)

singapore12.raw=read.table(paste(singaporedir, "/gencall.singapore12.raw.fam",
                                 sep=""), stringsAsFactors=FALSE, header=FALSE)
colnames(singapore12.raw) <- c("FamID", "IndID", "FatherID", "MotherID", "Sex",
                               "Pheno")

singapore3.raw=read.table(paste(singaporedir, "/gencall.singapore3.raw.fam",
                                sep=""), stringsAsFactors=FALSE, header=FALSE)
colnames(singapore3.raw) <- c("FamID", "IndID", "FatherID", "MotherID", "Sex",
                              "Pheno")

sanger1.raw$batch <- "sanger1"
sanger12.raw$batch <- "sanger12"
singapore12.raw$batch <- "singapore12"
singapore3.raw$batch <- "singapore3"

# sanger1 part of gencall.sanger12, so don't includ into the combined file
all_genotypes_raw <- rbind(sanger12.raw, singapore12.raw, singapore3.raw)
write.table(all_genotypes_raw,
            paste(qcdir, "/genotype_overview_sanger12_singapore123_raw.txt",
                  sep=""), sep="\t", quote=FALSE)

# get overview of samples with regards to diagnosis in the different batches
all_genotypes_bru <- merge(all_genotypes_raw, bru, by=1)
sanger1.bru <- merge(sanger1.raw, bru, by=1)

overview <- rbind(all_genotypes_bru[,c(7, 10)], sanger1.bru[,c(7,10)])
overview_bySex <- rbind(all_genotypes_bru[,c(7, 10,9)], sanger1.bru[,c(7,10,9)])
write.table(table(overview),
            paste(qcdir,
                  "/genotype_overview_sanger12_singapore123_diagnosis.txt",
                  sep=""),
            sep="\t", quote=FALSE)

### b)  remove geno IDs with 'questionable' status
### (decision to remove after discussion with London)
# Email April 15 by Paz: "10JP00134","10DK01008"->possibe sample mix up, PIHAT= 1
# both in singapore 3
remove <- c("10JP00134", "10DK01008")
# Email April 13 by Paz : "14AH01709", "14SW01745", "14RK03047", "14HM03269"
# -> no confirmed diagnosis in (sanger12, sanger12, singapore12,
remove <- c(remove, "14AH01709", "14SW01745", "14RK03047", "14HM03269")
# Email Rachel April 11: "10TW01294" -> genotyped twice, better genotyping rate
# in singapore3 (determined after QC)
remove <- c(remove, "10TW01294")
# Email March 8 by Paz:"10JH04320" -> possible DCM but no confirmed diagnosis
remove <- c(remove, "10JH04320")
# Email March 7 by Paz: "10NT00915"-> Not in database (singapore12),
# "14RK03047" -> no confirmed diagnosis  (singapore12),
# "10NJ00354" -> duplicate ID (singapore12)
remove <- c(remove, "10NT00915", "14RK03047", "10NJ00354")

remove <- unique(remove)
geno_remove <- sapply(remove, IDconvert, id_convert)
geno_remove_batch <- merge(geno_remove, all_genotypes_raw, by=1)

# Email Rachel April 11: "10TW01294" -> genotyped twice, better genotyping rate
# in singapore3 (determined after QC), remove once from list
geno_remove_batch <-
    geno_remove_batch[-intersect(which(geno_remove_batch$IndID=="10TW01294"),
                                 which(geno_remove_batch$batch=="singapore12")),]
geno_remove_batch <- geno_remove_batch[! duplicated(geno_remove_batch$IndID),]

### c)  write overview of genotypes after removal of preQC fails
all_genotypes_preQC <-
    all_genotypes_raw[!(all_genotypes_raw[,1] %in%
                        geno_remove_batch[geno_remove_batch[,1]!="10TW01294",1]),]
all_genotypes_preQC <-
    all_genotypes_preQC[-intersect(which(all_genotypes_preQC$IndID=="10TW01294"),
                                   which(all_genotypes_preQC$batch=="singapore12")),]
write.table(all_genotypes_preQC,
            paste(qcdir, "/genotype_overview_sanger12_singapore123_preQC.txt",
                  sep=""), sep="\t", quote=FALSE)

### d)  split genotypes to be removed by batch again and write to preQCfail.IDs
### in respective directory (duplicates have to be removed before using plink)
geno_remove_batch <- split(geno_remove_batch,
                           f=as.factor(geno_remove_batch$batch))
write.batch <-
    sapply(seq_along(geno_remove_batch), function(batch) {
               file <- paste(qcdir, "/gencall.", names(geno_remove_batch)[batch],
                             "/gencall.", names(geno_remove_batch)[batch],
                             ".preQCfail.IDs", sep="")
               write.table(geno_remove_batch[[batch]][,1:2],
                           file=file, col.names=FALSE, row.names=FALSE,
                           quote=FALSE, sep= "\t"))


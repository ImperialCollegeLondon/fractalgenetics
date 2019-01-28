###############################
### libraries and functions ###
###############################

library("crlmm")
library("illuminaio")
library("data.table")
library("plyr")


vmessage <- function(userinfo, display=TRUE, sep=" ") {
  if (display) message(paste(userinfo, collapse=sep))
}


generateEvokerFiles <- function(dataDir, dataname, dataset, fam, bim, mapData,
                                output, inputdata, SNPs2analyse, evokerDir,
                                thr=5e-3, sampleSheet=NULL, idatRed=NULL,
                                display=TRUE) {
    if (is.null(idatRed)) {
        vmessage(c("Read .idat files for samples in dataset and filter for ",
                   "SNPs in mapData and qc'ed .bim file"), display=display)
        # reduce sampleSheet, order samples according to .fam file order 
        # (required by evoker)
        sampleSheet <- sampleSheet[which(sampleSheet$Sample_ID %in% fam[,1]),]
        sampleSheet <- sampleSheet[match(fam[,1], sampleSheet$Sample_ID),]
        sampleSheet$fused <- paste(sampleSheet$SentrixBarcode_A,
                                   sampleSheet$SentrixPosition_A, sep="_")

        # read .idat files of both channels and create matrix of
        # nrSNP x nrsample*[sample_intensitesred + sample_intensitiss green]
        idatMatrix <- lapply(seq_along(sampleSheet$fused), function(x, dataDir) {
            red <- readIDAT(paste(dataDir,"/", sampleSheet$fused[x],
                                  "_Red.idat", sep=""))
            green <- readIDAT(paste(dataDir,"/", sampleSheet$fused[x],
                                    "_Grn.idat", sep=""))
            signals <- cbind(red$Quants[,1], green$Quants[,1])
            colnames(signals) <- rep(sampleSheet$Sample_ID[x],2)
            rownames(signals) <- red$MidBlock
            return(signals)}, dataDir=dataDir)
        idatMatrix <- do.call(cbind, idatMatrix)
        saveRDS(idatMatrix, file=paste(output, "/idatMatrix.rds", sep=""))

        # match the order of the intensity values to the order of the mapping
        # information and replace AdressA_ID with SNP ID
        idatRed <- idatMatrix[which(rownames(idatMatrix) %in%
                                    mapData$AddressA_ID),]
        mapDataRed <- mapData[which(mapData$AddressA_ID %in%
                                    rownames(idatMatrix)),]
        mapDataOrdered <- mapDataRed[match(as.numeric(rownames(idatRed)),
                                           mapDataRed$AddressA_ID),]
        rownames(idatRed) <- mapDataOrdered$Name

        # match the order of the intensity values to the order of the SNPs in
        # the automated qc'ed .bim file
        idatRed <- idatRed[which(rownames(idatRed) %in% bim[,2]),]
        bimRed <- bim[which(bim[,2] %in% rownames(idatRed)),]
        saveRDS(idatRed, file=paste(output, "/idatRed.rds", sep=""))
    } else {
        idatRed <- readRDS(file=paste(output, "/idatRed.rds", sep=""))
    }

    # extract intensity values of SNPs to be manually qc'ed
    vmessage(c("Filter for SNPs to be manually qc'ed"), display=display)
    idatForint2bnt <- idatRed[which(rownames(idatRed)  %in% SNPs2analyse$SNP),]
    idatForint2bnt <- merge(idatForint2bnt, by.x= 0, SNPs2analyse, by.y=1)
    colnames(idatForint2bnt)[1] <- c("SNP")

    system(paste("mkdir -p ", evokerDir, "/", dataset, "/Thr", thr,
                 "/prefiles", sep=""))
    write.table(idatForint2bnt[,1], paste(evokerDir, "/", dataset,
                                          "/prefiles/SNPForint2bnt.txt",sep=""),
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    system(paste("plink --noweb --bfile ", inputdata, " --extract ", evokerDir,
                 "/", dataset, "/prefiles/SNPForint2bnt.txt --make-bed --out ",
                 evokerDir, "/", dataset, "/Thr", thr, "/prefiles/", dataname,
                 ".evoker", sep="" ))

    # split per chromosome and order with increasing bp
    idatForint2bnt_list <- split(idatForint2bnt,
                                 f=as.factor(idatForint2bnt$CHR))
    idatForint2bnt_list <- lapply(idatForint2bnt_list, function(x)
                                  x[order(x$BP, decreasing=FALSE),])

    # write output files for evoker
    vmessage(c("Write evoker output files to", evokerDir), display=display)
    l_ply(seq(1, length(idatForint2bnt_list),1), function(x) {
        datacols <- ncol(idatForint2bnt_list[[x]])
        write.table(rbind(colnames(idatForint2bnt_list[[x]])[-c((datacols-2):datacols)],
                          as.matrix(idatForint2bnt_list[[x]][,-c((datacols-2):datacols)])),
                    file=paste(evokerDir, "/", dataset, "/Thr", thr, "/",
                               dataname, ".", x, ".int", sep=""), sep="\t",
                    quote=FALSE, col.names=FALSE, row.names=FALSE)
    })

    l_ply(seq(1, length(idatForint2bnt_list),1), function(x) {
        system.cmd  <- paste("plink --noweb --bfile ", evokerDir,  "/", dataset,
                             "/Thr", thr, "/prefiles/", dataname,
                             ".evoker --chr ", x, " --make-bed --out ",
                             evokerDir, "/", dataset, "/Thr", thr, "/",
                             dataname, ".", x, sep="")
        system(system.cmd)
    })

    l_ply(seq(1, length(idatForint2bnt_list),1), function(x) {
        system.cmd  <- paste("perl ~/GWAS/analysis/genotyping/int2bnt.pl  -i ",
                             evokerDir,  "/", dataset,  "/Thr", thr, "/",
                             dataname, ".", x,".int -f default -o ", evokerDir,
                             "/", dataset, "/Thr", thr, "/", dataname, ".", x,
                             ".bnt",  sep="")
        system(system.cmd)
    })
    return(list(idatRed=idatRed, idatForint2bnt=idatForint2bnt))
}



############
### data ###
############

### genotyping data and metadata
genodir <- '~/data/digital-heart/genotype'
evokerDir <- paste(genodir,"/QC/evoker", sep="")

thr=5e-3

### GWAS results
gwas_results_FD <- fread(paste("~/data/digital-heart/association/FD/",
                               "bgenie_slices_lm_pseudomt_all.csv", sep=""),
                         data.table=FALSE, stringsAsFactors=FALSE)
gwas_results_3D <- fread(paste("~/data/digital-heart/association/",
                               "HVOL/ManifoldLearning/Combined/",
                               "lm_mt_pcs_pvalue_genomeExX.csv", sep=""),
                               data.table=FALSE, stringsAsFactors=FALSE)

gwas_results_FD_thr <- gwas_results_FD[which(gwas_results_FD$P < thr),]
gwas_results_3D_thr <- gwas_results_3D[which(gwas_results_3D$P < thr),]
gwas_results <- rbind(gwas_results_FD_thr, gwas_results_3D_thr)

evokerSigSNP <- gwas_results[!duplicated(gwas_results$SNP), c(3,1,2,4)]
evokerSigSNP <- evokerSigSNP[order(evokerSigSNP[,2], evokerSigSNP[,3]),]
write.table(evokerSigSNP,
            paste(evokerDir, "/all/prefiles/allSNPs2analyse.txt", sep=""),
            col.names=FALSE, row.names=FALSE, quote=FALSE)

# 1. Sanger batch
dataDirSanger <- paste(genodir,
                       "/omnix_hhtmri_20141121/omnix_hhtmri_20141121.idat",
                       sep="")

sampleSheetSanger <-
    read.table(paste(dataDirSanger, "samplesheet.csv", sep=""),
               stringsAsFactors=FALSE)
names(sampleSheetSanger) <- c("Full_ID", "SentrixBarcode_A", "SentrixPosition_A",
                              "Sample_Plate", "Sample_Well", "gtc", "grn.idat",
                              "red.idat")
sampleSheetSanger$Sample_ID <- gsub("[0-9]*_[A-H0-9]*_","",
                                    sampleSheetSanger$Full_ID)

famSanger <-  read.table(paste(genodir,
                                "/QC/sanger12/gencall.sanger12.clean.related.fam",
                                sep="")
bimSanger <-  read.table(paste(genodir,
                                "/QC/sanger12/gencall.sanger12.clean.related.bim",
                                sep="")
mapDataSanger <- read.table(paste(genodir, "/omnix_hhtmri_20141121/",
                                  "humanomniexpress-12v1-1_a.csv",
                                  sep=""), skip=7, sep=",",
                            stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
mapDataSanger <- mapDataSanger[-which(duplicated(mapDataSanger$AddressA_ID)),]



# 2. Singapore
dataDirSingapore <- paste(genodir, "/omnix_NHCS_20160217/samplesheet_all",
                          sep="")

sampleSheetSingapore_HVOL_Plate01_03 <-
    read.table(paste(dataDirSingapore,
                     "/HVOL_Plate01-03_Omniexp_285samples.csv", sep=""),
                     stringsAsFactors=FALSE, header=TRUE, sep=",", skip=9)
sampleSheetSingapore_DCM_Plate01_03 <-
    read.table(paste(dataDirSingapore,
                     "/DCM_Plate01-03_Omniexp_246samples.csv", sep=""),
               stringsAsFactors=FALSE, header=TRUE, sep=",", skip=9)
sampleSheetSingapore_NAD_Plate01_03 <-
    read.table(paste(dataDirSingapore,
               "/NAD_Plate01-03_Omniexp_240samples.csv", sep=""),
               stringsAsFactors=FALSE, header=TRUE, sep=",", skip=9)
sampleSheetSingapore_HCM_Plate03_05 <-
    read.table(paste(dataDirSingapore,
               ,"/HCM_Plate03-05_Omniexp_258samples.csv", sep=""),
               stringsAsFactors=FALSE, header=TRUE, sep=",", skip=9)
sampleSheetSingapore_HCM_192samples_Batch01 <-
    read.table(paste(dataDirSingapore,
               "/HCM_192samples_batch01.csv", sep=""),
               stringsAsFactors=FALSE, header=TRUE, sep=",", skip=9)
sampleSheetSingapore_Heart_Failure_432samples <-
    read.table(paste(dataDirSingapore,
               "/Edmund_NHCS_Heart_Failure_432samples.csv", sep=""),
               stringsAsFactors=FALSE, header=TRUE, sep=",", skip=9)

sampleSheetsAll <- rbind(sampleSheetSingapore_HVOL_Plate01_03,
                         sampleSheetSingapore_DCM_Plate01_03,
                         sampleSheetSingapore_NAD_Plate01_03,
                         sampleSheetSingapore_HCM_Plate03_05,
                         sampleSheetSingapore_HCM_192samples_Batch01,
                         sampleSheetSingapore_Heart_Failure_432samples)
sampleSheetsAll$Full_ID <- sampleSheetsAll$Sample_ID
sampleSheetsAll$grn.idat <- paste(sampleSheetsAll$SentrixBarcode_A,
                                  sampleSheetsAll$SentrixPosition_A,
                                  "Grn.idat", sep="_")
sampleSheetsAll$red.idat <- paste(sampleSheetsAll$SentrixBarcode_A,
                                  sampleSheetsAll$SentrixPosition_A,
                                  "Red.idat", sep="_")

# map correct ID's (wrongly labeled in at genotyping facility)
IDmap <- read.table(paste(genodir,
                          "/omnix_NHCS_20160217/corrected_sample_id.txt",
                          sep=""), header=FALSE, stringsAsFactors=FALSE)
sampleSheetsAll$correctID <- sapply(sampleSheetsAll$Sample_ID, function(x) {
    if (x %in% IDmap[,2] ) return(IDmap[x == IDmap[,2],4])
    if (! x %in% IDmap[,2]) return(x)
})

sampleSheetsAll$Sample_ID <- gsub("[0-9]*_([A-Z0-9]*)_.*","\\1",
                                  sampleSheetsAll$correctID)

## HumanOmniExpress-24v1-0
famSingapore12 <- fread(paste(genodir, "/QC/gencall.singapore12/",
                              "gencall.singapore12.clean.related.fam", sep=""),
                        data.table=FALSE)
bimSingapore12 <- fread(file=paste(genodir, "/QC/gencall.singapore12/",
                                   "gencall.singapore12.clean.related.bim",
                                   sep=""), data.table=FALSE)
mapDataSingapore12 <-
    read.table(paste(genodir, "/omnix_NHCS_20160217/",
                     "HumanOmniExpress_24_sample_manifests/",
                     "HumanOmniExpress-24v1-0/humanomniexpress-24v1-0_a.csv",
                     sep=""),
               skip=7, sep=",", stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
mapDataSingapore12 <-
    mapDataSingapore12[-which(duplicated(mapDataSingapore12$AddressA_ID)),]

## HumanOmniExpress-24v1-0
famSingapore3 <- fread(paste(genodir, "/QC/gencall.singapore3/",
                             "gencall.singapore3.clean.related.fam", sep=""),
                       data.table=FALSE)
bimSingapore3 <- fread(paste(genodir,"/QC/gencall.singapore3/",
                             "gencall.singapore3.clean.related.bim", sep=""),
                       data.table=FALSE)
mapDataSingapore3 <-
    read.table(paste(genodir,"/omnix_NHCS_20160217/",
                     "HumanOmniExpress_24_sample_manifests/",
                     "HumanOmniExpress-24v1-1_A_NewVersion2015/",
                     "humanomniexpress-24-v1-1-a.csv", sep=""),
               skip=7, sep=",", stringsAsFactors=FALSE, fill=TRUE, header=TRUE)
mapDataSingapore3 <-
    mapDataSingapore3[-which(duplicated(mapDataSingapore3$AddressA_ID)),]

################
### analysis ###
################

sangerevoker <- generateEvokerFiles(dataset="gencall.sanger12",
  dataname="sanger12",
  fam=famSanger,
  bim=bimSanger,
  mapData=mapDataSanger,
  output=paste(genodir, "/QC/gencall.sanger12", sep=""),
  inputdata=paste(genodir, "/QC/gencall.sanger12/gencall.sanger12.clean.related",
                  sep=""),
  sampleSheet=sampleSheetSanger,
  dataDir=dataDirSanger,
  SNPs2analyse=evokerSigSNP,
  evokerDir=evokerDir)

singapore12evoker <- generateEvokerFiles(dataset="gencall.singapore12",
  dataname="singapore12",
  fam=famSingapore12,
  bim=bimSingapore12,
  mapData=mapDataSingapore12,
  output=paste(genodir, "/QC/gencall.singapore12", sep=""),
  inputdata=paste(genodir, "/QC/gencall.singapore12/",
                  "gencall.singapore12.clean.related", sep=""),
  sampleSheet=sampleSheetsAll,
  dataDir=dataDirSingapore,
  SNPs2analyse=evokerSigSNP,
  evokerDir=evokerDir)

singapore3evoker <- generateEvokerFiles(dataset="gencall.singapore3",
  dataname="singapore3",
  fam=famSingapore3,
  bim=bimSingapore3,
  mapData=mapDataSingapore3,
  output=paste(genodir, "/QC/gencall.singapore3", sep=""),
  inputdata=paste(genodir, "/QC/gencall.singapore3/",
                  "gencall.singapore3.clean.related", sep=""),
  sampleSheet=sampleSheetsAll,
  dataDir=dataDirSingapore,
  SNPs2analyse=evokerSigSNP,
  evokerDir=evokerDir)

genotypedEvokerSigSNPs <- unique(c(sangerevoker$idatForint2bnt[,1],
                                   singapore12evoker$idatForint2bnt[,1],
                                   singapore3evoker$idatForint2bnt[,1]))
write.table(genotypedEvokerSigSNPs,
            paste(evokerDir, "/all/prefiles/genotypedSNPs2analyseThr",
                  thr, ".txt", sep=""), col.names=FALSE, row.names=FALSE,
            quote=FALSE)

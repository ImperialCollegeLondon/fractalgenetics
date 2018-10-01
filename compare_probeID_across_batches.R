###################################################
###################################################
###                                             ###
###     Analyse different genotyping platforms  ###
###     * number of shared and unique SNPs      ###
###     * for concordance of SNP probes         ###
###                                             ###
###     by Hannah Meyer                         ###
###                                             ###
###################################################
###################################################


#################
### analysis  ###
#################

directory <- "/homes/hannah/GWAS/data/genotype"

### 1. SNP info data for the different platforms/batches
# Sanger batch 1 and 2
sanger <- read.table(paste(directory,
                           "omnix_hhtmri_20141121/HumanOmniExpress-12v1-1_A.csv",
                           sep=""), sep=",", header=TRUE, skip=7, fill=TRUE,
                     stringsAsFactors=FALSE)
sanger <- sanger[!is.na(sanger$AddressA_ID),]
sanger_chip <- "HumanOmniExpress-12v1-1_A"

# Singapore combined batch 1 and 2, and batch
singapore12 <- read.table(paste(directory,
                                "omnix_NHCS_20160217/Manifests/",
                                "HumanOmniExpress-24v1-0/",
                                "humanomniexpress-24v1-0_a.csv", sep=""),
                          sep=",", header=TRUE, skip=7, fill=TRUE,
                          stringsAsFactors=FALSE)
singapore12 <- singapore12[!is.na(singapore12$AddressA_ID),]
singapore12_chip <- "HumanOmniExpress-24v1-0"

singapore3 <- read.table(paste(directory, "/omnix_NHCS_20160217/Manifests/",
                               "HumanOmniExpress-24v1-1_A/",
                               "humanomniexpress-24-v1-1-a.csv", sep=""),
                         sep=",", header=TRUE, skip=7, fill=TRUE,
                         stringsAsFactors=FALSE)
singapore3 <- singapore3[!is.na(singapore3$AddressA_ID),]
singapore3_chip <-"HumanOmniExpress-24v1-1_A"

ID_overlap <- VennDiagram::calculate.overlap(
                                    list("Duke-NUS12" =singapore12[,"Name"],
                                         "Duke-NUS3"=singapore3[,"Name"],
                                          Sanger12=sanger[,"Name"]))
venn.plot <- VennDiagram::venn.diagram(list("Duke-NUS12" =singapore12[,"Name"],
                               "Duke-NUS3"=singapore3[,"Name"],
                               "Sanger12" =sanger[,"Name"]),
                          filename=paste(directory,
                                         "/Venn_genotyping_batches.pdf", sep=""),
                          imagetype="pdf",
                          category.names=
                              c("Duke-NUS12\nHumanOmniExpress-24v1-0",
                                "Duke-NUS3\nHumanOmniExpress-24v1-1_A",
                                "Sanger12\nHumanOmniExpress-12v1-1_A"),
                          cat.pos = c(-20, 0, 25),
                          cat.dist = c(0.08, 0.08, 0.02),
                          cat.col = c("darkorchid2","darkblue", "darkgreen"),
                          cat.cex=0.8
                          )

singapore <- merge(singapore12, singapore3, by="Name")

# merge based on rsIDs
all <- merge(sanger, singapore, by="Name")[, c("Name", "AlleleA_ProbeSeq",
                                               "AlleleA_ProbeSeq.x",
                                               "AlleleA_ProbeSeq.y")]

# all SNPs with same sequence ID
all_seq_index <- apply(all, 1, function(x) {
                           ifelse (length(unique(x[-1])) == 1, TRUE, FALSE)
                          })
all_seq <- all[all_seq_index,]

### -> SNPs common to all platforms have same probes, check!




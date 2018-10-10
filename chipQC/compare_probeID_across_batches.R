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


#############
### data  ###
#############

directory <- "~/data/genotype"

## 1. SNP info data for the different platforms/batches ####
# Sanger batch 1 and 2
sanger <- read.table(paste(directory,
                           "/omnix_hhtmri_20141121/HumanOmniExpress-12v1-1_A.csv",
                           sep=""), sep=",", header=TRUE, skip=7, fill=TRUE,
                     stringsAsFactors=FALSE)
sanger <- sanger[!is.na(sanger$AddressA_ID),]
sanger_chip <- "HumanOmniExpress-12v1-1_A"

# Singapore combined batch 1 and 2, and batch 3
singapore12 <- read.table(paste(directory,
                                "/omnix_NHCS_20160217/Manifests/",
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

#################
### analysis  ###
#################

ID_list <- list("Duke-NUS12"=singapore12$Name, "Duke-NUS3"=singapore3$Name,
                                          "Sanger12"=sanger$Name)

pdf(paste(directory, "/overlap_genotypeProbes.pdf", sep=""))
UpSetR::upset(UpSetR::fromList(ID_list),
              order.by = "freq",
              empty.intersections = "on", text.scale=1.2,
              #scale.sets="log10", scale.intersections="log10",
              # Include when UpSetR v1.4.1 is released
              # title="Overview quality control failures",
              main.bar.color="#1b9e77", matrix.color="#1b9e77",
              sets.x.label="Number of probes",
              mainbar.y.label="Number of common probes",
              sets.bar.color="#d95f02")
dev.off()

singapore <- merge(singapore12, singapore3, by="Name")

# merge based on rsIDs
all <- dplyr::select_(merge(sanger, singapore, by="Name"),
                      ~Name, ~AlleleA_ProbeSeq, ~AlleleA_ProbeSeq.x,
                      ~AlleleA_ProbeSeq.y)

# all SNPs with same sequence ID
all_seq_index <- apply(all, 1, function(x) {
                           ifelse (length(unique(x[-1])) == 1, TRUE, FALSE)
                          })
all_seq <- all[all_seq_index,]

### -> SNPs common to all platforms have same probes, check!
nrow(all_seq) == nrow(all)
# [1] TRUE



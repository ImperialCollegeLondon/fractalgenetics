#################
### libraries ###
#################

library("ggplot2")
library ("grid")
library("gridExtra")

#################
## functions ####
#################

################
## analysis ####
################
#args <- commandArgs(asValue=TRUE)
#inputdir = args$inputdir
#outputdir = args$outputdir

inputdir = "/Users/hannah/Documents/GWAS/data/genotype/MRI_genotype/imputed/gencall.combined/plots"
outputdir = "/Users/hannah/Documents/GWAS/data/genotype/MRI_genotype/imputed/gencall.combined/plots"
centromere <- read.table(file="/Users/hannah/Documents/GWAS/data/genotype/MRI_genotype/imputed/20142903_GRCh37_centromere", header=TRUE)
centromere <- centromere[c(1:9,12:24,10),]
chr_begin <-  as.matrix(read.table(file="/Users/hannah/Documents/GWAS/data/genotype/MRI_genotype/imputed/chr_begin.txt", header=FALSE))


# Display the maximum posterior probability of imputation per chunk per chromosome
data_files = dir(inputdir, pattern = "_data", full.names = TRUE, ignore.case = TRUE)
number_files = dir(inputdir, pattern = "_numbers", full.names = TRUE, ignore.case = TRUE)

for (i in 1:length(number_files)) {
    if (i == 23) {
      chr="X"
      chr_txt="X"
    } else if (i == 24) {
      chr="X PAR1"
      chr_txt="X_PAR1"
      X_PAR=c(60001,2699520)
    } else if (i == 25) {
      chr="X PAR2"
      chr_txt="X_PAR2"
      X_PAR=c(154931044,999999999)
    } else if (i < 10 ) {
      chr=i
      chr_txt=paste("0", i, sep="")      
    } else {
      chr=i
      chr_txt=i
    }
    file <- paste(outputdir, "/chr", chr_txt, sep ="")
    
    numbers <- read.table(number_files[i], header=FALSE, sep="\t", stringsAsFactors=FALSE)
    data_matrix <- read.table(data_files[i], header=FALSE, sep="\t", stringsAsFactors=FALSE)
  
   # chunks_wo_snp <- intersect(which(numbers[,4]==0), which(numbers[,3] == 0))
   chunks_wo_snp <- sort(c(which(numbers[,5]==0), which(is.na(numbers[,5]))))
   if (any(chunks_wo_snp)) {
         for (b in seq_along(chunks_wo_snp)) {
            if (length(which(data_matrix[((chunks_wo_snp[b]-1)*10+1):((chunks_wo_snp[b])*10),3]==0)) == 10) {
              numbers[((chunks_wo_snp[b])),5] <- 0
              next
            }
            if(chunks_wo_snp[b] != 1 && chunks_wo_snp[b] != nrow(numbers)) {
              data_matrix <- rbind(data_matrix[1:((chunks_wo_snp[b]-1)*10),], matrix(nr=10, nc=4, 0),data_matrix[(((chunks_wo_snp[b]-1)*10)+1): nrow(data_matrix),] )
            } else {
               if(chunks_wo_snp[b] == 1) {
                  data_matrix <- rbind(matrix(nr=10, nc=4, 0),data_matrix)
               }
               if(chunks_wo_snp[b] == nrow(numbers)) {
                  data_matrix <- rbind(data_matrix,matrix(nr=10, nc=4, 0))
               }
            }
        }
        data_list <-split(as.data.frame(data_matrix), gl(nrow(data_matrix)/10,10))
        data_list_red <- data_list[-chunks_wo_snp]
        chunks_wo_snp_str <- paste(chunks_wo_snp,collapse=", ")
        numbers_red <- numbers[-chunks_wo_snp,] 
    } else {
        chunks_wo_snp_str <- "None"
        data_list <-split(as.data.frame(data_matrix), gl(nrow(data_matrix)/10,10))
        numbers_red <- numbers
        data_list_red <- data_list
    }
    data_reordered <- as.matrix(do.call(cbind, data_list))
    if (ncol(data_reordered) != 4) {
        genotypes <- data_reordered[, seq(3,(ncol(data_reordered)-1),4)]
        total <- colSums(data_reordered[, seq(3,(ncol(data_reordered)-1),4)])
        concordance <- data_reordered[, seq(4,ncol(data_reordered),4)]
        genotypes_upper <- genotypes[8:10,]
        concordance_upper <- concordance[8:10,]
    }
    if (ncol(data_reordered) == 4) {
        genotypes <- data_reordered[, 3]
        total <- sum(data_reordered[, 3])
        concordance <- data_reordered[, 4]
        genotypes_upper <- genotypes[8:10]
        concordance_upper <- concordance[8:10]
    }
    chunk <- rep(seq(1, length(concordance)/10,1), each= 10)
    chunk_upper <- rep(seq(1, length(concordance)/10,1), each= 3)
    MPP <- as.factor(rep(c( "[0.7-0.8]","[0.8-0.9]", "[0.9-1.0]"), length(concordance)/10))
    data_fused <- as.data.frame(cbind(chunk, as.vector(t(t(genotypes)/total)*100), as.vector(concordance)))
    colnames(data_fused) <- c("chunk", "genotypes", "concordance")
    data_fused_upper <- data.frame(chunk_upper, as.vector(t(t(genotypes_upper)/total*100)), as.vector(concordance_upper), MPP)
    colnames(data_fused_upper) <- c(colnames(data_fused), "MPP")
    
    p <- vector(mode="list", length=nrow(numbers))
    q <- vector(mode="list", length=nrow(numbers))
    r <- vector(mode="list", length=nrow(numbers))
   
   if (any(chunks_wo_snp)) {
     total <- total[-chunks_wo_snp]
   }
    
    pdf(paste(file, ".qc_perChunk.pdf", sep=""), onefile = TRUE, paper = "a4r", width = 12, height = 6)
    for (a in 1:nrow(numbers_red)) {
        data_chunk <- rbind(c(0.0,0.0,data_list_red[[a]][1,3], data_list_red[[a]][1,4]), data_list_red[[a]])
        colnames(data_chunk) <- c("Start", "Maximum posterior probability","Number of genotypes", "Concordance [%]")
       # p[[a]] <-  ggplot() +
        p <- ggplot() +
        geom_point(data=data_chunk, aes(x = data_chunk[,1]+0.05, y=data_chunk[,4], size = data_chunk[,3]/total[a]*100), alpha = 0.4) +
        geom_step(data=data_chunk, aes(x = data_chunk[,2] , y=data_chunk[,4]), direction="vh") +
        geom_hline(aes(yintercept=numbers_red[a,5]), color="red") +
        annotate("text", label="Overall concordance", x=0.0, y=102, color="red", hjust = 0) +
        annotate("text", label=paste("SNPs total:",numbers_red[a,1]), x=0.0, y=80, color="black", hjust = 0) +
        annotate("text", label=paste("SNPs only ref:", numbers_red[a,2]), x=0.0, y=75, color="black", hjust = 0) + 
        annotate("text", label=paste("SNPs only sample:", numbers_red[a,3]), x=0.0, y=70, color="black", hjust = 0) + 
        annotate("text", label=paste("SNPs both ref and sample:", numbers_red[a,4]), x=0.0, y=65, color="black", hjust = 0) + 
        scale_x_continuous(breaks=seq(0.0,1,0.1)) +
        scale_size_continuous(name="Percentage of genotypes")+
        xlab("Maximum posterior probability") +
        ylab("Concordance [%]") +
        labs(title = expression("Percentage of SNPs within bins of 0.1 maximal posterior probability and respective concordance")) +
        element_wo_bg()
        print(p)
    }
    dev.off()
    
 #   pdf(paste(file, ".qc_perChrPercent.pdf", sep=""), onefile = TRUE, paper = "a4r", width = 12, height = 6)
 #   q <-  ggplot() +
 #   geom_point(data=data_fused, aes(x = chunk, y = concordance, size = genotypes), alpha=0.4) +
 #   xlab(paste("Chromosome", i, "[chunks]", sep=" ")) +
 #   ylab("Concordance [%]") +
 #   scale_size_continuous(name="Percentage of SNPs") +
 #   labs(title = expression("Percentage of SNPs within bins of 0.1 maximal posterior probability (vertical orientation of dots)\nacross chromosomal chunks")) +
 #   element_wo_bg()
 #   print(q)
 #   dev.off()
    if (i < 24) {
        pdf(paste(file, ".qc_perChrMPP.pdf", sep=""), onefile = TRUE, paper = "a4r", width = 12, height = 6)
        r  <- ggplot() +
        geom_bar(data = data_fused_upper, aes(x = chunk*3 + chr_begin[i]/1000000, y = genotypes, fill = MPP), stat="identity", position="dodge") +
        geom_line(data=numbers, aes(x=seq(1,nrow(numbers),1)*3+chr_begin[i]/1000000, y=numbers[,5]),color="red" ) +
        geom_line( aes(x=c(centromere[i,3]/1000000, centromere[i,4]/1000000), y=c(-2,-2)),color="darkblue" ) +
        annotate("text", label=paste("Centromere"), x=(centromere[i,3]/1000000 + centromere[i,4]/1000000)/2, y=-4, color="darkblue") + 
        xlab(paste("Chromosome", chr, "[Mbp]", sep=" ")) +
        ylab("Number of SNPs [%]") +
        labs(title = expression("Percentage of SNPs with maximum posterior probability (MPP) greater than 70%")) +
        scale_fill_manual(values = c("tan2", "navyblue","springgreen4")) +
        annotate("text", label=paste("Cumulative concordance [%]"), x=chr_begin[i]/1000000, y=102, color="red", hjust = 0) + 
        annotate("text", label=paste("Original chunks without SNPs:\n", chunks_wo_snp_str ), x=nrow(numbers)*3+chr_begin[i]/1000000, y=106, color="black", hjust = 1) + 
        element_wo_bg()
        print(r)
        dev.off()
    }
    if (i >= 24) {
        pdf(paste(file, ".qc_perChrMPP.pdf", sep=""), onefile = TRUE, paper = "a4r", width = 12, height = 6)
        r  <- ggplot() +
        geom_bar(data = data_fused_upper, aes(x = chunk*(X_PAR[1]/1000000+1), y = genotypes, fill = MPP), stat="identity", position="dodge") +
        geom_line(aes(x=c(X_PAR[1]/1000000+0.8, X_PAR[1]/1000000+1.4), y=unlist(c(numbers[5], numbers[5]))),color="red" ) +
        xlab(paste("Chromosome", chr, "[Mbp]", sep=" ")) +
        ylab("Number of SNPs [%]") +
        labs(title = expression("Percentage of SNPs with maximum posterior probability (MPP) greater than 70%")) +
        scale_fill_manual(values = c("tan2", "navyblue","springgreen4")) +
        annotate("text", label=paste("Cumulative concordance [%]"), x=chr_begin[i]/1000000, y=102, color="red", hjust = 0) + 
        annotate("text", label=paste("Original chunks without SNPs:\n", chunks_wo_snp_str ), x=X_PAR[1]/1000000+1, y=106, color="black", hjust = 0) + 
        element_wo_bg() +
        theme(axis.text.x = element_blank())
        print(r)
        dev.off()
    }
}

system(paste("PDFconcat --verbose --output ", outputdir,"/imputeQC.perChrMPP.UK10K1000Genomes.pdf ", outputdir, "/*perChrMPP.pdf", sep=""))
system(paste("PDFconcat --verbose --output ", outputdir,"/imputeQC.perChunk.UK10K1000Genomes.pdf ", outputdir, "/*perChrMPP.pdf", sep=""))

# overview of SNPS after each filtering step
overview_file <- read.table(paste(inputdir, "/SNPsPerChr.txt", sep=""), sep="\t", header=TRUE,stringsAsFactors=FALSE)
pdf(paste(outputdir, "/SNPsperChr.pdf", sep=""), onefile = TRUE, paper = "a4r", width = 14, height = 6)
plot_object <- data.frame(chr = rep(seq(1, nrow(overview_file)-1,1), each=3), SNP = as.vector(rbind(overview_file$SNPs.after.Imputation[-26],overview_file$SNPs.after.Imputation.filtering.INFO...0.4[-26],overview_file$SNPs.after.HWE.and.MAF.filtering[-26])),filter = factor(rep(c("Imputed SNPs","Imputed SNPs after imputation filtering (INFO > 0.4)","Imputed SNPs after HWE and MAF filtering"),25), levels=c("Imputed SNPs","Imputed SNPs after imputation filtering (INFO > 0.4)","Imputed SNPs after HWE and MAF filtering"), ordered=TRUE))
p <- ggplot() +
geom_bar(data= plot_object, aes(x= chr, y= SNP, fill= filter), stat="identity", position="dodge") +
#geom_line(aes(x=seq(1, nrow(overview_file)-1,1),y=overview_file$SNPs.after.Imputation[-26])) +
#geom_line(aes(x=seq(1, nrow(overview_file)-1,1),y=overview_file$SNPs.after.Imputation.filtering.INFO...0.4[-26]),col="blue") +
#geom_line(aes(x=seq(1, nrow(overview_file)-1,1),y=overview_file$SNPs.after.HWE.and.MAF.filtering[-26]), col="green") +
#annotate("text", label="Imputed SNPs", x=10, y=3.5e6,  hjust = 0) +
#annotate("text", label="Imputed SNPs after imputation filtering (INFO > 0.4)", x=10, y=3.3e6, color="blue", hjust = 0) +
#annotate("text", label="Imputed SNPs after HWE and MAF filtering", x=10, y=3.1e6, color="green", hjust = 0) +
xlab("Chromosomes") +
ylab("Number of SNPs") +
labs(title ="Overview of SNP numbers after imputation, imputation filtering and\nMAF/HWE filtering") +
element_wo_bg() +
scale_x_continuous(limits=c(0,26),breaks=1:25, labels=paste("chr",overview_file$Chr[-26], sep="")) +
scale_fill_brewer()+
theme( axis.text.x  = element_text(angle=90, hjust=1))+ 
theme(legend.position=c(1,3000000))
print(p)
dev.off()

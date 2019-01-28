#############
# libraries #
#############

options(import.path="/homes/hannah/RPackagesDevel/plinkQC")
options(bitmapType = 'cairo', device = 'pdf')

optparse <- modules::import_package('optparse')
plinkqc <- modules::import_package('plinkQC')


################
### Analysis ###
################

## command line arguments ####
option_list <- list(
    optparse$make_option(c("-d", "--directory"), action="store",
               dest="qcdir",
               type="character", help="Path to directory with QC data/raw data
               results [default: %default].", default=NULL),
    optparse$make_option(c("-N", "--name"), action="store", dest="alg",
               type="character", help="Prefix of plink-files i.e.
               prefix.bim/prefix.bed/prefix.fam [default: %default].",
               default=NULL),
    optparse$make_option(c("--check_sex"), action="store_true",
               dest="do.check_sex", help="Identify individuals with
               discordant sex information (genetic sex/assigned sex)
               [default: %default].", default=FALSE),
    optparse$make_option(c("--maleTh"), action="store", dest="maleTh",
               type="double", help="Threshold for homozygosity
               rates expected on X-chromosome for males [default: %default].",
               default=0.8),
    optparse$make_option(c("--femaleTh"), action="store", dest="femaleTh",
               type="double", help="Threshold for homozygosity
               rates expected on X-chromosome for females [default: %default].",
               default=0.2),
    optparse$make_option(c("--fixMixup"), action="store_true",
               dest="fixMixup", help="Should PEDSEX of individuals with
               mismatch between PEDSEX and externalSexSex when
               externalSexSex==SNPSEX automatically corrected: this will directly
               change the alg.bim/.bed/.fam files! [default: %default].",
               default=FALSE),
    optparse$make_option(c("--externalSex"), action="store",
               dest="externalSexFile",
               type="character", help="/path/to/file/with/external/sex/id
               containing sample IDs [externalSexID] and sex [externalSexSex]
               to double check if external and PEDSEX data (often processed at
               different centers) match. Only required if --fixMixup.
               [default: %default].", default=NULL),
    optparse$make_option(c("--externalSexID"), action="store",
               dest="externalSexID",
               type="character", help="Column identifier for column containing
               ID information in externalSex [default: %default].",
               default="IID"),
    optparse$make_option(c("--externalSexSex"), action="store",
               dest="externalSexSex",
               type="character", help="Column identifier for column containing
               sex information in externalSex [default: %default].",
               default="Sex"),
    optparse$make_option(c("--externalMale"), action="store", dest="externalMale",
               type="character", help="Identifier for 'male' in externalSexSex,
               e.g. 1 or 'male' [default: %default].", default='M'),
    optparse$make_option(c("--externalFemale"), action="store",
               dest="externalFemale", type="character", help="Identifier for
               'female' in externalSexSex file, e.g. 2 or 'female'
               [default: %default].", default='F'),
    optparse$make_option(c("--check_het_and_miss"), action="store_true",
               dest="do.check_het_and_miss", help="Identify individuals
               with outlying missing genotype or heterozygosity rates
               [default: %default].", default=FALSE),
    optparse$make_option(c("--imissTh"), action="store", dest="imissTh",
               type="double", help="Threshold for acceptable missing genotype
               rate per individual [default: %default].", default=0.03),
    optparse$make_option(c("--hetTh"), action="store", dest="hetTh",
               type="double", help="Threshold for acceptable heterozyosity rates;
               any sample outside het*standard deviation(heterozygosity) will
               be excluded [default: %default].", default=3),
    optparse$make_option(c("--check_relatedness"), action="store_true",
               dest="do.check_relatedness", help="Identify related
               individuals [default: %default].", default=FALSE),
    optparse$make_option(c("--filterRelated"), action="store_true",
               dest="filterRelated", help="Set flag if  related
               individuals should be removed from final dataset
               [default: %default].", default=FALSE),
    optparse$make_option(c("--relatednessTh"), action="store", dest="highIBDTh",
               type="double", help="Identity be descent threshold, i.e.
               threshold for filtering of relatedness [default: %default].",
               default=0.1875),
    optparse$make_option(c("--check_ancestry"), action="store_true",
               dest="do.check_ancestry", help="Identify genetic
               ancestry of individuals [default: %default].",
               default=FALSE),
    optparse$make_option(c("--europeanTh"), action="store", dest="europeanTh",
               type="double", help=" Scaling factor of radius to be drawn around
               centre of reference European samples, with study samples inside
               this radius considered to be of European descent and samples
               outside this radius of non-European descent. The radius is
               computed as the maximum Euclidean distance of reference Europeans
               to the centre of reference European samples [default: %default].",
               default=1.5),
    optparse$make_option(c("--prefixMerged"), action="store", dest="prefixMerged",
               help="Prefix of merged dataset (study and reference samples) used
               in plink --pca, resulting in prefixMerged.eigenvec [default:
               %default].", default=NULL),
    optparse$make_option(c("--refSamples"), action="store",
              dest="refSamplesFile", type="character", help="/path/to/file with
              sample identifiers [refSamplesIID] corresponding to reference IIDs
              in prefixMerged.eigenvec and population identifier [refSamplesPop]
              corresponding to population IDs [refColorsPop] in refColors.
              [default: %default].", default=NULL),
    optparse$make_option(c("--refColors"), action="store",
              dest="refColorsFile", type="character", help="/path/to/file with
              population IDs in column [refColorsPop] and corresponding color-
              code for PCA plot in column [refColorsColor] [default: %default].",
              default=NULL),
    optparse$make_option(c("--refSamplesIID"), action="store",
              dest="refSamplesIID", type="character", help="Column name of
              reference sample IDs in --refSamples [default: %default].",
              default="IID"),
    optparse$make_option(c("--refSamplesPop"), action="store",
              dest="refSamplesPop", type="character", help="Column name of
              reference sample population IDs in --refSamples [default:
              %default].", default="Pop"),
    optparse$make_option(c("--refColorsColor"), action="store",
              dest="refColorsColor", type="character", help="Column name of
              population colors in --refColors [default: %default].",
              default="Color"),
    optparse$make_option(c("--refColorsPop"), action="store",
              dest="refColorsPop", type="character", help="Column name of
              reference sample population in --refColors [default:
              %default].", default="Pop"),
    optparse$make_option(c("--studyColor"), action="store",
              dest="studyColor", type="character", help="Colour to be used for
              study population if PCA plot [default: %default].",
              default="#2c7bb6"),
    optparse$make_option(c("--lmissTh"), action="store", dest="lmissTh",
               type="double", help="Threshold for acceptable missing genotype
               rate per marker [default: %default].", default=0.01),
    optparse$make_option(c("--hweTh"), action="store", dest="hweTh",
               type="double", help="Threshold for considering deviation from
               Hardy-Weinberg equilibrium significant [default: %default].",
               default=1e-5),
    optparse$make_option(c("--mafTh"), action="store", dest="mafTh",
               type="double", help="Minor allele frequency treshold [default:
               %default].", default=0.01),
    optparse$make_option(c("--macTh"), action="store", dest="macTh",
               type="integer", help="Minor allele count treshold. If both --maf
               and --mac are supplied, --maf will be considered [default:
               %default].", default=20),
    optparse$make_option(c("--showProgress"), action="store_true",
               dest="verbose",
               default=FALSE, type="logical", help="If set, progress messages
               about analyses are printed to standard out [default: %default]."),
    optparse$make_option(c("-p", "--plot"), action="store_true",
               dest="plot", help="Plot QC results [default: %default].",
               default=FALSE),
    optparse$make_option(c("--plink"), action="store", dest="path2plink",
               help="Absolute path to where external plink software
               (https://www.cog-genomics.org/plink/1.9/) can be found. If not
               provided, assumed that PATH set-up works and plink will be found
               by system('plink'). [default: %default].", default=NULL),
    optparse$make_option(c("--debug"), action="store_true",
                        dest="debug", default=FALSE, type="logical",
                        help="If set, predefined arguments are used to test the
                        script [default: %default].")
)

args <- optparse$parse_args(optparse$OptionParser(option_list=option_list))

if (args$debug) {
    # testing parser
    args_vec <-
        c("--name=gencall.combined",
          "--directory=/homes/hannah/data/digital-heart/genotype/QC/combined",
          "--check_sex", "--fixMixup",
          "--maleTh=0.8","--femaleTh=0.2",
          "--externalSex=/homes/hannah/data/digital-heart/phenotype/2Dphenotype/20160209_All_BRU_family_format.txt",
          "--externalSexSex=Sex",
          "--externalSexID=Bru.Number",
          "--check_het_and_miss", "--imissTh=0.03", "--hetTh=3",
          "--check_relatednes", "--relatednessTh=0.1875",
          "--prefixMerge=gencall.combined.HapMapIII",
          "--refSamples=/homes/hannah/data/hmeyer/HapMap/HapMap_ID2Pop.txt",
          "--refColors=/homes/hannah/data/hmeyer/HapMap/HapMap_PopColors.txt",
          "--lmissTh=0.01", "--hweTh=1e-5", "--macTh=20",
          "--plink=/homes/hannah/bin/plink", "--plot", "--showProgress", "--check_ancestry")
    args <- optparse$parse_args(optparse$OptionParser(option_list=option_list),
                                arg=args_vec)
}


if (args$plot) pdf(paste(args$qcdir,"/", args$alg,".pdf",sep=""), width=8.3,
                   height=11.7, onefile=TRUE)
## Per-individual QC ####
if (args$verbose) message("per-individual QC")
SampleID <- data.table::fread(file=args$externalSexFile, header=TRUE,
                       quote="", stringsAsFactors=FALSE, na.strings=c("NA",""),
                       data.table=FALSE, fill=TRUE)
fail_individuals <-
    plinkqc$perIndividualQC(indir=args$qcdir, name=args$alg,
                        maleTh=args$maleTh,
                        femaleTh=args$femaleTh,
                        externalSex=SampleID,
                        fixMixup=args$fixMixup,
                        externalMale=args$externalMale,
                        externalFemale=args$externalFemale,
                        externalSexSex=args$externalSexSex,
                        externalSexID=args$externalSexID,
                        imissTh=args$imissTh,
                        hetTh=args$hetTh,
                        highIBDTh=args$highIBDTh,
                        prefixMergedDataset=args$prefixMerged,
                        do.run_check_ancestry=FALSE,
                        refSamplesFile=args$refSamplesFile,
                        refColorsFile=args$refColorsFile,
                        refSamplesIID=args$refSamplesIID,
                        refSamplesPop=args$refSamplesPop,
                        refColorsColor=args$refColorsColor,
                        refColorsPop=args$refColorsPop,
                        studyColor=args$studyColor,
                        path2plink=args$path2plink,
                        interactive=args$plot, verbose=args$verbose)

overview <- plinkqc$overviewPerIndividualQC(fail_individuals, interactive=args$plot)

## Per-marker QC ####
if (args$verbose) message("per-marker QC")
fail_markers <- plinkqc$perMarkerQC(indir=args$qcdir, name=args$alg,
                                    mafTh=args$mafTh, macTh=args$macTh,
                                    hweTh=args$hweTh, lmissTh=0.01,
                                    path2plink=args$path2plink,
                                    interactive=args$plot, verbose=args$verbose)
if (args$plot) dev.off()

## Create QCed dataset ####
if (args$verbose) message("Clean data")
plinkqc$clean_data(qcdir=args$qcdir, alg=args$alg, mafTh=args$mafTh,
                          macTh=args$macTh, hweTh=args$hweTh, lmissTh=0.01,
                          path2plink=args$path2plink,
                          filterRelated=args$filterRelated, verbose=args$verbose)

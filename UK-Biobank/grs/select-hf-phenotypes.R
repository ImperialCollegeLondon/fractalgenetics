#################
## libraries ####
#################
options(bitmapType = 'cairo', device = 'pdf')

modules::import_package('dplyr', attach=TRUE)
modules::import_package('tidyverse', attach=TRUE)
modules::import_package('bit64', attach=TRUE)
optparse <- modules::import_package('optparse')
data.table <- modules::import_package('data.table')

#################################
## parameters and input data ####
#################################
option_list <- list(
    make_option(c("-o", "--outdir"), action="store", dest="outdir",
               type="character", help="Path to output directory
               [default: %default].", default=NULL),
    make_option(c("-d", "--data"), action="store", dest="ukbdata",
               type="character", help="Path to ukb bulk tab-separated data file
               [default: %default].", default=NULL),
    make_option(c("--debug"), action="store_true",
               dest="debug", default=FALSE, type="logical",
               help="If set, predefined arguments are used to test the script
               [default: %default].")
)

args <- optparse$parse_args(OptionParser(option_list=option_list))

if (args$debug) {
    args <- list()
    args$outdir <- "~/data/ukbb/ukb-hrt"
    args$ukbdata <- "~/data/ukbb/ukb-hrt/rawdata/ukb40616.txt"
}

############
## data ####
############

## ukbb bulk data ####
bb <- data.table$fread(args$ukbdata, header=TRUE, sep="\t", fill=TRUE)

################
## analysis ####
################

# the important columns for phenotyping are:
#20002 - patient reported non-cancer diagnoses
#20004 - patient reported operations (CABG, triple bypass, angioplasty +/- stent)
#41270 - diagnoses codes recorded during hospital episodes (ICD10) - MISSING FIELD
#41271 - diagnoses codes recorded during hospital episodes (ICD9) - MISSING FIELD
#40001 - primary cause of death
#40002 - contributory/secondary cause of death
#41202 - main diagnoses codes recorded during all hospital episodes (ICD10)
#41203 - main diagnoses codes recorded during all hospital episodes (ICD9)
#41204 - secondary diagnosis codes recorded during all hospital episodes (ICD10)
#41205 - secondary diagnosis codes recorded during all hospital episodes (ICD9)
#41272 - operative procedures (OPCS4) - MISSING FIELD
#42000 - date of MI (date vs. NA)
#6150 - vascular problems diagnosed by doctor

#code to select relevant columns and subset data frame
select_eid <- bb %>%
  select(starts_with("f.eid"))
select_20002 <- bb%>%
  select(starts_with("f.20002."))
select_20004 <- bb%>%
  select(contains("f.20004."))
select_41202 <- bb %>%
  select(starts_with("f.41202."))
select_41203 <- bb %>%
  select(starts_with("f.41203."))
select_41204 <- bb %>%
  select(starts_with("f.41204."))
select_41205 <- bb %>%
  select(starts_with("f.41205."))
select_40001 <- bb %>%
  select(starts_with("f.40001."))
select_40002 <- bb %>%
  select(starts_with("f.40002."))
select_42000 <- bb %>%
  select(starts_with("f.42000."))
select_6150 <- bb %>%
  select(starts_with("f.6150."))
select_2443 <- bb %>%
  select(starts_with("f.2443."))

# bind columns to form focused data frames
df <- bind_cols(select_eid, select_40001, select_40002, select_41202,
                select_41204,select_41203, select_41205,select_20002,
                select_20004)
df_dr <- bind_cols(select_eid, select_6150)
df_diabetes <- bind_cols(select_eid, select_2443)


# column names as character string for search/filtering
#col_names <- colnames(select_df)

#[,"f.eid"] - use to identify which eid return from filtered search

#Can use uniqueN(x) to identify number of unique rows

## HEART FAILURE ####

#Heart failure ICD_10
I11.0 <- which(df == "I110", arr.ind=T)
I11.0_row <- I11.0[,1]
I13.0 <- which(df == "I130", arr.ind=T)
I13.0_row <- I13.0[,1]
I13.2 <- which(df == "I132", arr.ind=T)
I13.2_row <- I13.2[,1]
I25.5 <- which(df == "I255", arr.ind=T)
I25.5_row <- I25.5[,1]
I42.0 <- which(df == "I420", arr.ind=T)
I42.0_row <- I42.0[,1]
I42.5 <- which(df == "I425", arr.ind=T)
I42.5_row <- I42.5[,1]
I42.8 <- which(df == "I428", arr.ind=T)
I42.8_row <- I42.8[,1]
I42.9 <- which(df== "I429", arr.ind=T)
I42.9_row <- I42.9[,1]
I50.0 <- which(df == "I500", arr.ind=T)
I50.0_row <- I50.0[,1]
I50.1 <- which(df== "I501", arr.ind=T)
I50.1_row <- I50.1[,1]
I50.9 <- which(df == "I509", arr.ind=T)
I50.9_row <- I50.9[,1]
hf_icd10_total <- c(I11.0_row, I13.0_row, I13.2_row, I25.5_row, I42.0_row,
                    I42.5_row, I42.8_row, I42.9_row, I50.0_row, I50.1_row,
                    I50.9_row)
hf_icd10_total_unique <- unique(hf_icd10_total)

# Heart failure ICD_9
I425 <- which(df == "425", arr.ind=T)
I425_row <- I425[,1]
I4254 <- which(df == "4254", arr.ind=T)
I4254_row <- I4254[,1]
I4280 <- which(df == "4280", arr.ind=T)
I4280_row <- I4280[,1]
I428 <- which(df == "428", arr.ind=T)
I428_row <- I428[,1]
I4281 <- which(df == "4281", arr.ind=T)
I4281_row <- I4281[,1]
I4289 <- which(df == "4289", arr.ind=T)
I4289_row <- I4289[,1]
hf_icd9_total <- c(I4254_row, I4280_row, I4281_row, I4289_row, I425_row,
                   I428_row)
hf_icd9_total_unique <- unique(hf_icd9_total)

# Heart failure from interview
hf_interview1 <- which(df=="1076", arr.ind=T)
hf_interview1_row <- hf_interview1[,1]
hf_interview2 <- which(df=="1079", arr.ind=T)
hf_interview2_row <- hf_interview2[,1]
hf_interview_total <- c(hf_interview1_row, hf_interview2_row)
hf_interview_total_unique <- unique(hf_interview_total)

# HCM
I42.1 <- which(df == "I421", arr.ind=T)
I42.1_row <- I42.1[,1]
I42.2 <- which(df == "I422", arr.ind=T)
I42.2_row <- I42.2[,1]
hcm_interview <- which(df=="1588", arr.ind=T)
hcm_interview_row <- hcm_interview[,1]
hcm_total <- c(I42.1_row, I42.2_row, hcm_interview_row)
hcm_total_unique <- unique(hcm_total)

# Pooled heart failure IDs
hf_total <- c(hf_icd10_total_unique, hf_icd9_total_unique,
              hf_interview_total_unique)
hf_total_unique_with_hcm <- unique(hf_total)
hf_total_unique <- setdiff(hf_total_unique_with_hcm, hcm_total_unique)

# NEED TO CHANGE from I21, I22, I23 - READ THE NATURE GENETICS PAPER 
# https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19&nl=1

## CORONARY ARTERY DISEASE ####
# CAD ICD_10
I21.0 <- which(df == "I210", arr.ind=T)
I21.0_row <- I21.0[,1]
I21.1 <- which(df == "I211", arr.ind=T)
I21.1_row <- I21.1[,1]
I21.2 <- which(df == "I212", arr.ind=T)
I21.2_row <- I21.2[,1]
I21.3 <- which(df == "I213", arr.ind=T)
I21.3_row <- I21.3[,1]
I21.4 <- which(df == "I214", arr.ind=T)
I21.4_row <- I21.4[,1]
I21.9 <- which(df == "I219", arr.ind=T)
I21.9_row <- I21.9[,1]
I22.0 <- which(df == "I220", arr.ind=T)
I22.0_row <- I22.0[,1]
I22.1 <- which(df == "I221", arr.ind=T)
I22.1_row <- I22.1[,1]
I22.8 <- which(df == "I228", arr.ind=T)
I22.8_row <- I22.8[,1]
I22.9 <- which(df == "I229", arr.ind=T)
I22.9_row <- I22.9[,1]
# I23.1 <- which(df_icd10 == "I231", arr.ind=T)
# I23.1_row <- I23.1[,1]
# I23.2 <- which(df_icd10 == "I232", arr.ind=T)
# I23.2_row <- I23.2[,1]
# I23.3 <- which(df_icd10 == "I233", arr.ind=T)
# I23.3_row <- I23.3[,1]
# I23.6 <- which(df_icd10 == "I236", arr.ind=T)
# I23.6_row <- I23.6[,1]
# I23.8 <- which(df_icd10 == "I238", arr.ind=T)
# I23.8_row <- I23.8[,1]
# I24.1 <- which(df_icd10 == "I241", arr.ind=T)
# I24.1_row <- I24.1[,1]
I25.2 <- which(df == "I252", arr.ind=T)
I25.2_row <- I25.2[,1]
cad_icd10_total <- c(I21.0_row, I21.1_row, I21.2_row, I21.3_row, I21.4_row,
                     I21.9_row, I22.0_row, I22.1_row, I22.8_row, I22.9_row,
                     I25.2_row)
cad_icd10_total_unique <- unique(cad_icd10_total)

# CAD ICD_9
I410 <- which(df == "410", arr.ind=T)
I410_row <- I410[,1]
I411 <- which(df == "411", arr.ind=T)
I411_row <- I411[,1]
I412 <- which(df == "412", arr.ind=T)
I412_row <- I412[,1]
cad_icd9_total <- c(I410_row, I411_row, I412_row)
cad_icd9_total_unique <- unique(cad_icd9_total)

# CAD from nurse interview
#cad_interview1 <- which(df=="1074", arr.ind=T) - ANGINA ONLY
#cad_interview1_row <- cad_interview1[,1]
cad_interview2 <- which(df=="1075", arr.ind=T)
cad_interview2_row <- cad_interview2[,1]
cad_interview3 <- which(df=="1070", arr.ind=T)
cad_interview3_row <- cad_interview3[,1]
cad_interview4 <- which(df=="1095", arr.ind=T)
cad_interview4_row <- cad_interview4[,1]
cad_interview5 <- which(df=="1523", arr.ind=T)
cad_interview5_row <- cad_interview5[,1]
cad_interview_total <- c(cad_interview2_row, cad_interview3_row,
                         cad_interview4_row, cad_interview5_row)
cad_interview_total_unique <- unique(cad_interview_total)

# CAD from doctor interview
cad_dr <- which(df_dr == "1", arr.ind=T)
cad_dr_row <- cad_dr[,1]

# Pooled CAD
cad_total <- c(cad_icd10_total_unique, cad_icd9_total_unique,
               cad_interview_total_unique, cad_dr_row)
cad_total_unique <- unique(cad_total)

# HTN
# I10 <- which(df == "I10", arr.ind=T)
# I10_row <- I10[,]
# I11 <- which(df == "I11", arr.ind=T)
# I11_row <- I11[,]
# I12 <- which(df == "I12", arr.ind=T)
# I12_row <- I12[,]
# I13 <- which(df == "I13", arr.ind=T)
# I13_row <- I13[,]
# I15 <- which(df == "I15", arr.ind=T)
# I15_row <- I15[,]
# htn_icd10_total <- c(I10_row, I11_row, I12_row, I13_row, I15_row)
# htn_icd10_total_unique <- unique(htn_icd10_total)
#
# I401 <- which(df == "401", arr.ind=T)
# I401_row <- I401[,]
# I402 <- which(df == "402", arr.ind=T)
# I402_row <- I402[,]
# I403 <- which(df == "403", arr.ind=T)
# I403_row <- I403[,]
# I404 <- which(df == "404", arr.ind=T)
# I404_row <- I404[,]
# I405 <- which(df == "405", arr.ind=T)
# I405_row <- I405[,]
# htn_icd9_total <- c(I401_row, I402_row, I403_row, I404_row, I405_row)
# htn_icd9_total_unique <- unique(htn_icd9_total)
# 
# htn_interview1 <- which(df_interview=="1072", arr.ind=T)
# htn_interview1_row <- htn_interview1[,1]
# htn_interview2 <- which(df_interview=="1065", arr.ind=T)
# htn_interview2_row <- htn_interview2[,1]
# htn_interview_total <- c(htn_interview1_row, htn_interview2_row)
# htn_interview_total_unique <- unique(htn_interview_total)
# 
# htn_total <- c(htn_icd10_total_unique, htn_icd9_total_unique,
# htn_interview_total_unique)
# htn_total_unique <- unique(htn_total)
#
# T2 Diabetes
E11.0 <- which(df == "E110", arr.ind=T)
E11.0_row <- E11.0[,1]
E11.1 <- which(df == "E111", arr.ind=T)
E11.1_row <- E11.1[,1]
E11.2 <- which(df == "E112", arr.ind=T)
E11.2_row <- E11.2[,1]
E11.3 <- which(df == "E113", arr.ind=T)
E11.3_row <- E11.3[,1]
E11.4 <- which(df == "E114", arr.ind=T)
E11.4_row <- E11.4[,1]
E11.5 <- which(df == "E115", arr.ind=T)
E11.5_row <- E11.5[,1]
E11.6 <- which(df == "E116", arr.ind=T)
E11.6_row <- E11.6[,1]
E11.7 <- which(df == "E1107", arr.ind=T)
E11.7_row <- E11.7[,1]
E11.8 <- which(df == "E118", arr.ind=T)
E11.8_row <- E11.8[,1]
E11.9 <- which(df == "E119", arr.ind=T)
E11.9_row <- E11.9[,1]

t2dm_icd10_total <- c(E11.0_row, E11.1_row, E11.2_row, E11.3_row, E11.4_row,
                      E11.5_row, E11.6_row, E11.7_row, E11.8_row, E11.9_row)
t2dm_icd10_total_unique <- unique(t2dm_icd10_total)
t2dm_interview_total <- which(df=="1223", arr.ind=T)

t2dm_total <- c(t2dm_icd10_total_unique, t2dm_interview_total)
t2dm_total_unique <- unique(t2dm_total)

## All diabetes
E10.0 <- which(df== "E100", arr.ind=T)
E10.0_row <- E10.0[,1]
E10.1 <- which(df == "E101", arr.ind=T)
E10.1_row <- E10.1[,1]
E10.2 <- which(df== "E102", arr.ind=T)
E10.2_row <- E10.2[,1]
E10.3 <- which(df == "E103", arr.ind=T)
E10.3_row <- E10.3[,1]
E10.4 <- which(df == "E104", arr.ind=T)
E10.4_row <- E10.4[,1]
E10.5 <- which(df == "E105", arr.ind=T)
E10.5_row <- E10.5[,1]
E10.6 <- which(df == "E106", arr.ind=T)
E10.6_row <- E10.6[,1]
E10.7 <- which(df == "E107", arr.ind=T)
E10.7_row <- E10.7[,1]
E10.8 <- which(df == "E108", arr.ind=T)
E10.8_row <- E10.8[,1]
E10.9 <- which(df == "E109", arr.ind=T)
E10.9_row <- E10.9[,1]
E14.0 <- which(df == "E140", arr.ind=T)
E14.0_row <- E14.0[,1]
E14.1 <- which(df == "E141", arr.ind=T)
E14.1_row <- E14.1[,1]
E14.2 <- which(df == "E142", arr.ind=T)
E14.2_row <- E14.2[,1]
E14.3 <- which(df == "E143", arr.ind=T)
E14.3_row <- E14.3[,1]
E14.4 <- which(df == "E144", arr.ind=T)
E14.4_row <- E14.4[,1]
E14.5 <- which(df == "E145", arr.ind=T)
E14.5_row <- E14.5[,1]
E14.6 <- which(df == "E146", arr.ind=T)
E14.6_row <- E14.6[,1]
E14.7 <- which(df == "E147", arr.ind=T)
E14.7_row <- E14.7[,1]
E14.8 <- which(df == "E148", arr.ind=T)
E14.8_row <- E14.8[,1]
E14.9 <- which(df == "E149", arr.ind=T)
E14.9_row <- E14.9[,1]

diabetes_other <- c(E10.0_row, E10.1_row, E10.2_row, E10.3_row, E10.4_row,
                    E10.5_row, E10.6_row, E10.7_row, E10.8_row, E10.9_row,
                    E14.0_row, E14.1_row, E14.2_row, E14.3_row, E14.4_row,
                    E14.5_row, E14.6_row, E14.7_row, E14.8_row, E14.9_row)

diabetes_dr <- which(df_diabetes == "1", arr.ind=T)
diabetes_dr_row <- diabetes_dr[,1]

dm_total <- c(diabetes_dr_row, diabetes_other, t2dm_total_unique)
dm_total_unique <-unique(dm_total)

## ATRIAL FIBRILLATION ####
# I48 <- which(df_icd10 == "I48", arr.ind=T)
# I48_row <- I48[,1]
# I48.0 <- which(df_icd10 == "I480", arr.ind=T)
# I48.0_row <- I48.0[,1]
# I48.1 <- which(df_icd10 == "I481", arr.ind=T)
# I48.1_row <- I48.1[,1]
# I48.2 <- which(df_icd10 == "I482", arr.ind=T)
# I48.2_row <- I48.2[,1]
# I48.3 <- which(df_icd10 == "I483", arr.ind=T)
# I48.3_row <- I48.3[,1]
# I48.4 <- which(df_icd10 == "I484", arr.ind=T)
# I48.4_row <- I48.4[,1]
# I48.9 <- which(df_icd10 == "I489", arr.ind=T)
# I48.9_row <- I48.9[,1]
# I4273 <- which(df_icd9 == "4273", arr.ind=T)
# I4273_row <- I4273[,1]
# af_interview1 <- which(df_interview == "1471", arr.ind=T)
# af_interview1_row <- af_interview1[,1]
# af_interview2 <- which(df_interview == "1483", arr.ind=T)
# af_interview2_row <- af_interview2[,1]
# af_interview3 <- which(df_interview== "1524", arr.ind=T)
# af_interview3_row <- af_interview3[,1]
#
# af_total <- c(I48_row, I48.0_row, I48.1_row, I48.2_row, I48.3_row, I48.4_row,
# I48.9_row, I4273_row,af_interview1_row,af_interview2_row,af_interview3_row )
# af_total_unique <- unique(af_total)

# patients with heart failure and no recorded cad diagnosis (i.e. assumed NICM)
nicm_pop <- setdiff(hf_total_unique, cad_total_unique)

# NICM as per Aragam paper definition
aragam_nicm_total <- c(I42.0_row, I50.1_row, I4281_row)
aragam_nicm_total_unique_with_hcm <- unique(aragam_nicm_total)
aragam_nicm_total_unique <- setdiff(aragam_nicm_total_unique_with_hcm,
                                    hcm_total_unique)
aragam_nicm_pop <- setdiff(aragam_nicm_total_unique, cad_total_unique)

# NICM as per SLZ definition
sz_nicm_total <- c(I42.0_row, I50.1_row, I50.0_row, I50.9_row, hf_icd9_total,
                   hf_interview_total)
sz_nicm_total_unique_with_hcm <- unique(sz_nicm_total)
sz_nicm_total_unique <- setdiff(sz_nicm_total_unique_with_hcm, hcm_total_unique)
sz_nicm_pop <- setdiff(sz_nicm_total_unique, cad_total_unique)

# ICM
icm_pop <- hf_total_unique[hf_total_unique %in% cad_total_unique]

# list of patient IDs
hf_eid <- df[hf_total_unique,1]
nicm_eid <- df[nicm_pop,1]
icm_eid <- df[icm_pop,1]
cad_eid <- df[cad_total_unique,1]
aragam_nicm_eid <- df[aragam_nicm_pop,1]
sz_nicm_eid <- df[sz_nicm_pop,1]
dm_eid <- df[dm_total_unique,1]
eid <- df[,1]

write.table(hf_eid, file=file.path(args$outdir, "hf_eid.csv"), sep=",",
            quote=FALSE, col.names="IID", row.names=FALSE)
write.table(nicm_eid, file=file.path(args$outdir, "nicm_eid.csv"), sep=",",
            quote=FALSE, col.names="IID", row.names=FALSE)
write.table(icm_eid, file=file.path(args$outdir, "icm_eid.csv"), sep=",",
            quote=FALSE, col.names="IID", row.names=FALSE)
write.table(cad_eid, file=file.path(args$outdir, "cad_eid.csv"), sep=",",
            quote=FALSE, col.names="IID", row.names=FALSE)
write.table(aragam_nicm_eid, file=file.path(args$outdir, "aragam_nicm_eid.csv"),
            sep=",",
            quote=FALSE, col.names="IID", row.names=FALSE)
write.table(sz_nicm_eid, file=file.path(args$outdir, "sz_nicm_eid.csv"), sep=",",
            quote=FALSE, col.names="IID", row.names=FALSE)
write.table(dm_eid, file=file.path(args$outdir, "dm_eid.csv"), sep=",",
            quote=FALSE, col.names="IID", row.names=FALSE)
write.table(eid, file=file.path(args$outdir, "eid.csv"), sep=",",
            quote=FALSE, col.names="IID", row.names=FALSE)



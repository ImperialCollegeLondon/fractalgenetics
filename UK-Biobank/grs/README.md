# UKBB filter

Code for filtering UK biobank data for heart failure phenotypes.

```
bb <- fread("filename.tab", header=TRUE, sep="\t", fill=TRUE)
```
Change \<filename\> to the `.tab` UKBB data file.

## UK Biobank data fields searched
  Patient self- reported non-cancer diagnosis (datafield 20002) and operations (datafield 2004) at initial nurse-led interview.
  
  Primary (datafield 40001) and contributory/secondary (datafield 40002) causes of death.
  
  Primary (ICD-10 - datafield 41202; ICD-9 - datafield 41203) and secondary (ICD-10 41204; ICD-9 datafield 41205) diagnosis codes recorded during all hospital admission episodes. 
  
  Patient self-reported previous vascular problems diagnosed by a doctor (datafield 6150). 
  
## UK Biobank disease phenotype definitions
**All cause heart failure**  
  Any heart failure (as defined below) excluding individuals with HCM (ICD10 codes I42.1 and I42.2, and self reported HCM at nurse interview).
  
**Non-ischaemic cardiomyopathy (or heart failure without known coronary artery disease)**  
  Any heart failure (as defined below) excluding individuals with HCM (ICD10 codes I42.1 and I42.2, and self reported HCM at nurse interview) and coronary artery disease (defined below).
  
**Non-ischaemic cardiomyopathy (SZ definition)**  
  Non-ischaemic cardiomyopathy (as defined below) excluding individuals with HCM (ICD10 codes I42.1 and I42.2, and self reported HCM at nurse interview) and coronary artery disease (defined below). 

**Non-ischaemic cardiomyopathy (Aragam definition)**  
 Hospitalisation or death due to the following codes: ICD-10: I42.0 (dilated cardiomyopathy) and I50.1 (left ventricular failure), ICD-9: left heart failure (4281); excluding individuals with HCM (ICD10 codes I42.1 and I42.2, and self reported HCM at nurse interview) and coronary artery disease (defined below). 

**Coronary artery disease**  
  Coronary artery disease (as defined below). 
  
## Any heart failure
**ICD-10 (40001, 400002, 41202, 41204)**  
  I11.0 (I110) – Hypertensive heart disease with HF  
  I13.0 (I130)  – Hypertensive heart and renal disease with HF  
  I13.2 (I132) - Hypertensive heart and renal disease with HF and renal failure  
  I25.5 (I255) – Ischaemic CM  
  I42.0 (I420)  – dilated CM  
  I42.5 (I425) – restrictive CM  
  I42.8 (I428) – other CM  
  I42.9 (I429) – CM, unspecified  
  I50.0 (I500) – congestive HF  
  I50.1 (I501) – left ventricular failure  
  I50.9 (I509) – HF, unspecified 
  
**ICD-9 (41203, 41205)**  
  428 – HF  
  425 – CM  
  4254 – other primary CM  
  4280 – congestive HF  
  4281 – left HF  
  4289 – HF unspecified  
  
**Self-reported (20002)**  
  1076 – HF/pulmonary oedema  
  1079 - CM  

## non-ischaemic cardiomyopathy
**ICD-10 (40001, 400002, 41202, 41204)**  
  I42.0 (I420)  – dilated CM  
  I50.0 (I500) – congestive HF  
  I50.1 (I501) – left ventricular failure  
  I50.9 (I509) – HF, unspecified 
  
**ICD-9 (41203, 41205)**  
  428 – HF  
  425 – CM  
  4254 – other primary CM  
  4280 – congestive HF  
  4281 – left HF  
  4289 – HF unspecified  
  
**Self-reported (20002)**  
  1076 – HF/pulmonary oedema  
  1079 - CM  

## coronary artery disease
**ICD-10 (40001, 400002, 41202, 41204)**  
  I21.0 (I210) – acute anterior MI  
  I21.1 (I211) – acute inferior MI  
  I21.2 (I212) – acute other site MI  
  I21.3 (I213) – acute unspecified site MI  
  I21.4 (I214) – acute subendocardial MI  
  I21.9 (I219) – acute MI unspecified  
  I22 (I22) – subsequent MI  
  I22.0 (I220) – subsequent anterior mI  
  I22.1 (I221) – subsequent inferior MI  
  I22.8 (I228) – subsequent other site MI  
  I22.9 (I229) – subsequent unspecified site MI  
  I25.2 – old MI  
  
**ICD-9 (41203, 41205)**  
  410  
  411  
  412  
  
**Self-reported (20002, 20004)**  
  1075 - MI  
  1095 – CABG  
  1523 – Triple heart bypass  
  1070 – Coronary angioplasty +/- stenting  
  
**Vascular problems diagnosed by dr (6150)**  
  1 – heart attack  

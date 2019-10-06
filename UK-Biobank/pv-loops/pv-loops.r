# Pressure-volume (PV) loops for UK Biobank data (UKBB) and finite element model data (FEM)

# Define a function to rearrange the PV loop data

reformat.for.PV.loop<- function(new.y, LVEDP1, LVEDP2)
{
  if (is.null(nrow(new.y))==FALSE) {
    new.coords.display<- t(new.y[,c(1,1,2,2,2,4,3,5,1,6)])
    new.coords.display[2,]<- LVEDP1
    new.coords.display[4,]<- LVEDP2
    rownames(new.coords.display)<- c("LVESV","LVEDP1","LVEDV","LVEDP2","LVEDV","DBP","LVMSV","SBP","LVESV","ESP")
    return (new.coords.display)}
  if (is.null(nrow(new.y))==TRUE) {
    new.coords.display<- new.y[c(1,1,2,2,2,4,3,5,1,6)]
    new.coords.display[2]<- LVEDP1
    new.coords.display[4]<- LVEDP2
    names(new.coords.display)<- c("LVESV","LVEDP1","LVEDV","LVEDP2","LVEDV","DBP","LVMSV","SBP","LVESV","ESP")
    return (new.coords.display)}
}


# Read in dummy data
    dataI<- read.table("dummydata.txt")

# Create a matrix with the pressure and volume landmarks included for a PV loop
    LVEDP1<- 4 
    LVEDP2<- 8 
    terms<- c("LVESV", "LVEDV","LVMSV", "PWA_DIASTOLICBP", "PWA_CENTRALSYSTOLICBP", "PWA_ENDSYSTOLICPRESSURE", "TRABMASS")
    landmark.mat<- dataI[,terms]


# Prepare data for Biobank analysis
    terms.reg<- c("MEANGLOBALFD", "CONTRACTILITY", "SVR", "TRABMASS", "HRfromCMR")
    df.sc<- data.matrix(dataI[, c(terms, terms.reg)])
    colnames(df.sc)<- c(terms, terms.reg)

# Fit linear model
    fit<- list()
    for (i in 1:length(terms)) {eqn<- as.formula(paste(paste(terms[i], collapse = "+"), "~", paste(terms.reg, collapse = "+")));
    fit[[i]]<-lm(eqn, data=data.frame(df.sc))}


# Now predict the fitted data for the 'smooth' and 'rough' cases

      # Independent variable for the two cases
          new.df<- data.frame(MEANGLOBALFD = c(1.000, 1.16934), CONTRACTILITY=c(1.5,1.61), SVR=c(1440,1440), TRABMASS=c(0,25), HRfromCMR=c(60,60))
      # Predict the new values for the PV loop
          new.y<- matrix(0, ncol=length(terms), nrow=2, dimnames=list(c("Smooth","Rough"), terms))
          for (i in 1:length(terms)) {new.y[,i]<- predict(fit[[i]], newdata=new.df)}
      # Convert to display format for PV loop  
          new.coords.display<- reformat.for.PV.loop(new.y, LVEDP1, LVEDP2)

          
# Now provide the FEM data for comparison
          FEM.coords.display<- matrix(c(54.8, 4,   114.8, 5,   114.8, 50,   61.9, 90,   54.8, 88,
                                        79.6, 2.1, 159.8, 4.7, 159.8, 65.6, 87.4, 122.6, 79.6, 116.7), ncol=2, byrow=F)
          
          
# Plot the PV loops for UKBB and FEM
    par(mfrow=c(1,2))
    plot.PV.loop(new.coords.arr, 50, 170, 0, 120, new.coords.display, 4, "blue")
    plot.PV.loop(new.coords.arr, 50, 170, 0, 120, FEM.coords.display, 4, "red")



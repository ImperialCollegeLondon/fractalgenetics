# Pressure-volume (PV) loops for UK Biobank data (UKBB) and finite element model data (FEM)

# Install and load packages for session
    packages<-  c("pracma","sf","smoothr")
    inst <- packages %in% installed.packages()
    if(length(packages[!inst]) > 0) install.packages(packages[!inst])
    lapply(packages, require, character.only=TRUE)
    
# Define functions

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
    
    plot.PV.loop<- function(a, b, c, d, new.coords.display, line.thickness, col)
    {
      plot(0, type='n', xlim=c(a,b), ylim=c(c,d), xaxt='n', yaxt='n', bty='n', xlab="", ylab="")
      axis(side=1, at=seq(a,b,20), labels=seq(a,b,20), lwd=2)
      axis(side=2, at=seq(c,d,20), labels=seq(c,d,20), lwd=2)
      mtext("Volume (ml)", side=1, line=2.5)
      mtext("Pressure (mmHg)", side=2, line=2.5)
     
      for (i in 1:2) {coords<- matrix(new.coords.display[,i],ncol=2, byrow=T)[c(1:5,1),]
                      lines(pvloop(coords), type='l', col=col, lty=i, lwd=line.thickness)}
    }
    
    pvloop<- function(coords)
    {
      shape<- matrix(unlist(smooth(st_sfc(st_polygon(list(coords))), method='ksmooth', max_distance=2, bandwidth=30)),ncol=2)
      f1<- findpeaks(-shape[,1], zero="+")
      f2<- findpeaks(shape[,1], zero="+")
      shape1<- rbind(shape[1:f2[1,2],],
                     shape[f2[1,2]:f1[1,2],],
                     shape[f1[1,2]:nrow(shape),])
      return(shape1)
    }

# Read in dummy data
    dataI<- read.table("dummydata.txt")

# Create a matrix with the pressure and volume landmarks included for a PV loop
    LVEDP1<- 4 
    LVEDP2<- 8 
    terms<- c("LVESV", "LVEDV","LVMSV", "PWA_DIASTOLICBP", "PWA_CENTRALSYSTOLICBP", "PWA_ENDSYSTOLICPRESSURE")
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
    plot.PV.loop(50, 170, 0, 120, new.coords.display, 4, "blue")
    plot.PV.loop(50, 170, 0, 120, FEM.coords.display, 4, "red")



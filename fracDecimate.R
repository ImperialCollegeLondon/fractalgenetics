# Code to interpolate missing values and fit FDs to a set number of slices
# Tim Dawes May 2018
    
    #install.packages("stats")
    library(stats)
    is.this.an.FD.value<- function(vec) {sapply(vec, function(x) {!x %in% c("NA","NaN","Meagre blood pool", "Sparse myocardium", "FD measure failed")})}
    

# Change working directory to point at FD data  
    switch(Sys.info()[['sysname']],
           Windows= {setwd("C:/Users/[your login]/Desktop")},
           Darwin = {setwd("~/Desktop")})
    cat("Current working directory is:",getwd())

# Read in the FD data
    d<- read.csv(file="FD.csv", fill=TRUE)

    
# Pull out the columns to use for interpolation
    FR.all<- d[,10:29]
    rownames(FR.all)<- d$Folder

# Work out how many viable slices are present in each subject (takes about a minute to run)
    no.of.values<- rep(0, nrow(FR.all))
    for (i in 1:nrow(FR.all)) {if ((i/1000)==round(i/1000)) cat(",",i,sep=""); no.of.values[i]<- sum(is.this.an.FD.value(FR.all[i,]))}
    hist(no.of.values, col="blue", xlab="No.of.slices/subject", ylab="Frequency", main="Histogram of number of slices per subject")
    cut.off=3
    abline(v=cut.off, col="red", lty=2, lwd=4)
    
# Remove any subjects with fewer than a set number of FD values available for analysis
# The threshold for this is set by the variable "cut-off"
    
    cat("Excluding ",round(100*length(which(no.of.values<cut.off))/nrow(FR.all),1),"% subjects because they have <",cut.off," slices available.",sep="")
    FR.all<- FR.all[-which(no.of.values<3),]
    
    # Set up the output matrix
      interpNoSlices<- 10
      FRi<- matrix(0, nrow=nrow(FR.all), ncol=interpNoSlices, dimnames=list(rownames(FR.all), paste("Slice_",1:10,sep="")))
      head(FRi)
      
    # Loop over each subject and interpolate to a set number of slices
    # WARNING: loop takes about 55 seconds per 1,000 subjects
  
    for (i in 1:nrow(FR.all))
      {
        if (round(i/100)==(i/100)) {cat(",",i, sep="")}
        values<- is.this.an.FD.value(FR.all[i,])
        o<- which(values==TRUE)
        
        xs.orig<- which(values==TRUE)
        ys.orig<- as.numeric(as.character(unlist(FR.all[i,xs.orig])))
        
        # Fit the data to a template of 10 ventricular levels
            tmp<- ksmooth(xs.orig, ys.orig, kernel="normal", bandwidth=1.5, range.x=range(xs.orig), n.points=interpNoSlices)
            if(sum(is.na(tmp$y))==0) {FRi[i,]<- tmp$y} else {FRi[i,]<- NA}
            
        # Comment this bit in to view a single iteration and how it was dealt with (will slow runtime significantly if left uncommented during loop)
            #par(mfrow=c(2,1))
            #r<- range(na.omit(ys.orig)); a<- "Longitudinal position in RV"; b<- "Fractal dimension"
            #r2<- range(xs.orig)
            #plot(xs.orig, ys.orig, lwd=4, col="red", type="b", ylim=r, xlim=r2, xlab=a, ylab=b, main="Original Data")
            #plot(tmp$x, tmp$y, col="blue", type="b", lwd=4, xlim=r2, ylim=r, xlab=a, ylab=b, main="Interpolated Ten-Slice Model Data")
      }

      
      # Remove any subjects in which there were no points within the hernel width -> Nadaraya-Watson estimator to become 0/0 = NaN
        FRi<- FRi[-which(is.na(FRi[,1])==TRUE),]
      
      
      
# Output results to text file
      write.table(FRi, file="FD_Interpolated.txt", col.names=TRUE, row.names=TRUE)
      

      



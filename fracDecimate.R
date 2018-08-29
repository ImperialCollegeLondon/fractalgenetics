# Code to interpolate missing values and fit FDs to a set number of slices
# Tim Dawes May 2018

fracDecimate<- function (interpNoSlices=10, cut.off=3, data=NULL, filename=NULL,
    id.col.name="Folder", interactive=FALSE) {
    # Install stats package if not currently installed and then load
    if("stats" %in% rownames(installed.packages()) == FALSE) {
        install.packages("stats")
    }
    library(stats)

    # private function
    is.this.an.FD.value <- function(vec) {
        nonFD <- c("NA","NaN","Meagre blood pool", "Sparse myocardium",
            "FD measure failed")
        !(is.na(vec) | vec %in% nonFD)
    }

    # Error handling
    if ((interpNoSlices > 3 && interpNoSlices < 50) == FALSE) {
        stop("\n Interpolation template should be between 3 and 50 slices")
    }
    if (!is.numeric(interpNoSlices)) {
        stop("\n First parameter must be the number of slices to which the data are interpolated")
    }
    if (cut.off < 1 | cut.off > 20) {
        stop("\n Threshold for discarding subjects must be 1-20")
    }
    if (!is.null(filename) && !substr(filename, nchar(filename) - 3, nchar(filename)) == ".csv") {
        stop("\n File must be a .csv file")
    }
    if (is.null(data) && is.null(filename)) {
        stop("Either data or filename must be provided")
    }
    # Read in the FD data and pull out the columns to use for interpolation
    if (is.null(data)) {
        data <- read.csv(file=filename, fill=TRUE, na.strings=c("NaN", "NA"))
    }
    #FR.all<- data[,10:29]
    FR.all <- data[, grepl("Slice \\d{1,2}", colnames(data))]
    rownames(FR.all) <- data[,colnames(data) == id.col.name]

    # Work out how many viable slices are present in each subject
    no.of.values <- colSums(apply(FR.all,1, is.this.an.FD.value))

    if (interactive) {
        hist(no.of.values, col="blue", xlab="No.of.slices/subject",
             ylab="Frequency", main="Histogram of number of slices per subject")
        abline(v=cut.off, col="red", lty=2, lwd=4)
    }

    # Remove any subjects with fewer than a set number of FD values available for analysis
    # The threshold for this is set by the variable "cut-off"
    exclude <- round(100 * length(which(no.of.values < cut.off))/nrow(FR.all),1)
    cat("\n \n Excluding ", exclude, "% subjects because they have <", cut.off,
        " slices available. \n \n", sep="")
    FR.all <- FR.all[which(no.of.values >= cut.off),]

    # Set up the output matrix
    FRi <- matrix(0, nrow=nrow(FR.all), ncol=interpNoSlices,
                  dimnames=list(rownames(FR.all),
                  paste("Slice_", 1:10, sep="")))
    if (interactive) head(FRi)

    # interpolate to a set number of slices per subject
    interpolateSlices <- function(slices, interpNoSlices) {
        xs.orig <- which(is.this.an.FD.value(slices))
        ys.orig <- as.numeric(as.character(slices[xs.orig]))
        tmp <- ksmooth(xs.orig, ys.orig, kernel="normal", bandwidth=1.5,
                       range.x=range(xs.orig), n.points=interpNoSlices)
        if(any(is.na(tmp$y))) {
            return(rep(NA, interpNoSlices))
        } else {
            return(tmp$y)
       }
    }
    cat(" Interpolate slices to ", interpNoSlices, "for each of the",
        nrow(FR.all), "subjects\n")
    FRi <- t(apply(FR.all, 1, interpolateSlices, interpNoSlices=interpNoSlices))
    dimnames(FRi) <- list(rownames(FR.all), paste("Slice_", 1:10, sep=""))

    # Remove any subjects in which there were no points within the kernel width -> Nadaraya-Watson estimator to become 0/0 = NaN
    FRi <- FRi[!is.na(FRi[,1]),]
    return(FRi)
}

#' Convert cartesian coordinates to polar coordinates.
#' 
#' @param x  x-values [vector] of cartesian coordinates [double].
#' @param y  y-values [vector] of cartesian coordinates [double].
#' @return data.frame with r, the radius and theta, the
#' angle of the polar coordinate.
cart2pol <- function(x, y) {
    r <- sqrt(x^2 + y^2)
    r[r==0] <- 1e-4
    theta <- acos(x/r)
    theta[y < 0] <- -theta[y < 0]
    return(data.frame(r=r, theta=theta))
}

#' Convert polar coordinates to cartesian coordinates.
#' 
#' @param r  radius [vector] of polar coordinates [double].
#' @param theta  angle [vector] of polar coordinates in radians [double].
#' @return data.frame with x, the x-values and y, the y-values of the cartesian
#' coordinates.
pol2cart <- function(r, theta) {
    x <- r * cos(theta)
    y <- r * sin(theta)
    return(data.frame(x=x, y=y))
}

#' Interpolate radius of polar coordinates based on closest angle matches.
#' 
#' Conversion of image edges into polar coordinates yields angles depending on
#' the original dimensions of the image i.e. number of pixels in the edge. In 
#' order to make edges directly comparable across different images, this
#' function takes the angles and radii of the image edge and interpolates
#' radius to a specified angle. The interpolation uses the mean value of the 
#' two closest radii. The closest radii are chosen based on closest angles and
#' radii, a step necessary to ensure interpolation for folded edges works 
#' properly. The number of closest angles and radii specified can be set with n.
#' 
#' @param angle_match angle in radians [double] for which to interpolate radius.
#' @param angle_original angles [vector] in radians [double] used for
#' interpolation.
#' @param r_original radii [vector, double] used for interpolation.
#' @param n number of neighbours [n] to check for closest matching radius and
#' angles.
#' @return data.frame with theta, the angle (as specified in angle_match) and r,
#' the corresponding, interpolated radius.
interpolateRadius <- function(angle_match, angle_original, r_original,
                              n=2) {
    distance.angle <- abs(angle_original - angle_match)
    min.angle.pos <- order(distance.angle)
    distance.angle <- distance.angle[min.angle.pos]
    if (distance.angle[1] == 0) {
        return(data.frame(theta=angle_match, r=r_original[min.angle.pos[1]]))
    } else {
        min.angle <- angle_original[min.angle.pos[1]]
        min.r <-  r_original[min.angle.pos[1]]
        distance <- abs(r_original[min.angle.pos[2:n]] - min.r)
        min.r.pos <- order(distance)[1]
        r <- mean(c(min.r, r_original[min.angle.pos[min.r.pos]]))
        return(data.frame(theta=angle_match, r=r))
    }
}

#' Load segmentations of perimeter and edges of trabeculation.
#' 
#' @param IID Sample identifier [character]; has to be the same as identifier
#' of folder containing all images of samples.
#' @param filename path/to/trabeculation file name [character] eg 
#' directory/IID/Edge-Image-Slice-10-ED.png
#' @param sliceInfo [data.frame] containing at least columns IID with the sample
#' identifier and columns specifying the slices with eg Slice10. Uses
#' Slice column information to determine if particular slice in that individual
#' has been segemented and FD estimation was not possible. In the latter case, 
#' for these filenames loadSlices returns NULL.
#' @param nonseg [vector] whose entries [character] specify possible criteria
#' for missing segmentation values in sliceInfo. 
#' @param onlyBackground [logical] if True, only slices containing information
#' about tissue background ie Perimter-On-Slice images are read. In this case,
#' filename has to be of eg directory/IID/Perimter-On-Slice-10-ED.png.
#' @return named list with perimeter, a data.frame containing xvalues, yvalues,
#' intensity, IID and slice information of myocardial outline, and edges, a
#' data.frame containing xvalues, yvalues, intensity, IID and slice information
#' of trabeculation outline. 
loadSlices <- function(IID, filename, sliceInfo,
                       nonseg=c("Sparse myocardium",  "Meagre blood pool", 
                                "FD measure failed"),
                       onlyBackground=FALSE) {
    iid <- gsub(".*/(.*)/(.*)-Slice-.*", '\\1', filename)
    if (IID != iid) stop("Provided ID does not match slice ID")
    slice_id <- gsub(".*/(.*)/(.*)-Slice-(.*)-ED.png", '\\3', filename)
    slice_id <- paste("Slice_", slice_id, sep="")
    slice_value <- sliceInfo[sliceInfo$IID == IID,
                             colnames(sliceInfo) == slice_id]
    if (slice_value  %in% nonseg) return(NULL)
    if (onlyBackground) {
        mask_background <- imager::load.image(filename)
        background <- as.data.frame(mask_background)
        colnames(background) <- c('xvalues', 'yvalues', 'intensity')
        background$IID <- iid
        background$slice <- slice_id
        return(background)
    }
    fileBinary <- gsub("Edge-Image", 'Binary-mask', filename)
    mask_perimeter <- imager::load.image(fileBinary)
    perimeter <- as.data.frame(imager::boundary(mask_perimeter))[,-3]
    colnames(perimeter) <- c('xvalues', 'yvalues', 'intensity')
    perimeter$IID <- iid
    perimeter$slice <- slice_id

    mask_edges <- imager::load.image(filename)
    edges <- as.data.frame(imager::boundary(mask_edges))[,-3]
    colnames(edges) <- c('xvalues', 'yvalues', 'intensity')
    edges$IID <- iid
    edges$slice <- slice_id
    return(list(perimeter=perimeter, edges=edges))
}

#' Prepare images for registration.
#' 
#' Images loaded via loadSlices are mean centered and scale to maximum image
#' dimension of 1. Based on the mean centered and scaled images, the xvalues
#' and yvalue are converted to polar coordinates via cart2pol.
#' 
#' @param sliceData Output of loadSlices (with onlyBackground==FALSE) [list] ie
#' list with perimeter, a data.frame containing xvalues, yvalues,
#' intensity, IID and slice information of myocardial outline, and edges, a
#' data.frame containing xvalues, yvalues, intensity, IID and slice information
#' of trabeculation outline.
#' @return Named list with perimeter, edges and centering. Perimeter is a
#' data.frame containing centered and scaled xvalues, centered and scaled
#' yvalues, radius (r) and angle (theta) of xvalues and yvalues converted to
#' polar corrdinates, pixel intensity, IID and slice information of myocardial
#' outline.  Edges is a data.frame containing centered and scaled xvalues,
#' centered and scaled yvalues, radius (r) and angle (theta) of xvalues and
#' yvalues converted to polar corrdinates, pixel intensity, IID and slice
#' information of trabeculation outline. Centering is a list containing the
#' center coordinates (xcenter and ycenter) as well as the maximum pixel
#' distance (maxdist) used for centering and scaling of edges and perimeter.
prepOutlines <- function(sliceData) {
    perimeter <- sliceData$perimeter
    ycenter <- mean(perimeter$yvalues)
    xcenter <- mean(perimeter$xvalues)
    perimeter$yvalues  <- perimeter$yvalues - ycenter
    perimeter$xvalues  <- perimeter$xvalues - xcenter
    
    distance <- sqrt(perimeter$xvalues^2 + perimeter$yvalues^2)
    maxdist <- max(distance)
    perimeter$yvalues <- perimeter$yvalues/maxdist
    perimeter$xvalues <- perimeter$xvalues/maxdist
    
    polarPerimeter <- cart2pol(perimeter$xvalues, perimeter$yvalues)
    perimeter$r <- polarPerimeter$r
    perimeter$theta <- polarPerimeter$theta
    
    edges <- sliceData$edges
    edges$yvalues  <- edges$yvalues - ycenter
    edges$xvalues  <- edges$xvalues - xcenter
    edges$yvalues <- edges$yvalues/maxdist
    edges$xvalues <- edges$xvalues/maxdist
    
    polarEdges <- cart2pol(edges$xvalues, edges$yvalues)
    edges$r <- polarEdges$r
    edges$theta <- polarEdges$theta
    return(list(perimeter=perimeter, edges=edges,
                centering=list(xcenter=xcenter, ycenter=ycenter,
                               maxdist=maxdist)))
}

#' Find center of non-LV tissue as reference point for registration.
#' 
#' Loads all background images (ie Perimeter-on-Slice) for the individual (IID)
#' and finds the image that contains most non-LV tissue pixels (as specified by
#' intensity). From this image, it finds the center of non-LV tissue pixels and
#' converts it to a polar coordinate. The angle of this coordinate is
#' substracted from the desired reference center and returned as the rotational
#' angle (bg.rotation).
#' 
#' @param IID Sample identifier [character]; has to be the same as identifier
#' of folder containing all images of samples.
#' @param directory path/to/direcory  [character] containing the IID folder with
#' all background images for this individual eg
#' directory/IID/Perimeter-On-Slice-10-ED.png.
#' @param sliceInfo [data.frame] containing at least columns IID with the sample
#' identifier and columns specifying the slices with eg Slice10. Uses
#' Slice column information to determine if particular slice in that individual
#' has been segemented and FD estimation was not possible. In the latter case, 
#' for these filenames loadSlices returns NULL.
#' @param intensity Pixel intensity value [double] of the non-LV tissue; should
#' be precise to at least four decimal places.
#' @param nonseg [vector] whose entries [character] specify possible criteria
#' for missing segmentation values in sliceInfo. 
#' @param center Reference center angle in radians [double] for rotational
#' registration.
#' @return Named list with bg.rotation, the angle in radians for rotational 
#' registration and bg_max, a data.frame with xvalues, yvalues, intensity, IID
#' and slice information of background image.
findOrientation <- function(IID, directory, sliceInfo, intensity, center=pi/2,
                            nonseg=c("Sparse myocardium",  "Meagre blood pool", 
                                     "FD measure failed")) {
    d <- file.path(directory, IID)
    files <- list.files(d, pattern="Perimeter-On-Segmentation",
               full.names=TRUE)
    all_bg <- lapply(files, function(x){
        bg <- loadSlices(IID=IID, filename=x, sliceInfo=sliceInfo,
                         nonseg=nonseg, onlyBackground=TRUE)
        values <- as.data.frame(table(bg$intensity), stringsAsFactors=FALSE)
        values$Var1 <- as.numeric(values$Var1)
        if (any(round(values$Var1,4)==round(intensity,4))) {
            bg_count <- values$Freq[round(values$Var1,4)==round(intensity,4)]
        } else {
            bg_count <- 0
        }
        return(list(bg=bg, bg_count=bg_count))
    })
    all_counts <- sapply(all_bg, function(x) x$bg_count)
    bg_max <- all_bg[[which.max(all_counts)]]$bg
    bg_pixels <- bg_max[round(bg_max$intensity,4) == round(intensity,4),]
    
    polarBg <- cart2pol(bg_pixels$xvalues, bg_pixels$yvalues)
    bg.rotation <- center - mean(polarBg$theta)
    return(list(bg.rotation=bg.rotation, bg_max=bg_max))
}

#' Rotational and radial registration perimeters and edges.
#' 
#' Perimeter and edges dataframes prepared for registration, ie mean centered,
#' scaled and transformed to polar coordinates (via prepOutlines) are registered
#' to circular outline. This is done in two steps. Firstly, both perimeter and
#' edge coordinates angular coordinates are aligned (via interpolateRadius).
#' Secondly, the transformation of the perimeter coordinates to the closest
#' enclosing circle coordinates is computed. These coordinate-wise
#' transformations are then applied to the edge coordinates. If center is
#' specified, the radially registed edges are the rotated according to
#' rotational center angle (in radians).
#' 
#' @param  perimeter [data.frame] containing centered and scaled (xvalues),
#' centered and scaled yvalues (yvalues), radius (r) and angle (theta) of
#' xvalues and yvalues converted to polar corrdinates, pixel intensity
#' (intensity), individual indentifier (IID) and slice identifier (slice) of
#' myocardial outline.
#' @param edges [data.frame] containing centered and scaled xvalues (xvalues),
#' centered and scaled yvalues (yvalues), radius (r) and angle (theta) of
#' xvalues and yvalues converted to polar corrdinates, pixel intensity
#' (intensity), individual indentifier (IID) and slice identifier (slice) of
#' trabeculation outline. 
#' @param angles [vector] of angles in radians [double] used as coordinates for
#' registration.
#' @param targetRadius Radius [double] of circle to register to.
#' @param n number of neighbours [n] to check for closest matching radius and
#' angles.
#' @param rotation Angle in radians [double] for rotational registration.
#' @return Named list with transformed perimeter and edges. perimeter is a
#' [data.frame] containing centered and scaled (xvalues),
#' centered and scaled yvalues (yvalues), radius (r) and angle (theta) of
#' xvalues and yvalues converted to polar coordinates, transformed radius
#' (r.transformed) and rotated angle (theta.rotated) of the registered
#' coordinates, the pixel intensity (intensity), individual indentifier (IID)
#' and slice identifier (slice) of myocardial outline. edges is a [data.frame]
#' containing centered and scaled (xvalues), centered and scaled yvalues
#' (yvalues), radius (r) and angle (theta) of
#' xvalues and yvalues converted to polar coordinates, transformed radius
#' (r.transformed) and rotated angle (theta.rotated) of the registered
#' coordinates, the pixel intensity (intensity), individual indentifier (IID)
#' and slice identifier (slice) of the trabeculation outline. 
registerSlices <- function(perimeter, edges, angles=seq(-pi, pi, pi/1000),
                           targetRadius=1, n=2, rotation=NULL) {
    edges.transform <- do.call(rbind,
                               lapply(angles, interpolateRadius,
                                      angle_original=edges$theta,
                                      r_original=edges$r,
                                      n=n))
    edges.transform$IID <- unique(edges$IID)
    edges.transform$slice <- unique(edges$slice)
    
    perimeter.transform <- do.call(rbind,
                                   lapply(angles, interpolateRadius,
                                          angle_original=perimeter$theta,
                                          r_original=perimeter$r,
                                          n=n))
    perimeter.transform$IID <- unique(perimeter$IID)
    perimeter.transform$slice <- unique(perimeter$slice)
    transformVector <- targetRadius - perimeter.transform$r
    perimeter.transform$r.transformed <- targetRadius

    edges.transform$r.transformed <- edges.transform$r + transformVector

    if (!is.null(rotation)) {
        perimeter.transform$theta.rotated <- perimeter.transform$theta - rotation
        edges.transform$theta.rotated <- edges.transform$theta - rotation
    }
    
    return(list(perimeter=perimeter.transform, 
                edges=edges.transform))
}

#' Load, align and register all segmented slices per individual.
#' 
#' @param IID Sample identifier [character]; has to be the same as identifier
#' of folder containing all images of samples.
#' @param directory path/to/direcory  [character] containing the IID folder with
#' all background images for this individual eg
#' directory/IID/Perimeter-On-Slice-10-ED.png.
#' @param sliceInfo [data.frame] containing at least columns IID with the sample
#' identifier and columns specifying the slices with eg Slice10. Uses
#' Slice column information to determine if particular slice in that individual
#' has been segemented and FD estimation was not possible. In the latter case, 
#' for these filenames loadSlices returns NULL.
#' @param intensity Pixel intensity value [double] of the non-LV tissue; should
#' be precise to at least four decimal places.
#' @param angles [vector] of angles in radians [double] used as coordinates for
#' registration.
#' @param targetRadius Radius [double] of circle to register to.
#' @param n number of neighbours [n] to check for closest matching radius and
#' angles.
#' @param nonseg [vector] whose entries [character] specify possible criteria
#' for missing segmentation values in sliceInfo. 
#' @param center Reference center angle in radians [double] for rotational
#' registration.
#' @return List with list entry for each segmented slice. Each list entry
#' contains a named list with the output of registerSlices i.e. a list with
#' transformed perimeter and edges. perimeter is a [data.frame] containing
#' centered and scaled (xvalues),
#' centered and scaled yvalues (yvalues), radius (r) and angle (theta) of
#' xvalues and yvalues converted to polar coordinates, transformed radius
#' (r.transformed) and rotated angle (theta.rotated) of the registered
#' coordinates, the pixel intensity (intensity), individual indentifier (IID)
#' and slice identifier (slice) of myocardial outline. edges is a [data.frame]
#' containing centered and scaled (xvalues), centered and scaled yvalues
#' (yvalues), radius (r) and angle (theta) of
#' xvalues and yvalues converted to polar coordinates, transformed radius
#' (r.transformed) and rotated angle (theta.rotated) of the registered
#' coordinates, the pixel intensity (intensity), individual indentifier (IID)
#' and slice identifier (slice) of the trabeculation outline. 
processSlices <- function(IID, directory, sliceInfo, intensity=0.941176,
                          angles=seq(-pi, pi, pi/1000),
                          targetRadius=1, n=2,
                          nonseg=c("Sparse myocardium",  "Meagre blood pool", 
                                   "FD measure failed"), 
                          center=pi/2, rotation=FALSE, verbose=TRUE) {
    if (verbose) message("Processing file from ", IID, "\n")
    if (rotation) {
        background <- findOrientation(IID=IID, directory=directory,
                                      sliceInfo=sliceInfo, intensity=intensity,
                                      nonseg=nonseg, center=center)
        rc=background$bg.rotation
    } else {
        rc=NULL
    }
    filelist <- list.files(file.path(directory, IID),
                           pattern="Edge-Image-Slice-.*-ED.png",
                           full.names=TRUE)
    if (length(filelist) == 0) {
        stop("No edge images in specified directory")
    }
    results <- lapply(filelist, function(x, iid, rc) {
        images <- loadSlices(IID=iid, filename=x, sliceInfo=sliceInfo,
                             nonseg=nonseg)
        
        if (!is.null(images)) {
            outlines <- prepOutlines(images)
            outlines.registered <- registerSlices(perimeter=outlines$perimeter,
                                                  edges=outlines$edges,
                                                  angles=angles,
                                                  targetRadius=targetRadius,
                                                  n=n, rotation=rc)
            return(list(original=outlines, registered=outlines.registered))
        } else {
            return(NULL)
        }
    }, iid=IID, rc=rc)
    results <- results[sapply(results, function(x) !is.null(x))]
    return(results)
}
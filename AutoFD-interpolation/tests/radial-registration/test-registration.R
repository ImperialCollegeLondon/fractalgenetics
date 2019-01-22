## Example of running radial registration on a
## set of target images for testIID. 
## Target images needed are images containing non-LV tissue
## to determine center of background pixels for rotation, 
## binary mask images of myocard region of interest (for outline),
## and trabeculation edges. 
## 
## The different steps of the analyses are visualised
## as radial plots of the original, transformed,
## rotated and rotated+tranformed images.
## 
## Test data can be found in AutoFD_interpolation/tests/radial-registration/data

###############################
## libraries and functions ####
###############################

options(import.path=".")
modules::import_package('plotrix')
rr <- modules::import('radial-registration')

############
## data ####
############

# Run from or setwd to analysis directory for paths to work
directory <- 'tests/radial-registration/data'
ID <- 'testIID'
bg_intensity <-  0.9411765
nonseg <-  c("Sparse myocardium", "Meagre blood pool",  "FD measure failed")

# FD measurements or failed segmentations for testIID
slices <- data.table::fread('tests/radial-registration/data/sliceInfo.csv')

################
## analysis ####
################

# Get image with maximum amount of background pixels
background <- rr$findOrientation(IID=ID, directory=directory,
                                 sliceInfo=slices, intensity=bg_intensity)
background.matrix <- reshape2::acast(background$bg_max, xvalues ~ yvalues,
                                     value.var = 'intensity')
background.matrix <- t(apply(background.matrix, 1, rev))
max_bg_slice <- unique(background$bg_max$slice)

# Process all valid, segmented slies of testIID
processed <- rr$processSlices(IID=ID, directory=directory, sliceInfo=slices,
                           n=2)
names(processed) <- lapply(processed, function(x)
    unique(x$registered$perimeter$slice))

# Find slice corresponing to maximum background slice (for direct comparison of
# rotation)
pos <- which(names(processed) == max_bg_slice)
perimeter.reg <- processed[[pos]]$registered$perimeter
perimeter <- processed[[pos]]$original$perimeter
edges.reg <- processed[[pos]]$registered$edges
edges <- processed[[pos]]$original$edges

# Depict image of background slice and myocardial and trabeculation outlines
# before and after step registration, including intermediate steps.
layout(matrix(c(1, 2, 3, 4, 5, 1, 6, 7, 8, 9), 2, 5, byrow = TRUE))
image(background.matrix)
radial.plot(perimeter.reg$r, radial.pos=perimeter.reg$theta, labels="",
            rp.type = "s", point.symbols = 20, show.grid.labels = FALSE,
            radial.lim=c(0, 1), point.col='#99d8c9', start=pi/2)
radial.plot(perimeter.reg$r,
            radial.pos=perimeter.reg$theta.rotated, labels="",
            rp.type = "s", point.symbols = 20, start=pi/2,
            radial.lim=c(0, 1), point.col='#66c2a4',
            show.grid.labels = FALSE)
radial.plot(perimeter.reg$r.transformed,
            radial.pos=perimeter.reg$theta, labels="",
            rp.type = "s", point.symbols = 20, start=pi/2,
            radial.lim=c(0, 1), point.col='#2ca25f',
            show.grid.labels = FALSE)
radial.plot(perimeter.reg$r.transformed,
            radial.pos=perimeter.reg$theta.rotated, labels="",
            rp.type = "s", point.symbols = 20, start=pi/2,
            radial.lim=c(0, 1), point.col='#006d2c',
            show.grid.labels = FALSE)

radial.plot(edges.reg$r, radial.pos=edges.reg$theta, labels="",
            rp.type = "s", point.symbols = 20, show.grid.labels = FALSE,
            radial.lim=c(0, 1), point.col='#9ebcda', main = 'original',
            start=pi/2)
radial.plot(edges.reg$r,
            radial.pos=edges.reg$theta.rotated, labels="", rp.type = "s",
            radial.lim=c(0, 1), point.symbols = 20, show.grid.labels = FALSE,
            point.col='#8c96c6', main = 'rotated', start=pi/2 )
radial.plot(edges.reg$r.transformed,
            radial.pos=edges.reg$theta, labels="", rp.type = "s",
            radial.lim=c(0, 1), point.symbols = 20, show.grid.labels = FALSE,
            point.col='#8856a7', main = 'transformed', start=pi/2)
radial.plot(edges.reg$r.transformed,
            radial.pos=edges.reg$theta.rotated, labels="", rp.type = "s",
            radial.lim=c(0, 1), point.symbols = 20, show.grid.labels = FALSE,
            point.col='#810f7c', main = 'rotated and transformed', start=pi/2)




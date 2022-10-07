# FD processing
R code for post-processing of FD values output by
[AutoFD](https://github.com/ImperialCollegeLondon/fractalgenetics/tree/master/automated-fractal-analysis).

1. interpolate.R (fracDecimate function)
    Interpolation of FD values to a common number of slices across individuals.
1. summaryFD.R (summaryStatistics function)
    Summary statistics of FD values per individual: mean global FD, mean and max
    Apical FD, mean and max Basal FD. This is an R version of the matlab
    function [FDstatististics](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/automated-fractal-analysis/pft_JC_FDStatistics.m).
    The function has been tested for equivalent output in
    [tests](https://github.com/ImperialCollegeLondon/fractalgenetics/tree/master/fractal-analysis-processing/tests).
1. radial-registration.R
    Collection of functions for co-registration of myocardial and trabeculation
    outline to closest enclosing circle. Functions include transformation of
    images to polar coordinates and interpolation of image data in polar space.

## Installation
 fracDecimate.R, summaryFD.R and radial-registration.R can simply be sourced to access
 the relevant functions. Alternatively,
 the [modules package](https://github.com/klmr/modules) offers flexible and tidy
 integration of R source files. Application of both functions via modules in
 [ukbb-fd](UK-Biobank/phenotypes/preparePheno.r).

A test file `FD.csv` is included in the repo, but any csv file can be called by
the function.

## Running the code
### fracDecimate
    output <- fracDecimate(interpNoSlices, cut.off, data, filename, id.col.name,
    interactive, verbose) [Enter]
with:
interpNoSlices = A number between 3 and 50 indicating the number of slices
    you want the data interpolated to (default: 10)

cut.off = A number between 1 and 20 indicating the minimum number of
    datapoints per patient for that patient to be included in the analysis
    (default: 3)

filename = The filename of the FD data exported from AutoFD (default: FD.csv)
    or
data = [N x MaxNrSlices] data matrix with FD measurements per slice for all
    MaxNrSlices across N individuals. Slice columns will be
    automatically extracted and need to contain "Slice" in their name.

id.col.name = column name of sample id column. If sample IDs are in input
    rownames specify id.col.names = "rownames" (default: Folder)

interactive = If TRUE, histogram of FD is plotted (default: FALSE)

verbose = If TRUE, progress messages are printed (default: TRUE)

usage is e.g., `output <- fracDecimate(interpNoSlices=9, filename="FD.csv")`

The output is a [N x interpNoSlices] matrix with the interpolated data. Rownames
are the subject IDs. Column names are: "Slice_x" where x is a number from 1 to
interpNoSlices.

Notes to fracDecimate
- Code will assume that any FD values which are not "NA", "NaN",
  "Meagre blood pool", "Sparse myocardium" or "FD measure failed" are FD values.
- Fitting is by a kernel regression estimate with a bandwidth of 1.5 slices

### summaryFD.R
    summaryStatistics(data, discard, NaN.val)
with:
data = [N x NrSlices] matrix of FD value for N individuals and NrSlices

discard = if TRUE, first and last slice are excluded for computing summary
    statistics (default: FALSE).

NaN.val = vector with character strings of accepted NaN values (default:
    c("Meagre blood pool", "Sparse myocardium" or "FD measure failed"))

The output is a [N x 6] matrix with summary statistics for N individuals:
NrSlices used for computation, mean global FD, mean and max Apical FD, mean and
max Basal FD

usage is e.g., `summaryStatistics_output <- t(apply(output, 1, summaryStatistics))`

### radial-registration.R
A detailed example of how to run radial co-registration of myocard and
trabeculation outline based on intermediate files generated in
[AutoFD](https://github.com/ImperialCollegeLondon/fractalgenetics/tree/master/automated-fractal-analysis) can be found
[here](https://github.com/ImperialCollegeLondon/fractalgenetics/blob/master/fractal-analysis-processing/tests/radial-registration/test-registration.R).

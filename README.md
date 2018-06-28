# AutoFD_interpolation
R code for post-processing of FD values output by [AutoFD](https://github.com/UK-Digital-Heart-Project/AutoFD).

## Requisites

R code will install the 'stats' package, if it is not currently installed

## Installation
Open fracDecimate.R 

A test file `FD.csv` is included in the repo, but any csv file can be called by the function.

## Running the code
1. Run the code, which defines a function called `fracDecimate`

2. In the console type: output.txt<- fracDecimate(a,b,c) [Enter], where:

    a = A number between 3 and 50 indicating the number of slices you want the data interpolated to (suggest: 10)

    b = A number between 1 and 20 indicating the minimum number of datapoints per patient for that patient to be included in the analysis (suggest: 3)

    c = The filename of the FD data exported from AutoFD (suggest: FD.csv)
    
    eg `output.txt<- fracDecimate(10,3,"FD.csv")`

The output is a text file with the interpolated data. Rownames are the subject IDs. Column names are: "Slice_x" where x is a number from 1 to a.

## Notes
- Code will assume that any FD values which are not "NA", "NaN", "Meagre blood pool", "Sparse myocardium" or "FD measure failed" are FD values.
- Interpolation takes ~55seconds per 1,000 subjects (Intel Xeon CPU 2.40GHz, 32Gb RAM, R 3.4.0, RStudio 1.1.447)
- Function includes code (commented out by default) to draw a single example of the original data and the interpolated fit. To understand the fitting process better this can be un-commented. Suggest running only a single subject rather than the whole loop.
- Fitting is by a kernel regression estimate with a bandwidth of 1.5 slices




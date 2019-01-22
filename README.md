# AutoFD
Matlab script for automated fractal analysis of segmented cardiac images.

## Prerequisites
- Matlab.
- Windows or Linux.
- GhostScript.
- Possibly XPDF.

## Input data

Compressed NIfTI files (NII.GZ) of end-diastolic left ventricular short axis stacks.

## Method

The approach is based on our [FracAnalyse](https://github.com/UK-Digital-Heart-Project/fracAnalyse) software, 
but uses pre-exisiting image segmentations to determine a region of interest within the myocardium for fractal analysis.  

![FD images](https://github.com/UK-Digital-Heart-Project/AutoFD/blob/master/FDworkflow.png)

## Installation
Clone this repo to a folder in your MATLAB workspace then add all directories to the path:

```addpath(genpath('folder')); savepath;```

## Usage
Put the input data into a top-level folder with sub-folders for each subject containing:
  * Grayscale main image  ```sa_ED.nii.gz```
  * Segmentation ```seg_sa_ED.nii.gz``` or ```seg_sa_ED.gipl```.

The labels are Background  = 0, Blood Pool  = 1, Myocardium  = 2, Other = 3 or 4.

Run the script ```pft_FractalDimensionCalculationOnMultipleFolders```

Dialogue boxes are used to confirm the following: 

  * Stacks are stored from Apex to Base (default for the Test Data) or vice-versa.
  * Interpolate images x4 in each direction in-plane (default for the Test Data), or to 0.25 mm pixels.
  * Minimum blood pool pixel count (default = 64 for the sample Test Data, optimised by means of a balanced probability calculation comparing manual and automated workflows) and percentage of blood pool connected to myocardium (default = 50%). Refer to the [Processing Flowchart](https://github.com/UK-Digital-Heart-Project/AutoFD/blob/master/Processing%20Flowchart.pdf) for details.
  * Keep or discard the end slices in the calculation of the summary statistics (the default is keep).
  
## Test data
Data for 7 subjects (A to G) are held in [Test Data](https://github.com/UK-Digital-Heart-Project/AutoFD/blob/master/Test%20data.zip). Unzip and select as the top level folder for analysis. The acquisition order must be Apex to Base.
The preferred interpolation is x4. For the other parameters, either accept or edit the processing defaults.

## Outputs
Outputs are written to a new sub-folder ```Automated FD Calculation Results```.  Each subject's folder will contain intermediate images (see figure) and box-counting statistics.

Fractal dimension values are output to ```Summary-Auto-FD-v0.csv```. If you run the script more than once, new results will be appended.

## Citation

Tokarczuk P, Dawes TJW, O'Regan DP (2018, April 25). UK-Digital-Heart-Project/AutoFD: AutoFD (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.1228963

[![DOI](https://zenodo.org/badge/129865649.svg)](https://zenodo.org/badge/latestdoi/129865649)



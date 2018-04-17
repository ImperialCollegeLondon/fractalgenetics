# AutoFD
Matlab script for automated fractal analysis of segmented cardiac images.

## Prerequisites
- Matlab.
- Windows or Linux.
- GhostScript.

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
Put the input data into a top-level folder, containing one or more sub-folders, each containing:
  * Grayscale main image  ```sa_ED.nii.gz```
  * Segmentation ```seg_sa_ED.nii.gz``` or ```seg_sa_ED.gipl```.

The labels are Background  = 0, Blood Pool  = 1, Myocardium  = 2, Other = 3 or 4.

Run the script ```pft_FractalDimensionCalculationOnMultipleFolders```

Dialogue boxes are used to confirm the following: 

  * Stacks are stored from Base to Apex (default) or vice-versa.
  * Interpolate images to 0.25 mm pixels (default) or x4 in each direction in-plane.
  * Minimum blood pool pixel count (default = 50) and a connection percentage (of the blood pool perimeter
  to the myocardium). Refer to the [Processing Flowchart](https://github.com/UK-Digital-Heart-Project/AutoFD/blob/master/Processing%20Flowchart.pdf) for details.
  * Keep or discard the end slices in the calculation of the summary statistics (the default is Keep).

## Test data
Data for 7 subjects (A to G) is held in [Test Data](https://github.com/UK-Digital-Heart-Project/AutoFD/tree/master/Test%20data). Unzip  and select as the top level folder for analysis. Choose acquisition order of Apex to Base.

## Outputs
Outputs are written to a new sub-folder ```Automated FD Calculation Results```.  Each subject's folder will show the intermediate images (see figure) and box-counting.

Fractal dimensional values are output to ```Summary-Auto-FD-v0.csv```. If you run the script more than once, new results will be appended.

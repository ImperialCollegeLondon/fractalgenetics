# AutoFD
Matlab script for automated fractal analysis of segmented cardiac images.

## Prerequisites
- Matlab.
- Windows or Linux.
- GhostScript

## Input data

Compressed NIfTI files (NII.GZ) of end-diastolic left ventricular short axis stacks.

## Method

The approach is based on our [FracAnalyse](https://github.com/UK-Digital-Heart-Project/fracAnalyse) software, 
but uses pre-exisiting image segmentations to determine a region of interest within the myocardium for fractal analysis.  

![FD images](https://github.com/UK-Digital-Heart-Project/AutoFD/blob/master/FD_workflow.png)

## Installation
Clone this repo to a folder in your MATLAB workspace then add all directories to the path:

```addpath(genpath('folder')); savepath;```

## Usage
Put the input data into a top-level folder, containing one or more sub-folders, each containing:
  * Grayscale main image  ```sa_ED.nii.gz```
  * Segmentation ```seg_sa_ED.nii.gz``` or ```seg_sa_ED.gipl```.

The labels are Background  = 0, Blood Pool  = 1, Myocardium  = 2, Other Heart = 3 or 4.

Run;
```pft_FractalDimensionCalculationOnMultipleFolders```;

There are dialog boxes for 

  * Stacks are stored from Base to Apex (the default) or vice-versa.
  * Interpolate images to 0.25 mm pixels (the default) or x4 in each direction in-plane.
  * Minimum blood pool pixel count (the default is 50) and a connection percentage (of the blood pool perimeter
  to the myocardium): the default is 50.0. Refer to the Processing Flowchart for details.
  * Keep or discard the end slices in the calculation of the summary statistics (the default is Keep =
  Do Not Discard).

## Test data
- See the archive in Examples; both normal execution and most of the error conditions in the flowchart are illustrated.

## Outputs
- A new sub-folder will be created under the original top-level folder, with results grouped by study underneath.
- Each study will have its own audit trail, consisting of images documenting the process of the FD calculation,
  and some small CSV files, slice by slice.
- A master CSV file will be created in the original top-level folder, as well as a backup.
  If you run the script more than once, new results will be appended to the master file.

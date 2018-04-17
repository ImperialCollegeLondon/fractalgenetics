# AutoFD
Matlab script for automated fractal analysis of segmented cardiac images.

## Prerequisites
- Matlab.
- Windows or Linux.
- A working installation of GhostScript on either platform.

## Input data

Compressed NIfTI files (NII.GZ).

## Method

The approach is based on our [FracAnalyse](https://github.com/UK-Digital-Heart-Project/fracAnalyse) software, 
but uses pre-exisiting image segmentations to determine a region of interest within the myocardium for fractal analysis.  

![FD images](https://github.com/UK-Digital-Heart-Project/AutoFD/blob/master/FD_workflow.png)

## Installation
- Unpack the repo to a folder in your MATLAB workspace.
- Add the main directory (containing pft_FractalDimensionCalculationOnMultipleFolders.m) and its sub-folders to the path.
- Be sure to know where you have your installation of GhostScript, in case you are prompted for it (this should only happen once).

## Usage
- Organize your input data into a top-level folder, containing one or more sub-folders, each containing:
  - A grayscale image stack called "sa_ED.nii.gz" (short-axis, end-diastole);
  - A segmentation stack called "seg_sa_ED.nii.gz" or "seg_sa_ED.gipl".
    This latter may be padded in the slice direction with zeros, but the script will squeeze those out so that the
    grayscale image and the segmentation match.
    The labelling scheme is:
    - Background  = 0;
    - Blood Pool  = 1;
    - Myocardium  = 2;
    - Other Heart = 3 or 4.
- Check whether your stacks are stored from Base to Apex (the default) or vice-versa.
- Decide whether to interpolate your images to 0.25 mm pixels (the default) or x4 in each direction in-plane.
- Decide on a minimum blood pool pixel count (the default is 50) and a connection percentage (of the blood pool perimeter
  to the myocardium): the default is 50.0. Refer to the Processing Flowchart for details.
- Decide whether to keep or discard the end slices in the calculation of the summary statistics (the default is Keep =
  Do Not Discard).

## Test data
- See the archive in Examples; both normal execution and most of the error conditions in the flowchart are illustrated.

## Outputs
- A new sub-folder will be created under the original top-level folder, with results grouped by study underneath.
- Each study will have its own audit trail, consisting of images documenting the process of the FD calculation,
  and some small CSV files, slice by slice.
- A master CSV file will be created in the original top-level folder, as well as a backup.
  If you run the script more than once, new results will be appended to the master file.

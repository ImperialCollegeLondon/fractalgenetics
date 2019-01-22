%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pft_FractalDimensionCalculationOnMultipleFolders                                                                              %
%                                                                                                                               %
% A function to process all the sub-folders within a top-level folder.                                                          %
%                                                                                                                               %
% PFT - 15. 11. 2016.                                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear the workspace as usual
clear all
close all
clc

fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Locate a batch folder
if ispc
  Username = getenv('Username');
  Home = fullfile('C:', 'Users', Username, 'Desktop');
elseif isunix
  [ Status, CmdOut ] = system('whoami');
  Home = fullfile('home', CmdOut, 'Desktop');
end  
  
TopLevelFolder = uigetdir(Home, 'Select a top-level folder with data folders inside');

if (TopLevelFolder == 0)
  msgbox('No selection - quitting !', 'Exit');
  return;
end

% Locate all the folders beneath
Listing = dir(TopLevelFolder);
Entries = { Listing.name  };
Folders = [ Listing.isdir ];
Entries = Entries(Folders);

SingleDot = strcmp(Entries, '.');
Entries(SingleDot) = [];
DoubleDot = strcmp(Entries, '..');
Entries(DoubleDot) = [];

Results = strcmp(Entries, 'Automated FD Calculation Results - x4');

if ~isempty(Results)
  Entries(Results) = [];
end

Results = strcmp(Entries, 'Automated FD Calculation Results - 0.25 mm pixels');

if ~isempty(Results)
  Entries(Results) = [];
end

Entries = Entries';

if isempty(Entries)
  msgbox('No sub-folders found.', 'Exit');
  return;
end

NDIRS = length(Entries);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fetch the acquisition slice order (as implemented in the segmentation files)
  AcquisitionOrder = pft_GetAcquisitionOrder;
% AcquisitionOrder = 'Base to Apex';

% Select the type of interpolation
  InterpolationType = pft_GetInterpolationType;
% InterpolationType = 'Imresize - 0.25 mm pixels - cubic';

% Set the default perimeter type - there is no choice here, and if the default cannot be created, then the brute-force Ansatz is applied
  PerimeterType = 'Out from blood pool';

% Fetch the blood pool threshold parameters - 60-65 (64) pixels is optimum for Genscan, 38 for UKBB, according to TJWD's balanced/maximum probability study
  [ MinimumPixelCount, ConnectedPercentage ] = pft_GetBloodPoolThresholdParameters;
% MinimumPixelCount = 38;
% ConnectedPercentage = 50.0;

% Ask whether to trim data for summary FD statistics
  Ans = questdlg('Discard end slices for summary statistics ?', 'Processing decision', 'Yes', 'No', 'No');
% Ans = 'No';

switch Ans
  case { '', 'No' }
    DiscardEndSlices = 'No';
  case 'Yes'
    DiscardEndSlices = 'Yes';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select the output Excel sheet and back it up straightaway
SummaryFile       = fullfile(TopLevelFolder, 'Summary-Auto-FD-v0.csv');
SummaryBackupFile = fullfile(TopLevelFolder, 'Summary-Auto-FD-v0-Backup.csv');

if (exist(SummaryFile, 'file') ~= 2)
  Head = [ 'Folder,', ...
           'Slice order in segmentation,', 'Interpolation,', 'Default perimeter drawn,', ...
           'BP threshold (pixels),', 'Connection threshold (per cent),', ...
           'Original resolution / mm,', 'Output resolution / mm,', ...
           'Slices present,', ...
           'FD - Slice 1,', 'Slice 2,', 'Slice 3,', 'Slice 4,', 'Slice 5,', ...
           'Slice 6,', 'Slice 7,', 'Slice 8,', 'Slice 9,', 'Slice 10,', ...
           'Slice 11,', 'Slice 12,', 'Slice 13,', 'Slice 14,', 'Slice 15,', ...
           'Slice 16,', 'Slice 17,', 'Slice 18,', 'Slice 19,', 'Slice 20,', ...
           'End slices discarded for statistics,', ...
           'Slices evaluated,', 'Slices used,', ...
           'Mean global FD,', ...
           'Mean basal FD,', 'Mean apical FD,', ...
           'Max. basal FD,', 'Max. apical FD' ];     
                 
  fid = fopen(SummaryFile, 'at');
  fprintf(fid, '%s\n', Head);
  fclose(fid);
end

copyfile(SummaryFile, SummaryBackupFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are signalling error conditions for the o/p CSV file
MeagreBloodPool  = -111;
SparseMyocardium = -222;
NoROICreated     = -333;
FDMeasureFailed  =  0.0;    % Signal that an attempt was made, but failed - this will be excluded from the FD statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process all the suitable folders
switch InterpolationType    
  case 'Imresize - (x4 x4) - cubic'
    if (exist(fullfile(TopLevelFolder, 'Automated FD Calculation Results - x4'), 'dir') ~= 7)
      mkdir(TopLevelFolder, 'Automated FD Calculation Results - x4');
    end
  case 'Imresize - 0.25 mm pixels - cubic'
    if (exist(fullfile(TopLevelFolder, 'Automated FD Calculation Results - 0.25 mm pixels'), 'dir') ~= 7)
      mkdir(TopLevelFolder, 'Automated FD Calculation Results - 0.25 mm pixels');
    end
end   

h1 = waitbar(0, sprintf('Processed 0 of %1d folders', NDIRS), 'Units', 'normalized', 'Position', [0.225 0.45 0.2 0.1]);

set(h1, 'Name', 'Progress - folders');

for n = 1:NDIRS
    
  switch InterpolationType    
    case 'Imresize - (x4 x4) - cubic'
      if (exist(fullfile(TopLevelFolder, 'Automated FD Calculation Results - x4', Entries{n}), 'dir') ~= 7)
        mkdir(fullfile(TopLevelFolder, 'Automated FD Calculation Results - x4'), Entries{n});
      end
    case 'Imresize - 0.25 mm pixels - cubic'
      if (exist(fullfile(TopLevelFolder, 'Automated FD Calculation Results - 0.25 mm pixels', Entries{n}), 'dir') ~= 7)
        mkdir(fullfile(TopLevelFolder, 'Automated FD Calculation Results - 0.25 mm pixels'), Entries{n});
      end
  end      
  
  SourceFolder = fullfile(TopLevelFolder, Entries{n});
  
  switch InterpolationType    
    case 'Imresize - (x4 x4) - cubic'
      TargetFolder = fullfile(TopLevelFolder, 'Automated FD Calculation Results - x4', Entries{n});
    case 'Imresize - 0.25 mm pixels - cubic'
      TargetFolder = fullfile(TopLevelFolder, 'Automated FD Calculation Results - 0.25 mm pixels', Entries{n});
  end
  
  [ ImageStack, SegmentationStack, BinaryMask, PerimeterStack, Conditions, OriginalResolution ] = ...
  pft_ExtractMatchedAndShiftedImages(SourceFolder, AcquisitionOrder, MinimumPixelCount, ConnectedPercentage);

  if isempty(ImageStack)
    rmdir(TargetFolder, 's');
    Data = sprintf('%s, %s %s', Entries{n}, repmat('  ,', [1, 35]), '  ');
    fid = fopen(SummaryFile, 'at');
    fprintf(fid, '%s\n', Data);
    fclose(fid);
    waitbar(n/NDIRS, h1, sprintf('Processed %1d of %1d folders', n, NDIRS));
    continue;
  end

  [ NR, NC, NP ] = size(ImageStack);
  
  switch InterpolationType    
    case 'Imresize - (x4 x4) - cubic'
      OutputResolution = OriginalResolution/4.0;
    case 'Imresize - 0.25 mm pixels - cubic'
      OutputResolution = 0.25;
  end
  
  h2 = waitbar(0, sprintf('Processed 0 of %1d slices', NP), 'Units', 'normalized', 'Position', [0.525 0.45 0.2 0.1]);
  
  set(h2, 'Name', 'Progress - slices');
  
  FD = NaN(1, 20);
  FractalDimensions = repmat({ 'NaN' }, [1, 20]);
  
  for p = 1:NP
      
    switch Conditions{p}
        
      case 'Meagre blood pool'
        Wzor = ImageStack(:, :, p);        
        Segmentation = SegmentationStack(:, :, p); 
        
        pft_WriteOriginallySizedImages(Wzor, Segmentation, p, TargetFolder);
          
        FD(p) = MeagreBloodPool;
        FractalDimensions{p} = 'Meagre blood pool';  
        pft_WriteAllBlankImages(p, TargetFolder, 'Meagre blood pool');
        
      case 'Sparse myocardium' 
        Wzor = ImageStack(:, :, p);        
        Segmentation = SegmentationStack(:, :, p); 
        
        pft_WriteOriginallySizedImages(Wzor, Segmentation, p, TargetFolder);
          
        FD(p) = SparseMyocardium;
        FractalDimensions{p} = 'Sparse myocardium';  
        pft_WriteAllBlankImages(p, TargetFolder, 'Sparse myocardium');
        
      case 'No ROI created' 
        Wzor = ImageStack(:, :, p);        
        Segmentation = SegmentationStack(:, :, p); 
        
        pft_WriteOriginallySizedImages(Wzor, Segmentation, p, TargetFolder);
          
        FD(p) = NoROICreated;
        FractalDimensions{p} = 'No ROI created';
        pft_WriteAllBlankImages(p, TargetFolder, 'No ROI created');        
        
      case 'OK'  
        Wzor = ImageStack(:, :, p);        
        Mask = BinaryMask(:, :, p);
        Segmentation = SegmentationStack(:, :, p);        
        Perimeter = PerimeterStack(:, :, p);
        
        pft_WriteOriginallySizedImages(Wzor, Segmentation, p, TargetFolder);
        pft_WriteOriginallySizedMask(Mask, p, TargetFolder);
    
        s = regionprops(Mask, 'BoundingBox');
        
        Wzor = imcrop(Wzor, s(1).BoundingBox);
        Mask = imcrop(Mask, s(1).BoundingBox);
        Segmentation = imcrop(Segmentation, s(1).BoundingBox);
        Perimeter = imcrop(Perimeter, s(1).BoundingBox);
        
        [ Wzor, Mask, Segmentation, Perimeter ] = ...
        pft_InterpolateImages(Wzor, Mask, Segmentation, Perimeter, OriginalResolution, InterpolationType);  
    
        pft_WriteInputImages(Wzor, Mask, Segmentation, Perimeter, p, TargetFolder);
    
        try
          FD(p) = pft_JC_FractalDimensionCalculation(Wzor, Mask, p, TargetFolder);
          FractalDimensions{p} = sprintf('%.9f', FD(p));
        catch
          FD(p) = FDMeasureFailed;
          Conditions{p} = 'FD measure failed';
          FractalDimensions{p} = 'FD measure failed';
          pft_WriteBlankOutputImages(p, TargetFolder, 'FD measure failed');
        end       
 
    end
       
    waitbar(p/NP, h2, sprintf('Processed %1d of %1d slices', p, NP));
  
  end
  
  waitbar(1, h2, sprintf('Processed %1d of %1d slices', NP, NP));
  
  pause(0.1);
  
  close(h2);
  
  waitbar(n/NDIRS, h1, sprintf('Processed %1d of %1d folders', n, NDIRS));  
  
  % Extract and process the FD values for the current stack, trimmed to the number of slices present
  StackFD = FD(1:NP);
    
  switch DiscardEndSlices
    case 'No'
      S = pft_JC_FDStatistics(StackFD, false);
    case 'Yes'
      S = pft_JC_FDStatistics(StackFD, true);
  end  
  
  % Write out the formatted o/p as text 
  FormattedFDOutput = '';
  
  for c = 1:19
    switch Conditions{c}
      case { 'Meagre blood pool', 'Sparse myocardium', 'No ROI created', 'FD measure failed' }
        FormattedFDOutput = [ FormattedFDOutput sprintf('%s,', FractalDimensions{c}) ];
      case 'OK'
        if isnan(FD(c))
          FormattedFDOutput = [ FormattedFDOutput 'NaN,' ];
        else
          FormattedFDOutput = [ FormattedFDOutput sprintf('%s,', FractalDimensions{c}) ];
        end
    end
  end  
  
  switch Conditions{20}
    case { 'Meagre blood pool', 'Sparse myocardium', 'No ROI created', 'FD measure failed' }
      FormattedFDOutput = [ FormattedFDOutput sprintf('%s', FractalDimensions{20}) ];
    case 'OK'
      if isnan(FD(20))
        FormattedFDOutput = [ FormattedFDOutput 'NaN' ];
      else
        FormattedFDOutput = [ FormattedFDOutput sprintf('%s', FractalDimensions{20}) ];
      end
  end
        
  Data = sprintf('%s,%s,%s,%s,%1d,%.2f,%.9f,%.9f,%1d,%s,%s,%1d,%1d,%.9f,%.9f,%.9f,%.9f,%.9f', ...
                 Entries{n}, ...
                 AcquisitionOrder, InterpolationType, PerimeterType, ...                 
                 MinimumPixelCount, ConnectedPercentage, ...             
                 OriginalResolution, OutputResolution, ...
                 NP, ...
                 FormattedFDOutput, ...
                 DiscardEndSlices, S.SlicesEvaluated, S.SlicesUsed, ...
                 S.MeanGlobalFD, ...
                 S.MeanBasalFD, S.MeanApicalFD, ...
                 S.MaxBasalFD, S.MaxApicalFD);
             
  fid = fopen(SummaryFile, 'at');
  fprintf(fid, '%s\n', Data);
  fclose(fid);
            
end

waitbar(1, h1, sprintf('Processed %1d of %1d folders', NDIRS, NDIRS));

pause(0.1);
  
close(h1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Home, James !
msgbox('Done !', 'Quit');





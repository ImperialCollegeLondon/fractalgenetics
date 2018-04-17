function [ ImageStack, SegmentationStack, BinaryMask, PerimeterStack, Conditions, OriginalResolution ] = ...
         pft_ExtractMatchedAndShiftedImages(Folder, AcquisitionOrder, MinimumPixelCount, ConnectedPercentage)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:   Folder              - Self-explanatory.                                                                                     %
%           AcquisitionOrder    - Base-to-Apex or Apex-to-Base.                                                                         %
%           MinimumPixelCount   - A check for a meagre blood pool.                                                                      %
%           ConnectedPercentage - A check for a myocardium detached from the blood pool (i.e., a count of "B" pixels touching "M".      %
%                                                                                                                                       %
% Outputs:  ImageStack          - Possibly trimmed in the 3rd. dimension.                                                               %
%           SegmentationStack   - Possibly trimmed in the 3rd. dimension.                                                               %
%           BinaryMask          - A stack of areas covering the blood pool and extending beyond it into the myocardium.                 %
%           PerimeterStack      - A stack of boundaries to act as a starting point for the FD calculation.                              %
%           Conditions          - For each slice, one of:                                                                               %
%                                                           'OK',                                                                       %
%                                                           'Meagre blood pool',                                                        %
%                                                           'Sparse myocardium',                                                        %
%                                                           'No ROI created'.                                                           %
%           OriginalResolution  - In mm.                                                                                                %
%                                                                                                                                       %
%                                                                                                                                       %
% PFT - 23. 03. 2018.                                                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create some default return values, in case of an early exit
ImageStack = [];
SegmentationStack = [];
BinaryMask = [];
PerimeterStack = [];
Conditions = {};
OriginalResolution = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the end-diastole image stack - o/p variables will be empty on failure
if (exist(fullfile(Folder, 'sa_ED.nii.gz'), 'file') ~= 2)
  return;
end

X = load_nii(fullfile(Folder, 'sa_ED.nii.gz'));
A = X.img;

% Convert A to int16 format if necessary - we know that this works
if ~isa(A, 'int16')
  A = int16(A);
end

AX0 = X.hdr.hist.qoffset_x;
AY0 = X.hdr.hist.qoffset_y;
AZ0 = X.hdr.hist.qoffset_z;

VoxelSizeA = X.hdr.dime.pixdim(2:4);

[ AR, AC, M ] = size(A);

AR = int32(AR);
AC = int32(AC);
M  = int32(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the end-diastole segmentation stack - o/p variables will be empty on failure
if (exist(fullfile(Folder, 'seg_sa_ED.nii.gz'), 'file') ~= 2)
  if (exist(fullfile(Folder, 'seg_sa_ED.gipl'), 'file') ~= 2)  
    return;
  else
    if ispc
      a = fullfile('Windows', 'convert');
      b = fullfile(Folder, 'seg_sa_ED.gipl');
      c = fullfile(Folder, 'seg_sa_ED.nii.gz');
      
      Cmd = sprintf('%s "%s" "%s" 1>nul 2>nul', a, b, c);
    elseif isunix
      a = fullfile('Linux', 'convert');
      b = fullfile(Folder, 'seg_sa_ED.gipl');
      c = fullfile(Folder, 'seg_sa_ED.nii.gz');
            
      Cmd = sprintf('%s "%s" "%s" > /dev/null', a, b, c);
    end
             
    system(Cmd); 
  end
end

FileExists = false;

while (FileExists == false)
  pause(0.5);
  if (exist(fullfile(Folder, 'seg_sa_ED.nii.gz'), 'file') == 2)
    FileExists = true;
  end;  
end

Y = load_nii(fullfile(Folder, 'seg_sa_ED.nii.gz'));
B = Y.img;

% Convert B to int16 format if necessary - we know that this works
if ~isa(B, 'int16')
  B = int16(B);
end

BX0 = Y.hdr.hist.qoffset_x;
BY0 = Y.hdr.hist.qoffset_y;
BZ0 = Y.hdr.hist.qoffset_z;

VoxelSizeB = Y.hdr.dime.pixdim(2:4);

[ BR, BC, N ] = size(B);

BR = int32(BR);
BC = int32(BC);
N  = int32(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The sign convention is not the obvious one here and disagrees with usage elsewhere (doubtless due to the image orientation)
Shift = - round((AZ0 - BZ0)/VoxelSizeA(3));
Shift = int32(Shift);

[ ImageStack, SegmentationStack ] = pft_RealignImages(A, B, Shift);

% If necessary, flip the image and segmentation in the slice direction - put the base first and the apex last
if strcmpi(AcquisitionOrder, 'Apex to Base')
  ImageStack = flip(ImageStack, 3);
  SegmentationStack = flip(SegmentationStack, 3);
end

% Extract the image dimensions again
[ AR, AC, M ] = size(ImageStack);
[ BR, BC, N ] = size(SegmentationStack);

AR = int32(AR);
AC = int32(AC);
M  = int32(M);

BR = int32(BR);
BC = int32(BC);
N  = int32(N);

% Report the original image resolution
OriginalResolution = VoxelSizeA(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the Conditions array with a default value to be overwritten
Conditions = repmat({ 'OK' }, [1, 20]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look for planes in the segmentation stack CONTAINING a blood pool, as well as myocardium:                                             %
% 0 = Background, 1 = LV blood, 2 = myocardium, 4 = remaining heart - this is the original convention                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Background = 0;
Blood      = 1;
Myocardium = 2;
Heart      = 4;

% The Other category has a use here - to pick up "mis-labelled" heart regions - and another elsewhere (during interpolation)
Other = 3; 

% Correct any possible mis-labelling of the heart - e.g., as "3" rather than "4"
SegmentationStack(SegmentationStack == Other) = Heart;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now, check for a meagre blood pool 
for n = 1:N
  Temp = SegmentationStack(:, :, n);
  Pool = (Temp == Blood);
  
  C = sum(Pool(:));
  
  if (C < MinimumPixelCount)
    Conditions{n} = 'Meagre blood pool';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where the blood pool is adequate in size, attempt to create a standard blood pool perimeter entirely within the myocardium
BinaryMask     = false([BR, BC, N]);
PerimeterStack = false([BR, BC, N]);

PerimeterFound = false([1, N]);

MaxRad = 8;

for n = 1:N

  if strcmpi(Conditions{n}, 'OK')                   
      
    Temp = SegmentationStack(:, :, n);
    Pool = (Temp == Blood);
    Wall = (Temp == Myocardium);
      
    for Radius = MaxRad:-1:1
        
      Area = imdilate(Pool, strel('disk', Radius, 4));
      Edge = bwperim(Area, 8);
      OutOfBounds = Edge & ~Wall;
    
      if isempty(find(OutOfBounds == true, 1, 'first'))
        PerimeterFound(n) = true;
        BinaryMask(:, :, n) = Area;
        PerimeterStack(:, :, n) = Edge;
        break;
      end
      
    end
  
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, where the blood pool was adequate but no default perimeter could be constructed,                                             %
% attempt to generate the improvised alternative                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:N
    
  if strcmpi(Conditions{n}, 'OK') && (PerimeterFound(n) == false)
      
    Temp = SegmentationStack(:, :, n);
    Pool = (Temp == Blood);
    Wall = (Temp == Myocardium);        
      
    [ Total, Adjacent ] = pft_CountBloodPoolBorderPixels(Pool, Wall);
    
    T = double(Total);
    A = double(Adjacent);
    
    if ((A/T)*100.0 < ConnectedPercentage)
      Conditions{n} = 'Sparse myocardium';      
      PerimeterFound(n) = false;                % Not needed, but included for emphasis
      continue;
    end
    
    if (bweuler(Pool, 8) ~= 1)
      Pool = bwconvhull(Pool, 'union');
    end
    
    Hull = bwconvhull(Wall, 'union');
    
    Area = Hull & Pool;
    Area = imdilate(Area, strel('disk', 1, 4));
    
    Edge = bwperim(Area, 8);
    
    if ~isempty(find(Edge == true, 1, 'first'))
      PerimeterFound(n) = true;
      Conditions{n} = 'OK';                     % Not needed, but included for emphasis
      BinaryMask(:, :, n) = Area;
      PerimeterStack(:, :, n) = Edge;
    else
      PerimeterFound(n) = false;
      Conditions{n} = 'No ROI created';         % Not needed, but included for emphasis
    end
    
  end
  
end

end








function FD = pft_JC_FractalDimensionCalculation(CroppedImage, BinaryMask, Slice, OutputFolder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This commented block is Jiashen Cai's code: it appears to give identical results to the next section in a few test cases -    %
% but let's keep it here, ready to be reinstated if needed.                                                                     %                                                                                       
% Level set threshold:                                                                                                          %
% PreThresholdImage = regionfill(CroppedImage, ~BinaryMask);                                                                    %
% Or possibly:                                                                                                                  %
% PreThresholdImage = roifill(CroppedImage, ~BinaryMask);                                                                       %
% PFT - 21/03/2018.                                                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adjust the cropped image for better contrast (JC's technique) and exclude any part of the image outside the selected ROI
AdjustedCroppedImage = imadjust(CroppedImage);

PreThresholdImage = zeros(size(AdjustedCroppedImage), class(AdjustedCroppedImage));
PreThresholdImage(BinaryMask) = AdjustedCroppedImage(BinaryMask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The BFE/LSE step follows.                                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The choice of sign below is JC's; it actually makes no difference because the Sobel operator is symmetric across a boundary
ThresholdImage = ~pft_JC_LevelSetDetection(PreThresholdImage);
ThresholdImage = ThresholdImage.*BinaryMask;

% The next section is adapted from JC's code, with some extra care taken in placing the best "equivalent moments" ellipse
LargestArea = bwconvhull(bwpropfilt(imbinarize(ThresholdImage, 0.5), 'Area', 1, 'largest'));

Stats = regionprops(LargestArea, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');

x0 = Stats.Centroid(1);
y0 = Stats.Centroid(2);
aa = 1.05*Stats.MajorAxisLength/2.0;    % Expand the axis by 5 P.C.
bb = 1.05*Stats.MinorAxisLength/2.0;    % Expand the axis by 5 P.C.
Theta = (pi/180.0)*Stats.Orientation;
cc = cos(Theta);
ss = sin(Theta);

[ rows, cols ] = size(LargestArea);

[ xx, yy ] = meshgrid(1:cols, 1:rows);

XX = cc*(xx - x0) - ss*(yy - y0);
YY = ss*(xx - x0) + cc*(yy - y0);

EquivalentEllipse = (XX.^2/aa^2 + YY.^2/bb^2 <= 1.0);

Delta       = round(0.025*(aa + bb));                           
LargestArea = imdilate(LargestArea, strel('disk', Delta, 4));   % Expand the largest area by roughly the same amount as the ellipse

Surround = LargestArea | EquivalentEllipse;

ThresholdImage = ThresholdImage.*Surround;

% Edge detection - the two-step operation below is JC's method also
[ ~, Threshold ] = edge(ThresholdImage, 'sobel');
EdgeImage = edge(ThresholdImage, 'sobel', 0.5*Threshold);

% Box-counting with FD computation
FD = pft_JC_bxct(EdgeImage, Slice, OutputFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These o/p images will be created and written upon successful completion of the FD calculation.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now output some more images
a = double(CroppedImage);
a = uint8(255.0*normalize01(a));

e = double(EdgeImage);
e = uint8(255.0*normalize01(abs(e)));

t = zeros(size(a), class(a));
t(Surround) = a(Surround);

% Also the level-set computed edge image (JC's calculation)
FileName = sprintf('Edge-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(e, PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end  

% This is a poor choice of name by JC, but let's keep things consistent with his nomenclature
FileName = sprintf('Thresholded-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(imadd(e, t), PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

end


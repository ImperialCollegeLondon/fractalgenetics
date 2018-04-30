function pft_WriteOriginallySizedImages(Wzor, Segmentation, Slice, OutputFolder)

% Output the available i/p images straightaway
a = double(Wzor);
b = double(Segmentation);

a = uint8(255.0*normalize01(a));
b = uint8(60.0*b);

% Output the original uncropped image
FileName = sprintf('Original-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(a, PathName);
  pause(0.05);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

% Output the original uncropped segmentation
FileName = sprintf('Original-Segmentation-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(b, PathName);
  pause(0.05);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

end


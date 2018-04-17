function pft_WriteInputImages(CroppedImage, BinaryMask, Segmentation, Perimeter, Slice, OutputFolder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out some i/p images to the BFE/LSE step before those have a chance to fail, so that we only lose the outputs,           %
% which, of course, are never created.                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output the available i/p images straightaway
a = double(CroppedImage);
b = double(BinaryMask);
c = double(Segmentation);
d = double(Perimeter);

a = uint8(255.0*normalize01(a));
b = uint8(255.0*normalize01(b));
c = uint8(60.0*c);
d = uint8(255.0*normalize01(abs(d)));

% Output the cropped image
FileName = sprintf('Cropped-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(a, PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

% Output the binary mask
FileName = sprintf('Binary-Mask-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(b, PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

% Show where the outer perimeter was drawn (within the image plane)
FileName = sprintf('Perimeter-On-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(imadd(a, d), PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

% Show where the outer perimeter was drawn (within the segmentation plane)
FileName = sprintf('Perimeter-On-Segmentation-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(imadd(c, d), PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

% Show the segmentation on its own - useful for writing the technical paper
FileName = sprintf('Segmentation-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(c, PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

% And finally the perimeter
FileName = sprintf('Perimeter-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(d, PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

end


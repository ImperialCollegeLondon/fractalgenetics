function pft_WriteOriginallySizedMask(Mask, Slice, OutputFolder)

% Output a binary mask image before FD processing
a = double(Mask);
a = uint8(255.0*normalize01(a));

% Output the binary mask at the original uncropped resolution
FileName = sprintf('Original-Mask-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(a, PathName);
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

end


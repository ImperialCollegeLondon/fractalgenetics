function pft_WriteAllBlankImages(Slice, OutputFolder, Caption)

Black = zeros([256, 256], 'uint8');

iptsetpref('ImshowBorder', 'tight');

f = figure('Name', 'Blank o/p', 'MenuBar', 'none', 'NumberTitle', 'off', 'Visible', 'off');
axis off, box off, grid off, hold on

% Output the cropped image - jog the frame into position for the first (invisible) screen capture
imshow(Black, [0, 255]);
pause(0.1);
imshow(Black, [0, 255]);
text(8, 16, sprintf('Slice %1d', Slice), 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 32, 'No cropped image', 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 48, Caption, 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');

G = getframe(f);
X = G.cdata;

FileName = sprintf('Cropped-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(X, PathName, 'png');
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end

clf(f);

pause(0.1);

% Show where the outer perimeter was drawn (within the image plane)
imshow(Black, [0, 255]);
text(8, 16, sprintf('Slice %1d', Slice), 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 32, 'No perimeter on image', 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 48, Caption, 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');

G = getframe(f);
X = G.cdata;

FileName = sprintf('Perimeter-On-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(X, PathName, 'png');
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

clf(f);

pause(0.1);

% Show where the outer perimeter was drawn (within the segmentation plane)
imshow(Black, [0, 255]);
text(8, 16, sprintf('Slice %1d', Slice), 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 32, 'No perimeter on segmentation', 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 48, Caption, 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');

G = getframe(f);
X = G.cdata;

FileName = sprintf('Perimeter-On-Segmentation-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(X, PathName, 'png');
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

clf(f);

pause(0.1);

% Show the segmentation on its own - useful for writing the technical paper
imshow(Black, [0, 255]);
text(8, 16, sprintf('Slice %1d', Slice), 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 32, 'No segmentation', 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 48, Caption, 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');

G = getframe(f);
X = G.cdata;

FileName = sprintf('Segmentation-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(X, PathName, 'png');
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

clf(f);

pause(0.1);

% And finally the perimeter - that completes the i/p images
imshow(Black, [0, 255]);
text(8, 16, sprintf('Slice %1d', Slice), 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 32, 'No perimeter', 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 48, Caption, 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');

G = getframe(f);
X = G.cdata;

FileName = sprintf('Perimeter-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(X, PathName, 'png');
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

clf(f);

pause(0.1);

% Now for the the level-set computed edge image (JC's calculation) - one of two outputs
imshow(Black, [0, 255]);
text(8, 16, sprintf('Slice %1d', Slice), 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 32, 'No edge image', 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 48, Caption, 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');

G = getframe(f);
X = G.cdata;

FileName = sprintf('Edge-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(X, PathName, 'png');
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end  

clf(f);

pause(0.1);

% This is a poor choice of name by JC, but let's keep things consistent with his nomenclature
imshow(Black, [0, 255]);
text(8, 16, sprintf('Slice %1d', Slice), 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 32, 'No thresholded image', 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
text(8, 48, Caption, 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');

G = getframe(f);
X = G.cdata;

FileName = sprintf('Thresholded-Image-Slice-%1d-ED.png', Slice);
PathName = fullfile(OutputFolder, FileName);
FileWritten = false;
while (FileWritten == false)
  imwrite(X, PathName, 'png');
  pause(0.1);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

clf(f);

pause(0.1);

% Don't forget to delete the figure
delete(f);

end




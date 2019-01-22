function pft_WriteBlankOutputImages(Slice, OutputFolder, Caption)

Black = zeros([256, 256], 'uint8');

iptsetpref('ImshowBorder', 'tight');

f = figure('Name', 'Blank o/p', 'MenuBar', 'none', 'NumberTitle', 'off', 'Visible', 'off');
axis off, box off, grid off, hold on

% Now for the the level-set computed edge image (JC's calculation) - one of two outputs - 
% jog the frame into position for the first (invisible) screen capture
imshow(Black, [0, 255]);
pause(0.05);
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
  pause(0.05);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end  

clf(f);

pause(0.05);

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
  pause(0.05);
  if (exist(PathName, 'file') == 2)
    FileWritten = true;
  end
end 

clf(f);

pause(0.05);

% Don't forget to delete the figure
delete(f);

end




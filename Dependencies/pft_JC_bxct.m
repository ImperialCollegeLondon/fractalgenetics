function FD = pft_JC_bxct(EdgeImage, Slice, OutputFolder)

% Padding image to have equal dimensions
s = regionprops(EdgeImage, 'BoundingBox');

if (s(1).BoundingBox(3)*s(1).BoundingBox(4) > 0.25*size(EdgeImage, 1)*size(EdgeImage, 2))
  EdgeImage = imcrop(EdgeImage, s(1).BoundingBox);
end

[height, width] = size(EdgeImage);

if height > width
  padWidth = height-width;
  padHeight = 0;
  dimen = height;
else
  padHeight = width-height;
  padWidth = 0;
  dimen = width;
end

padI = padarray(EdgeImage, [mod(padHeight, 2), mod(padWidth, 2)], 'pre');
padI = padarray(padI, [floor(padHeight/2), floor(padWidth/2)]);

[nBoxA, boxSizeA] = initGrid(padI, dimen);
[nBoxB, boxSizeB] = initGrid(flip(padI, 1), dimen);
[nBoxC, boxSizeC] = initGrid(flip(padI, 2), dimen);
[nBoxD, boxSizeD] = initGrid(flip(flip(padI, 1), 2), dimen);

nBox = min([nBoxA; nBoxB; nBoxC; nBoxD]); % Maximal efficient covering
boxSize = boxSizeA;

totalBoxSizes = numel(boxSize);

p  =   polyfit(log(boxSize), log(nBox), 1);
FD = - p(1);

% Output the raw numbers for re-formatting
FileName = sprintf('Box-Count-Slice-%1d.csv', Slice);
PathName = fullfile(OutputFolder, FileName);

fid = fopen(PathName, 'wt');

fprintf(fid, 'Box size, Box count, Ln(Box size), Ln(Box count), p(1), p(2)\n');

for n = 1:totalBoxSizes
  fprintf(fid, '%d, %d, %.9f, %.9f, %.9f, %.9f\n', boxSize(n), nBox(n), log(boxSize(n)), log(nBox(n)), p(1), p(2));
end

fclose(fid);

% Add the polynomial fit as well
FileName = sprintf('Polynomial-Fit-Slice-%1d.csv', Slice);
PathName = fullfile(OutputFolder, FileName);

fid = fopen(PathName, 'wt');

fprintf(fid, 'Ln(Box size), Ln(Box count)\n');

% PFT - 13-02-2017
% fprintf(fid, '%.9f, %.9f\n', log(boxSize(1)), polyval(p, log(boxSize(1))));
% fprintf(fid, '%.9f, %.9f\n', log(boxSize(npts)), polyval(p, log(boxSize(npts))));

xx = log([2.0, 0.45*dimen]);
yy = polyval(p, xx);

fprintf(fid, '%.9f, %.9f\n', xx(1), yy(1));
fprintf(fid, '%.9f, %.9f\n', xx(2), yy(2));

fclose(fid);

% Box Count Plot
bcFig = figure('Name', 'Box Count Plot', 'MenuBar', 'none', 'NumberTitle', 'off', 'Visible', 'off');
pause(0.05);
set(0, 'CurrentFigure', bcFig);
loglog(boxSizeA, nBoxA, 's-');
xlabel('r, box size (pixels)'); 
ylabel('n(r), box count');

x1 = linspace(2, 0.45*dimen);
y1 = exp(polyval(p, log(x1)));
hold on;
loglog(x1, y1);

legend('Box Count', 'Regression');

FileName = sprintf('Box-Count-Plot-Slice-%1d-ED.png', Slice);

export_fig(gca, fullfile(OutputFolder, FileName), '-png', '-m3');

pause(0.05);

delete(bcFig);

    function [nBox, boxSize] = initGrid(Image, dimen)
        
        startBoxSize = floor(0.45*dimen);
        curBoxSize   = startBoxSize;
        
        nBox    = zeros(1, (startBoxSize-2+1));
        boxSize = (startBoxSize:-1:2);
        
        for sizeCount = 1:(startBoxSize-2+1)
            curBoxSize = boxSize(sizeCount);
            
            for macroY = 1:ceil(dimen/curBoxSize)
                for macroX = 1:ceil(dimen/curBoxSize)
                    boxYinit = (macroY-1)*curBoxSize+1;
                    boxXinit = (macroX-1)*curBoxSize+1;
                    boxYend = min(macroY*curBoxSize,dimen);
                    boxXend = min(macroX*curBoxSize,dimen);
                    
                    boxFound = false;
                    for curY = boxYinit:boxYend
                        for curX = boxXinit:boxXend
                            if Image(curY,curX)
                                boxFound = true;
                                nBox(sizeCount) = nBox(sizeCount) + 1;
                                break;
                            end
                        end
                        
                        if boxFound == true
                            break;
                        end
                    end                    
                end
            end
        end    
    end

end
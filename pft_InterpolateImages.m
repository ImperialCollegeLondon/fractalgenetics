function [ W, M, S, P ] = pft_InterpolateImages(Wzor, Mask, Segmentation, Perimeter, OriginalResolution, InterpolationType)

% Extract the common (original) image size
[ MR, MC ] = size(Wzor);

% Re-assign values in the Segmentation image to make the interpolation and rounding more straightforward
Background = 0;
Blood      = 1;
Myocardium = 2;
Heart      = 4;

% The Other category is used here as a temporary swap value, to avoid the creation of interpolated pixels with a value of 3,
% since these would then be unassigned to a tissue type within the traditional "bit-binary" labelling scheme
Other = 3;

switch InterpolationType
    
  case 'Imresize - (x4 x4) - cubic'     
      
    W = imresize(Wzor, 4, 'cubic');
    M = imresize(Mask, 4, 'cubic');
    
    Segmentation(Segmentation == Heart) = Other;
    S = imresize(Segmentation, 4, 'cubic');
    S(S == Other) = Heart;
     
    P = imresize(Perimeter, 4, 'cubic');
    
  case 'Imresize - 0.25 mm pixels - cubic'
      
    NewPixelSize  = 0.25;
    Magnification = OriginalResolution/NewPixelSize;
    
    NR = uint16(round(Magnification*MR));
    NC = uint16(round(Magnification*MC));
    
    W = imresize(Wzor, [NR, NC], 'cubic');
    M = imresize(Mask, [NR, NC], 'cubic');
    
    Segmentation(Segmentation == Heart) = Other;
    S = imresize(Segmentation, [NR, NC], 'cubic');
    S(S == Other) = Heart;
        
    P = imresize(Perimeter, [NR, NC], 'cubic');
    
end

end
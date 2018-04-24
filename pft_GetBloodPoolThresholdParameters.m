function [ MinimumPixelCount, ConnectedPercentage ] = pft_GetBloodPoolThresholdParameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pft_BloodPoolThresholdParameters                                                                      %
%                                                                                                       %
% Input limits to test for a scanty blood pool, or one which is poorly connected to the myocardium.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Options.Resize = 'off';
Options.WindowStyle = 'modal';
Options.Interpreter = 'tex';

Prompt = { 'Blood pool pixel count threshold: ', 'Minimum connected percentage: ' };

Starts = { '60', '50.0' };  % 60-65 pixels is best for Genscan, 38 for Biobank

Layout = zeros(2, 2, 'int16');
Layout(:, 1) = 1;
Layout(:, 2) = 60;

Answers = inputdlg(Prompt, 'Blood pool thresholding parameters', Layout, Starts, Options);

Amended = false;

if (length(Answers) == length(Starts))
  MinimumPixelCount = int32(str2double(Answers{1}));
  ConnectedPercentage = str2double(Answers{2});  
  
  if ~isnumeric(MinimumPixelCount) 
    MinimumPixelCount = int32(str2double(Starts{1}));
    Amended = true;
  elseif isnan(MinimumPixelCount) || isinf(MinimumPixelCount)
    MinimumPixelCount = int32(str2double(Starts{1}));
    Amended = true;
  end
  
  if ~isnumeric(ConnectedPercentage) 
    ConnectedPercentage = str2double(Starts{2});
    Amended = true;
  elseif isnan(ConnectedPercentage) || isinf(ConnectedPercentage)
    ConnectedPercentage = str2double(Starts{2});
    Amended = true;
  end  
else
  MinimumPixelCount = int32(str2double(Starts{1}));
  ConnectedPercentage = str2double(Starts{2});
  Amended = true;
end

if (MinimumPixelCount < 10)
  MinimumPixelCount = 10;
  Amended = true;
elseif (MinimumPixelCount > 1000)
  MinimumPixelCount = 1000;
  Amended = true;
end

if (ConnectedPercentage < 25.0)
  ConnectedPercentage = 25.0;
  Amended = true;
elseif (ConnectedPercentage > 100.0)
  ConnectedPercentage = 100.0;
  Amended = true;
end

if (Amended == true)
  beep;
  Warning = { 'Input amended:', sprintf('Threshold = %d pixels', MinimumPixelCount), sprintf('Connected percentage = %.2f', ConnectedPercentage) };
  Title   =   'Error correction';
  h = warndlg(Warning, Title, 'modal');                
  uiwait(h);
  delete(h);
end
  
end


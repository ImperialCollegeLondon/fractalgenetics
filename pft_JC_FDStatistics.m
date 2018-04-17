function Statistics = pft_JC_Statistics(FD, Discard)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2017 Cai Jiashen.                                                                           %
% Edited for clarity - PFT - 07/03/2018.                                                                    %
% The i/p array FD will always be large enough that there will not be a need to check whether it's empty.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create some default outputs in case of an early return
Statistics = struct('SlicesEvaluated',  0,   'SlicesUsed',  0, ...
                    'MeanMeanGlobalFD', 0.0, ...
                    'MeanBasalFD',      0.0, 'MaxBasalFD',  0.0, ...
                    'MeanApicalFD',     0.0, 'MaxApicalFD', 0.0);
                
% Suppress zero values for the purposes of averaging, but avoid making an FP comparison with 0.0
THRESHOLD = 0.1;

FD(FD < THRESHOLD) = NaN;

% Count the slices from the beginning, not from the first non-NaN value
Statistics.SlicesEvaluated = numel(FD);

n = Statistics.SlicesEvaluated;
q = floor(n/2);
r = rem(n, 2);

if (Discard == false)
  Statistics.SlicesUsed = n;
  LG = 1;           % Lower global index
  UG = n;           % Upper global index
  LB = 1;           % Lower basal index
  UB = q;           % Upper basal index
  LA = 1 + q + r;   % Lower apical index
  UA = n;           % Upper apical index  
elseif (Discard == true)
  Statistics.SlicesUsed = n - 2;
  LG = 2;           % Lower global index
  UG = n - 1;       % Upper global index
  LB = 2;           % Lower basal index
  UB = q;           % Upper basal index
  LA = 1 + q + r;   % Lower apical index
  UA = n - 1;       % Upper apical index
end

% Trim the data if necessary, but in any case, assign the global, basal and apical arrays
GlobalFD = FD(LG:UG);
BasalFD  = FD(LB:UB);
ApicalFD = FD(LA:UA);

% Calculate the statistics
Statistics.MeanGlobalFD = nanmean(GlobalFD);

Statistics.MeanApicalFD = nanmean(ApicalFD);
Statistics.MaxApicalFD  = nanmax(ApicalFD);  

Statistics.MeanBasalFD  = nanmean(BasalFD);
Statistics.MaxBasalFD   = nanmax(BasalFD);

end

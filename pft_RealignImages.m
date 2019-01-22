function [ AA, BB ] = pft_RealignImages(A, B, S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:   A, B - 3-D arrays, with dimensions 1, 2 assumed the same but 3 different;                                       %
%           i.e., A has size [ R, C, M ] and B has size [ R, C, N ].                                                        %
%                                                                                                                           %
%           S - an integer shift (of B to the right w.r.t. A).                                                              %
%                                                                                                                           %
% Outputs:  AA, BB - the overlap of images A, B in the third direction.                                                     %
%                                                                                                                           %
% Null results are returned if the arrays are shifted past each other with no overlap.                                      %
%                                                                                                                           %
% PFT - 15. 11. 2016.                                                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ AR, AC, M ] = size(A);    % Deliberate notation - we expect dimensions 1, 2 to be the same but 3 to be different
[ BR, BC, N ] = size(B);    % Deliberate notation - we expect dimensions 1, 2 to be the same but 3 to be different

% YOU ARE HERE
AR = int32(AR);
AC = int32(AC);
M  = int32(M);

% YOU ARE HERE
BR = int32(BR);
BC = int32(BC);
N  = int32(N);

Returning = false;

if (AR ~= BR) 
  fprintf('Incompatible 1st. dimensions between two input arrays.\n');  % Shouldn't happen if checked elsewhere
  Returning = true;
end

if (AC ~= BC) 
  fprintf('Incompatible 2nd. dimensions between two input arrays.\n');  % Shouldn't happen if checked elsewhere
  Returning = true;
end

if (Returning == true)
  AA = [];
  BB = [];
  return;
end

if (S == 0)
  Lo = 1;
  Hi = min(M, N);
  AA = A(:, :, Lo:Hi);
  BB = B(:, :, Lo:Hi);
elseif (S > 0)
  LoA = 1 + S;
  HiA = min(M, N + S);
  LoB = 1;
  HiB = min(M - S, N);
  if (LoA > M)
    fprintf('Array B right-shifted past array A - no overlap.\n');
    AA = [];
    BB = [];
  else
    AA = A(:, :, LoA:HiA);
    BB = B(:, :, LoB:HiB);    
  end
else
  LoB = 1 - S;
  HiB = min(M - S, N);
  LoA = 1;
  HiA = min(M, N + S);
  if (LoB > N)
    fprintf('Array A right-shifted past array B - no overlap.\n');
    AA = [];
    BB = [];
  else
    AA = A(:, :, LoA:HiA);
    BB = B(:, :, LoB:HiB);
  end
end

end


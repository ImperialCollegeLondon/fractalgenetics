function [ Total, Adjacent ] = pft_CountBloodPoolBorderPixels(B, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A function to count the total number of border pixels in a BW blood-pool array, and those touching myocardium, with a connectivity of 4.  %
%                                                                                                                                           %
% PFT - 19. 03. 2018.                                                                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the perimeter of the blood pool
P = bwperim(B, 8);

% Isolate the bright pixels
[ FR, FC ] = find(P == true);

% Count how many there are - and by the time this function is called, there will certainly be some
Total = numel(FR);

% Set a starting value for the number of adjacent pixels, but quit if the blood pool is too small
Adjacent = 0;

% Now count those pixels adjacent to the myocardium, with a connectivity of 8
[ NR, NC ] = size(P);

for n = 1:Total
  Row = FR(n);
  Col = FC(n);
  
  lr = Row - 1;
  ur = Row + 1;
  lc = Col - 1;
  uc = Col + 1;

  if (1 <= lr) && (lr <= NR) && (1 <= ur) && (ur <= NR) && (1 <= lc) && (lc <= NC) && (1 <= uc) && (uc <= NC)
    if (M(Row, lc) || M(Row, uc) || M(lr, Col) || M(ur, Col)) || M(lr, lc) || M(lr, uc) || M(ur, lc) || M(ur, uc)
      Adjacent = Adjacent + 1;
    end
  end
end

end


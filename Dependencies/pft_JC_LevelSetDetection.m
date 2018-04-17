function contourI = pft_JC_LevelSetDetection(Image)

sigma   = 4;
epsilon = 3;

Img = double(Image);
A   = 255;
Img = A*normalize01(Img);
nu  = 0.001*A^2;            % Coefficient of arc length term

iter_outer = 100;           % Outer iteration for level set evolution
iter_inner = 10;            % Inner iteration for level set evolution

timestep = 0.1;
mu = 1;                     % Coefficient for distance regularization term (regularize the level set function)

c0 = 1;

% Initialize level set function
initialLSF = c0*ones(size(Img));
[xSize, ySize] = size(Img);
initialLSF(round(0.20*xSize):round(0.80*ySize),round(0.20*xSize):round(0.80*ySize)) = - c0;
u = initialLSF;

b = ones(size(Img));                                        % Initialize the bias field

K    = fspecial('gaussian', round(2*sigma)*2 + 1, sigma);   % Gaussian kernel
KI   = conv2(Img, K, 'same');                               % Apparently unused
KONE = conv2(ones(size(Img)), K, 'same');

[row, col] = size(Img);
N = row*col;
  
for n=1:iter_outer
  [u, b, C] = lse_bfe(u, Img, b, K, KONE, nu, timestep, mu, epsilon, iter_inner);    
end

contourI = (u > 0);

end


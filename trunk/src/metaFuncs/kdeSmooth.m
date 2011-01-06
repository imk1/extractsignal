function x = kdeSmooth(x,halfWidth,kernelType,normFlag)
% Smoothes x using kernel density smoothing
% function [x] = kdesmooth(x,halfWidth,kernelType,normFlag)
% x: input vector
% halfWidth: kernel halfWidth
% kernelType: name of smoothing kernel
% normFlag: if set to true kernel smoothing is equivalent to smooth averaging. If not set, it is
%           equivalent to a smooth sum.

sizeX = size(x);
assert((min(sizeX)==1),'input vector x cannot be a matrix');

isRowVector = (sizeX(1)==1);

% Convert x to column vector
if isRowVector
    x = x(:);
end

% Generate kernel
kval = generateKernel(kernelType,halfWidth,normFlag);
normKval = generateKernel(kernelType,halfWidth,1);

% Compute convolution
denom = isfinite(x); % Check for NaNs or Infs in x
hasNaN = any(~denom);
if hasNaN
    x(~denom) = 0; % Replace NaN and Infs with 0
end
% The division normalizes for missing values represented by NaNs or Infs
x = conv2(x,kval,'same')./conv2(double(denom),normKval,'same');

% Reset original NaN values to NaN
x(~denom) = NaN;
x = single(x);

if isRowVector
    x = x';
end

end

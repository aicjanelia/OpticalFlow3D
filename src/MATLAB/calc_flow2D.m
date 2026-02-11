function [vx,vy,rel] = calc_flow2D(images ,xySig, tSig, wSig)
% The calc_flow2D function calculates optical flow velocities for a single
% z-slice over time. Surrounding images in time are necessary to perform
% the calculations. To peform calculations on an entire timelapse, see
% parse_flow.m
%
% This script uses the convention that (0,0) is located in the upper-left
% corner of an image. This is inline with conventions used in other
% programs (e.g., ImageJ/FIJI), but note that it means that positive
% y-velocities point down, which can be non-intuitive. In MATLAB, rows are
% the first dimension, so the code adopts the convention that dimension 1 =
% y and dimension 2 = x.
%
% USAGE: [vx,vy,rel] = calc_flow2D(images,xySig, tSig, wSig)
%
% INPUTS:
% images = 3D matrix with dimensions N_Y, N_X, N_T
%           N_T should be odd as only the central timepoint will be analyzed.
%           N_T must be greater than or equal to 2*3*tSig+1.
% xySig  = sigma value for smoothing in all spatial dimensions. Default 3.
%           Larger values remove noise but remove spatial detail.
% tSig   = sigma value for smoothing in the temporal dimension. Default 1.
%           Larger values remove noise but remove temporal detail.
% wSig   = sigma value for Lucas-Kanade neighborhood. Default is 4.
%           Larger values include a larger neighboorhood in the
%           Lucas-Kanade constraint and will smooth over small features.
%
% OUTPUTS:
% vx     = Velocity in the x direction, reported as pixels/frame
% vy     = Velocity in the y direction, reported as pixels/frame
% rel    = Reliability, the smallest eigenvalue of (A'wA)
%           This is a measure of confidence in the linear algebra solution
%           and can be used to mask the velocities for downstream analysis.
%
% Change Log:
% 2025/02/03 Rachel M. Lee - Created function
% 2025/03/20 Rachel M. Lee - Updated for MATLAB 2024b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the images are 2D + time
if length(size(images)) ~= 3
    error('ERROR: Input image must be a 3D matrix with dimensions N_Y, N_X, N_T')
end

% Implement default values if they are not specified
if ~exist('xySig','var') || isempty(xySig)
    disp('Using default xyzSig value of 3...')
    xySig = 3;
end
if ~exist('tSig','var') || isempty(tSig)
    tSig = 1;
    disp('Using default tSig value of 1...')
end
if ~exist('wSig','var') || isempty(wSig)
    wSig = 4;
    disp('Using default wSig value of 4...')
end

% Check image size against tSig
Nt =  size(images,3);
if Nt < 6*tSig+1
    % The kernel size is 3*tSig in forward and back time (or 6*tSig total)
    % There is also a central pixel, so need at least 6*tSig+1 images in t
    error('ERROR: Input images will lead to edge effects. N_T must be >= 6*tSig+1')
end
if ~mod(Nt,2) % enforce odd number
    error('ERROR: Input images must have an odd number of timepoints. Only the central time point is analyzed')
end
NtSlice = ceil(Nt/2);

images = double(images); % floating point calculations

%% Set up filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Common spatial terms
x = -ceil(3*xySig):ceil(3*xySig);
xySig2 = xySig/4;
y = -ceil(3*xySig2):ceil(3*xySig2);
fderiv = exp(-x.*x/2/xySig/xySig)/sqrt(2*pi)/xySig;
fsmooth = exp(-y.*y/2/xySig2/xySig2)/sqrt(2*pi)/xySig2;
gderiv = x/xySig/xySig;
gsmooth = 1;

%   Build y-gradient filter kernels (along first dimension)
yFil1 = (fderiv.*gderiv)';
xFil1 = (fsmooth.*gsmooth);
%   Build x-gradient filter kernels (along second dimension)
yFil2 = (fsmooth.*gsmooth)';
xFil2 = (fderiv.*gderiv);

%   Build t-gradient filter kernels (t = third dimension)
t = -ceil(3*tSig):ceil(3*tSig);
fx = exp(-x.*x/2/xySig/xySig)/sqrt(2*pi)/xySig;
ft = exp(-t.*t/2/tSig/tSig)/sqrt(2*pi)/tSig;
gx = 1;
gt = t/tSig/tSig;
yFil3 = (fx.*gx)';
xFil3 = yFil3';
tFil3 = permute(ft.*gt,[3 1 2]);

% Structure tensor -- Lucas Kanade neighborhood filter
wRange = -ceil(3*wSig):ceil(3*wSig);
gw = exp(-wRange.*wRange/2/wSig/wSig)/sqrt(2*pi)/wSig;
yFil4 = (gw)';
xFil4 = (gw);

clear gderiv gsmooth gt gw gx ft fx fsmooth fderiv x y t gw wRange

%% Spatial and Temporal Gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial gradients require only the frame of interest, while  the temporal
% gradient requires N_T >= 2*3*tSig+1. Keep only the relevant slice of dtI
% after it is calculated.

% dt (dimension = 3)
% In two steps to save memory.
dtI = imfilter(images, tFil3, 'replicate');
dtI = dtI(:,:,NtSlice);
images = images(:,:,NtSlice); % now we only need the slice of interest, decrease memory usage.
dtI = imfilter(imfilter(dtI, yFil3, 'replicate'), xFil3, 'replicate');
clear xFil3 yFil3 tFil3

% dy (dimension = 1)
dyI = imfilter(imfilter(images, yFil1, 'replicate'), xFil1, 'replicate'); % Filtering to calculate the gradient
clear xFil1 yFil1

% dx (dimension = 2)
dxI = imfilter(imfilter(images, yFil2, 'replicate'), xFil2, 'replicate');
clear xFil2 yFil2

clear images % No longer needed for calculations; clear to save memory

%% Structure Tensor Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following calculations are for the individual elements of the
% matrices required for the optical flow calculation, incorporating
% Gaussian weighting into the Lucas-Kanade constraint.

% Time componenents
wdtx = imfilter(imfilter(dxI.*dtI, yFil4, 'replicate'), xFil4, 'replicate');
wdty = imfilter(imfilter(dyI.*dtI, yFil4, 'replicate'), xFil4, 'replicate');
clear dtI

% Spatial Components
wdxy = imfilter(imfilter(dxI.*dyI, yFil4, 'replicate'), xFil4, 'replicate');
wdx2 = imfilter(imfilter(dxI.*dxI, yFil4, 'replicate'), xFil4, 'replicate');
clear dxI
wdy2 = imfilter(imfilter(dyI.*dyI, yFil4, 'replicate'), xFil4, 'replicate');
clear dyI
clear xFil4 yFil4


%% Optical Flow calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation is v = (A' w A)^-1 A' w b
% A = -[dxI dyI]
% b = [dtI]
% w multiplication is incorporated in the structure tensor inputs above
% A' w b = -[wdtx wdty]  (minus sign because of negative sign on A)
% (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]

% Determinant
determinant = (wdx2.*wdy2) - (wdxy.^2);

% A^-1 = [a b ; c d]^-1 = (1/det(A))[d - b; -c a]
vx = ((determinant + eps).^-1).*((wdy2.*-wdtx)+(-wdxy.*-wdty));
vy = ((determinant + eps).^-1).*((-wdxy.*-wdtx)+(wdx2.*-wdty));

clear wdtx wdty wdxy

%% Eigenvalues for Reliability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve det(A^T w A - lamda I) = 0
% (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]

trace = wdx2 + wdy2;

clear wdx2 wdy2

L1 = (trace + sqrt(trace.^2 - 4*determinant))/2;
L2 = (trace - sqrt(trace.^2 - 4*determinant))/2;
rel = real(min(L1, L2));

clear L1 L2

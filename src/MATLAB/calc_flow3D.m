function [vx,vy,vz,rel] = calc_flow3D(images ,xyzSig, tSig, wSig)
% The calc_flow3D function calculates optical flow velocities for a single
% z-stack of images. Surrounding z-stacks in time are necessary to perform
% the calculations. To peform calculations on an entire timelapse, see
% parse_flow.m
%
% This script uses the convention that (0,0) is located in the upper-left
% corner of an image. This is inline with conventions used in other
% programs (e.g., ImageJ/FIJI), but note that it means that positive
% y-velocities point down, which can be non-intuitive. In MATLAB, rows are
% the first dimension, so the code adopts the convention that dimension 1 =
% y and dimension 2 = x.
% Z = 0 is the first tif image in the multi-tif file, which is dimension =
% 3 in MATLAB. Positive z-velocity is from the first slice towards the last
% slice. This may or may not correspond to positive z in the lab space
% depending on how the images were acquired.
%
% USAGE: [vx,vy,vz,rel] = calc_flow3D(images ,xyzSig, tSig, wSig)
%
% INPUTS:
% images = 4D matrix with dimensions N_Y, N_X, N_Z, N_T
%           N_T should be odd as only the central timepoint will be analyzed.
%           N_T must be greater than or equal to 3*tSig+1.
% xyzSig = sigma value for smoothing in all spatial dimensions. Default 3.
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
% vz     = Velocity in the z direction, reported as pixels/frame
% rel    = Reliability, the smallest eigenvalue of (A'wA)
%           This is a measure of confidence in the linear algebra solution
%           and can be used to mask the velocities for downstream analysis.
%
% Change Log:
% 2025/01/27 Rachel M. Lee - Created function
% 2025/03/20 Rachel M. Lee - Updated for MATLAB 2024b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the images are 3D + time
if length(size(images)) ~= 4
    error('ERROR: Input image must be a 4D matrix with dimensions N_Y, N_X, N_Z, N_T')
end

% Implement default values if they are not specified
if ~exist('xyzSig','var') || isempty(xyzSig)
    disp('Using default xyzSig value of 3...')
    xyzSig = 3;
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
Nt =  size(images,4);
if Nt < 6*tSig+1
    % The kernel size is 3*tSig in forward and back in time (or 6*tSig total)
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
x = -ceil(3*xyzSig):ceil(3*xyzSig);
xySig2 = xyzSig/4;
y = -ceil(3*xySig2):ceil(3*xySig2);
fderiv = exp(-x.*x/2/xyzSig/xyzSig)/sqrt(2*pi)/xyzSig;
fsmooth = exp(-y.*y/2/xySig2/xySig2)/sqrt(2*pi)/xySig2;
gderiv = x/xyzSig/xyzSig;
gsmooth = 1;

%   Build y-gradient filter kernels (along first dimension)
yFil1 = (fderiv.*gderiv)';
xFil1 = (fsmooth.*gsmooth);
zFil1 = permute(xFil1,[3 1 2]);
%   Build x-gradient filter kernels (along second dimension)
yFil2 = (fsmooth.*gsmooth)';
xFil2 = (fderiv.*gderiv);
zFil2 = permute(yFil2,[3 2 1]);
%   Build z-gradient filter kernels (along third dimension)
yFil3 = (fsmooth.*gsmooth)';
xFil3 = (fsmooth.*gsmooth);
zFil3 = permute(fderiv.*gderiv,[3 1 2]);

%   Build t-gradient filter kernels (t = fourth dimension)
t = -ceil(3*tSig):ceil(3*tSig);
fx = exp(-x.*x/2/xyzSig/xyzSig)/sqrt(2*pi)/xyzSig;
ft = exp(-t.*t/2/tSig/tSig)/sqrt(2*pi)/tSig;
gx = 1;
gt = t/tSig/tSig;
yFil4 = (fx.*gx)';
xFil4 = yFil4';
zFil4 = permute(xFil4,[3 1 2]);
tFil4 = permute(ft.*gt, [4, 1, 3, 2]);

% Structure tensor -- Lucas Kanade neighborhood filter
wRange = -ceil(3*wSig):ceil(3*wSig);
gw = exp(-wRange.*wRange/2/wSig/wSig)/sqrt(2*pi)/wSig;
yFil5 = (gw)';
xFil5 = (gw);
zFil5 = permute(gw,[3 1 2]);

clear gderiv gsmooth gt gw gx ft fx fsmooth fderiv x y t gw wRange

%% Spatial and Temporal Gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial gradients require only the frame of interest, while  the temporal
% gradient requires N_T >= 2*3*tSig+1. Keep only the relevant slice of dtI
% after it is calculated.

% dt (dimension = 4)
% In two steps to save memory
dtI = imfilter(images, tFil4, 'replicate');
images = images(:,:,:,NtSlice);
dtI = dtI(:,:,:,NtSlice);
dtI = imfilter(imfilter(imfilter(dtI, yFil4, 'replicate'), xFil4, 'replicate'), zFil4, 'replicate');
clear xFil4 yFil4 zFil4 tFil4

% dy (dimension = 1)
dyI = imfilter(imfilter(imfilter(images, yFil1, 'replicate'), xFil1, 'replicate'), zFil1, 'replicate'); % Filtering to calculate the gradient
clear xFil1 yFil1 zFil1

% dx (dimension = 2)
dxI = imfilter(imfilter(imfilter(images, yFil2, 'replicate'), xFil2, 'replicate'), zFil2, 'replicate');
clear xFil2 yFil2 zFil2

% dz (dimension = 3)
dzI = imfilter(imfilter(imfilter(images, yFil3, 'replicate'), xFil3, 'replicate'), zFil3, 'replicate');
clear xFil3 yFil3 zFil3

clear images % No longer needed for calculations; clear to save memory


%% Structure Tensor Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following calculations are for the individual elements of the
% matrices required for the optical flow calculation, incorporating
% Gaussian weighting into the Lucas-Kanade constraint.

% Time componenents
wdtx = imfilter(imfilter(imfilter(dxI.*dtI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
wdty = imfilter(imfilter(imfilter(dyI.*dtI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
wdtz = imfilter(imfilter(imfilter(dzI.*dtI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
clear dtI

% Spatial Components
wdx2 = imfilter(imfilter(imfilter(dxI.*dxI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
wdxy = imfilter(imfilter(imfilter(dxI.*dyI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
wdxz = imfilter(imfilter(imfilter(dxI.*dzI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
clear dxI
wdyz = imfilter(imfilter(imfilter(dzI.*dyI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
wdy2 = imfilter(imfilter(imfilter(dyI.*dyI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
clear dyI
wdz2 = imfilter(imfilter(imfilter(dzI.*dzI, yFil5, 'replicate'), xFil5, 'replicate'), zFil5,'replicate');
clear dzI
clear xFil5 yFil5 zFil5


%% Optical Flow calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation is v = (A' w A)^-1 A' w b
% A = -[dxI dyI dzI]
% b = [dtI]
% w multiplication is incorporated in the structure tensor inputs above
% A' w b = -[wdtx wdty wdtz]  (minus sign because of negative sign on A)
% (A' w A) = [a=wdx2 b=wdxy c=wdxz ; d=wdxy e=wdy2 f=wdyz ; g=wdxz h=wdyz i=wdz2]
% (A' w A)^-1 = 1/determinant * [A D G ; B E H ; C F I];

% Determinant
determinant = (wdx2.*wdy2.*wdz2) + (2*wdxy.*wdxz.*wdyz) - (wdy2.*wdxz.^2) - (wdz2.*wdxy.^2) - (wdx2.*wdyz.^2);

% using inverse notation from here: https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices
%     A = wdy2.*wdz2 - wdyz.*wdyz; %ei-fh
%     B = wdyz.*wdxz - wdxy.*wdz2; %fg-di
%     C = wdxy.*wdyz - wdy2.*wdxz; %dh-eg
%     D = wdxz.*wdyz - wdxy.*wdz2; %ch-bi
%     E = wdx2.*wdz2 - wdxz.*wdxz; %ai-cg
%     F = wdxy.*wdxz - wdx2.*wdyz; %bg-ah
%     G = wdxy.*wdyz - wdxz.*wdy2; %bf-ce
%     H = wdxz.*wdxy - wdx2.*wdyz; %cd-af
%     I = wdx2.*wdy2 - wdxy.*wdxy; %ae-bd
% 
%     vx = ((determinant + eps).^-1).*(A.*wdtx + D.*wdty + G.*wdtz);
%     vy = ((determinant + eps).^-1).*(B.*wdtx + E.*wdty + H.*wdtz);
%     vz = ((determinant + eps).^-1).*(C.*wdtx + F.*wdty + I.*wdtz);

% Velocity calculations
vx = -((determinant + eps).^-1).*((wdy2.*wdz2 - wdyz.*wdyz).*wdtx + (wdxz.*wdyz - wdxy.*wdz2).*wdty + (wdxy.*wdyz - wdxz.*wdy2).*wdtz);
vy = -((determinant + eps).^-1).*((wdyz.*wdxz - wdxy.*wdz2).*wdtx + (wdx2.*wdz2 - wdxz.*wdxz).*wdty + (wdxz.*wdxy - wdx2.*wdyz).*wdtz);
vz = -((determinant + eps).^-1).*((wdxy.*wdyz - wdy2.*wdxz).*wdtx + (wdxy.*wdxz - wdx2.*wdyz).*wdty + (wdx2.*wdy2 - wdxy.*wdxy).*wdtz);

clear determinant
clear wdtx wdty wdtz

%% Eigenvalues for Reliability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve det(A^T w A - lamda I) = 0
% (A' w A) = [a=wdx2 b=wdxy c=wdxz ; d=wdxy e=wdy2 f=wdyz ; g=wdxz h=wdyz i=wdz2]
% det(A^T w A - lambda I) = (a-lambda)(e-lambda)(i-lambda) + 2*b*c*f - c^2(e-lambda) - b^2(i-lamda) - f^2(a-lambda)

% Assemble the structure tensor
w = NaN*ones(3,3,size(wdx2,1),size(wdx2,2),size(wdx2,3));
w(1,1,:,:,:) = wdx2;
clear wdx2
w(1,2,:,:,:) = wdxy;
w(1,3,:,:,:) = wdxz;
w(2,1,:,:,:) = wdxy;
clear wdxy
w(2,2,:,:,:) = wdy2;
clear wdy2
w(2,3,:,:,:) = wdyz;
w(3,1,:,:,:) = wdxz;
clear wdxz
w(3,2,:,:,:) = wdyz;
clear wdyz
w(3,3,:,:,:) = wdz2;
clear wdz2

% Eigenvalues
rel = pageeig(w); % 3 values
rel = real(squeeze(min(rel))); % Take the minimum
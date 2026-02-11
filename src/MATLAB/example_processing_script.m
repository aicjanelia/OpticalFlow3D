%% Example 1: Process a single experiment saved as a sequence of tif files

% Start with a clear workspace
clc
clear
close all

%%%% Set up information about the experiment files
% Folder containing sequence of tif files.
imDir = 'C:\fullpath\to\image\files\';
% File name of images of interest. Insert * as wildcards for values that
% change across the sequence.
imName = 'scan_CamA_ch0_CAM1_stack*_488nm_*msec_*msecAbs_000x_000y_000z_0000t_decon'; 
% There is one tif per timepoint, so this is a "SequenceT" type movie
fileType = 'SequenceT';
% This is 3D data (a z-stack)
spatialDimensions = 3;

%%%% Set up optical flow parameters
% Spatial smoothing (voxels)
xyzSig = 3;
% Temporal smoothing (frames)
tSig = 1;
% Lucas-Kanade neighborhood size (voxels)
wSig = 4;

%%%% Run the optical flow
% process_flow takes care of parsing file names, etc.
% It will display progress in the Command Window as it progresses
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 2: Process multiple folders of experiments
% Each experiment was saved as a single tif containing all Z & T using ImageJ

% Start with a clear workspace
clc
clear
close all

%%%% Set up information about the experiment files
% List all folders to be processed.
% All .tif files in each folder will be processed.
toProcess = {'C:\fullpath\to\image\files\1\',...
    'C:\fullpath\to\image\files\2\'};
% There is one tif for the entire time lapse, so this is a "OneTif" type movie
fileType = 'OneTif';
% This is 3D data (a z-stack)
spatialDimensions = 3;

%%%% Set up optical flow parameters
% Spatial smoothing (voxels)
xyzSig = 3;
% Temporal smoothing (frames)
tSig = 1;
% Lucas-Kanade neighborhood size (voxels)
wSig = 4;

%%%% Loop through folders
for kk = 1:length(toProcess)

    imDir = toProcess{kk};
    disp(['FOLDER: ' imDir])
    
    % Get the list of tifs in this folder
    tifList = dir([imDir filesep '*.tif']);
    
    for jj = 1:length(tifList) % Loop through all tifs
        
        imName = tifList(jj).name(1:end-4);
        disp(['FILE: ' imName])

        % Run the optical flow
        process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

        disp(' ') % Formatting for the command window

    end

    disp(' ') % Formatting for the command window

end






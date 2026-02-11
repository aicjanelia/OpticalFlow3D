% This script runs the optical flow over the data presented in the
% publication. It does not create visualizations or run any analyses.

% Start with a clean workspace
clc, clear, close all

%%% Set up optical flow parameters
% Spatial smoothing (voxels)
xyzSig = 3;
% Temporal smoothing (frames)
tSig = 1;
% Lucas-Kanade neighborhood size (voxels)
wSig = 4;

%% 2D data in the 'OneTif' format

%%%% Set up information about the experiment files
% List all folders to be processed.
% All .tif files in this folder will be processed.
toProcess = {'X:\Force Project\PublicationData\SpinningDisk_Myosin\2D_cropped',...
    'X:\Force Project\PublicationData\SpinningDisk_Myosin\2D_yslice1175'};
% There is one tif for the entire time lapse, so this is a "OneTif" type movie
fileType = 'OneTif';
% This is 2D data (a single slice)
spatialDimensions = 2;

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

%% 2D data in the 'SequenceT' format

%%%% Set up information about the experiment files
% Folder containing sequence of tif files.
imDir = 'X:\Force Project\PublicationData\SpinningDisk_Myosin\2D_zslice7'; 
% File name of images of interest. Insert * as wildcards for values that
% change across the sequence.
imName = '20241107_U2OS_SGRLC_100Xoil_15mintimelapse_01_25plaser-slice7_t*'; 
% There is one tif per timepoint, so this is a "SequenceT" type movie
fileType = 'SequenceT';
% This is 2D data (a single slice)
spatialDimensions = 2;

disp(['FOLDER: ' imDir])
disp(['FILE: ' imName])

%%%% Run the optical flow
% process_flow takes care of parsing file names, etc.
% It will display progress in the Command Window as it progresses
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

disp(' ')

%% 3D data in the 'OneTif' format

%%%% Set up information about the experiment files
% List all folders to be processed.
% All .tif files in this folder will be processed.
toProcess = {'X:\Force Project\PublicationData\SpinningDisk_Myosin\3D_cropped',...
    'X:\Force Project\PublicationData\SpinningDisk_Myosin\3D_timelapse'};
% There is one tif for the entire time lapse, so this is a "OneTif" type movie
fileType = 'OneTif';
% This is 3D data (a z-stack)
spatialDimensions = 3;

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

%% 3D data in the 'SequenceT' format

%%%% Set up information about the experiment files
% Folder containing sequence of tif files.
imDir = 'X:\Force Project\PublicationData\MOSAIC_Actin\deskew_after_decon'; 
% File name of images of interest. Insert * as wildcards for values that
% change across the sequence.
imName = 'scan_Cam1_ch0_tile0_t*_deskew_after_decon'; 
% There is one tif per timepoint, so this is a "SequenceT" type movie
fileType = 'SequenceT';
% This is 3D data (a z-stack)
spatialDimensions = 3;

disp(['FOLDER: ' imDir])
disp(['FILE: ' imName])

%%%% Run the optical flow
% process_flow takes care of parsing file names, etc.
% It will display progress in the Command Window as it progresses
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

disp(' ')

%% 3D data in the 'SequenceT' format - multichannel

%%%% Set up information about the experiment files
% Folder containing sequence of tif files.
imDir = 'X:\Force Project\PublicationData\LLSM_twoChannels\scan6'; 
% File name of images of interest. Insert * as wildcards for values that
% change across the sequence.
imName = 'scan_CamA_ch0_CAM1_stack*_488nm_*msec_*msecAbs_000x_000y_000z_0000t_decon'; 
% There is one tif per timepoint, so this is a "SequenceT" type movie
fileType = 'SequenceT';
% This is 3D data (a z-stack)
spatialDimensions = 3;

disp(['FOLDER: ' imDir])
disp(['FILE: ' imName])

%%%% Run the optical flow
% process_flow takes care of parsing file names, etc.
% It will display progress in the Command Window as it progresses
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

disp(' ')

% NOW THE SECOND CHANNEL
imName = 'scan_CamA_ch1_CAM1_stack*_560nm_*msec_*msecAbs_000x_000y_000z_0000t_decon'; 
disp(['FOLDER: ' imDir])
disp(['FILE: ' imName])
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

disp(' ')

%% 3D data in the 'SequenceT' format - multichannel

%%%% Set up information about the experiment files
% Folder containing sequence of tif files.
imDir = 'X:\Force Project\PublicationData\LLSM_twoChannels\scan1'; 
% File name of images of interest. Insert * as wildcards for values that
% change across the sequence.
imName = 'scan_CamA_ch0_CAM1_stack*_488nm_*msec_*msecAbs_000x_000y_000z_0000t_decon'; 
% There is one tif per timepoint, so this is a "SequenceT" type movie
fileType = 'SequenceT';
% This is 3D data (a z-stack)
spatialDimensions = 3;

disp(['FOLDER: ' imDir])
disp(['FILE: ' imName])

%%%% Run the optical flow
% process_flow takes care of parsing file names, etc.
% It will display progress in the Command Window as it progresses
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

disp(' ')

% NOW THE SECOND CHANNEL
imName = 'scan_CamA_ch1_CAM1_stack*_560nm_*msec_*msecAbs_000x_000y_000z_0000t_decon'; 
disp(['FOLDER: ' imDir])
disp(['FILE: ' imName])
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

disp(' ')

%% Data with artificial translation

%%%% Set up information about the experiment files
% Folder containing sequence of tif files.
imDir = 'X:\Force Project\PublicationData\SpinningDisk_Myosin\3D_translationOnly'; 
% File name of images of interest. Insert * as wildcards for values that
% change across the sequence.
imName = 'manualTranslation_t*'; 
% There is one tif per timepoint, so this is a "SequenceT" type movie
fileType = 'SequenceT';
% This is 3D data (a z-stack)
spatialDimensions = 3;

disp(['FOLDER: ' imDir])
disp(['FILE: ' imName])

%%%% Run the optical flow
% process_flow takes care of parsing file names, etc.
% It will display progress in the Command Window as it progresses
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

%% 3D data in the 'OneTif' format - SoRa

%%%% Set up information about the experiment files
% List all folders to be processed.
% All .tif files in this folder will be processed.
toProcess = {'X:\Force Project\PublicationData\SoRa_twoChannels'};
% There is one tif for the entire time lapse, so this is a "OneTif" type movie
fileType = 'OneTif';
% This is 3D data (a z-stack)
spatialDimensions = 3;

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

%% 2D data in the 'OneTif' format - SoRa

%%%% Set up information about the experiment files
% List all folders to be processed.
% All .tif files in this folder will be processed.
toProcess = {'X:\Force Project\PublicationData\SoRa_twoChannels\2Dslice'};
% There is one tif for the entire time lapse, so this is a "OneTif" type movie
fileType = 'OneTif';
% This is 2D data (a single slice)
spatialDimensions = 2;


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

%% 2D data in the 'OneTif' format - Chad

%%%% Set up information about the experiment files
% List all folders to be processed.
% All .tif files in this folder will be processed.
toProcess = {'X:\Force Project\ChadDiffusionData'};
% There is one tif for the entire time lapse, so this is a "OneTif" type movie
fileType = 'OneTif';
% This is 2D data (a single slice)
spatialDimensions = 2;


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

%% 3D data in the 'SequenceT' format

%%%% Set up information about the experiment files
% Folder containing sequence of tif files.
imDir = 'Y:\instruments\ForceProject\MOSAIC_deskew_after_decon\fused'; 
% File name of images of interest. Insert * as wildcards for values that
% change across the sequence.
imName = 'RotatedCropped_fused_tp_*_ch_0'; 
% There is one tif per timepoint, so this is a "SequenceT" type movie
fileType = 'SequenceT';
% This is 3D data (a z-stack)
spatialDimensions = 3;


disp(['FOLDER: ' imDir])
disp(['FILE: ' imName])

%%%% Run the optical flow
% process_flow takes care of parsing file names, etc.
% It will display progress in the Command Window as it progresses
process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)

disp(' ')


%%
disp('Processing Complete')
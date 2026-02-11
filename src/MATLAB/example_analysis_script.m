% This script assumes that process_flow has already been run on the
% datasets described in the Paths section. See example_processing_script if
% you have not yet calculated the raw optical flow field.

% Start with a clear workspace
clc
clear
close all

%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set up information about the experiment files
% Folder containing sequence of tif files.
imDir = 'C:\fullpath\to\image\files\';
% What is the name of the image sequence you are analyzing?
% If it was a SequenceT dataset, list your regular expression with all * removed.
% (This is also the name of the subfolder inside OpticalFlow3D)
flowName = '20241107_U2OS_SGRLC_100Xoil_15mintimelapse_01_25plaser';
% This is the directory the flow output is saved in. If you have not
% rearranged the folder structure since processing, this is correct.
flowDir = [imDir filesep 'OpticalFlow3D' filesep flowName];
% Where should any analysis figures or visualizations be saved?
% This default will make a subfolder next to the images.
saveDir = [imDir '\AnalysisOutput\'];

%%%% Metadata
xyscale = 0.0425177; % um/pixel
zscale = 0.3; % um/pixel
tscale = 30/60; % minutes/frame
tRange = 4:18; % which frames should be analyzed? Indexing begins at 1.

%%%% Analysis Parameters
% Reliability thresholding removes background regions and uncertain flow
% vectors. In this example a percentile threshold us used, but other
% thresholding methods could be applied as appropriate!
relPer = 90; % Percentile threshold for reliability

%%%% Visualization Parameters
% Vectors are plotted using "quiver". These parameters control appearance.
gap = 2; % Plot every gap-th vector to keep plot from being overcrowded.
qscale = 1; % 1 or 0 plot vectors at their exact scale, other values scale 
% the vectors appropriately (e.g., 2 would show the vectors at twice their
% real length).
zSlice = 7; % z-slice to use for 2D visualizations

%% START PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up

disp(datetime('now'))

% If the directory for saving does not already exist, make it.
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

% Get Image Sizes to format plots appropriately
relFile = [flowName '_rel_t' num2str(4,'%04u') '.tiff'];
meta = imfinfo([flowDir filesep relFile]);
Nx = meta(1).Width;
Ny = meta(2).Height;
Nz = length(meta);

% Set up matrices to help with the quiver plots
x = (0:Nx-1)*xyscale;
y = (0:Ny-1)*xyscale;
z = (0:Nz-1)*zscale;
[X,Y,Z] = meshgrid(x,y,z);
clear x y z
% These will be useful if gap ~=1
xD = X(1:gap:end);
yD = Y(1:gap:end);
zD = Z(1:gap:end);

%% EXAMPLE AT ONE TIME POINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exampleFrame = 4;

%% Reliability thresholding
    
% Get the appropriate file name for this time point
relFile = [flowName '_rel_t' num2str(exampleFrame,'%04u') '.tiff'];
% Load the reliability tif
rel = TIFFvolume([flowDir filesep relFile],Nz);
% Determine the threshold
relThresh = prctile(rel(:),relPer);
% Apply the threshold
relMask = rel > relThresh;

figure(1)
imshow(relMask(:,:,zSlice))
title('Reliability Mask')
disp(['Reliability threshold = ' num2str(relThresh)])

%% Load in velocities
% Each file is first loaded in, then masked by the reliability
% thresholding. Values are converted from pixels/frame (as they are
% saved in the tif files) to real units. This conversion is important
% to account for any anisotropy in z.

vx = TIFFvolume([flowDir filesep flowName '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
vx = vx.*relMask./relMask*xyscale/tscale;

vy = TIFFvolume([flowDir filesep flowName '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
vy = vy.*relMask./relMask*xyscale/tscale;

vz = TIFFvolume([flowDir filesep flowName '_vz_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
vz = vz.*relMask./relMask*zscale/tscale;

%% Useful Quantities
% Magnitude, theta, and phi are useful for characterizing the flow.
% They could be plotted as quiver plots, saved as their own .tif files,
% or can be used for statistical measurements.

mag = sqrt(vx.^2 + vy.^2 + vz.^2);
theta = atan2(vy,vx);
phi = atan(vz./sqrt(vx.^2+vy.^2));

%% Example Quiver Figure

% Here the visualization will be based on the 2nd z-slice
vD = theta(:,:,zSlice);

% Apply the optional gap parameter to clean up the visualization
vxD = vx(1:gap:end);
vyD = vy(1:gap:end);    
vD = vD(1:gap:end);

% Color is based on binned values of theta. Here 15 bins are used, and
% a circular colormap is chosen from colorcet. See https://colorcet.com/gallery.html
bins = linspace(-pi-eps,pi+eps,15);
cmap = colorcet('C8');
cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

figure(2)
set(gcf,'Position',[100 100 560*2 420*2])
% Loop through the bins to set the color appropriately
for jj = 1:length(bins)-1     
    slice = (vD>=bins(jj)) & (vD<bins(jj+1));
    j = quiver(xD(slice),yD(slice),vxD(slice),vyD(slice),0,...
        'Color',cmap(jj,:),'LineWidth',0.25);
    hold on
    % Set vector scaling as necessary
    if qscale ~= 0
        hU = get(j,'UData') ;
        hV = get(j,'VData') ;
        set(j,'UData',qscale*hU,'VData',qscale*hV)
    end
end
hold off
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'ydir','reverse')
set(gca,'Color','k')
set(gca,'XTick',[],'YTick',[])
% This sets the view to the full image; other axis limits could be
% appropriate.
xlim([0 (Nx-1)*xyscale])
ylim([0 (Ny-1)*xyscale])

% Draw the figure before saving
drawnow;
% Save the figure
quivA = getframe(gca);
imwrite(quivA.cdata,[saveDir filesep 'ExampleThetaQuiver.png'])

%% Statistics: Magnitude Distribution

binEdges = linspace(0,0.5,30);

figure(3)
histogram(mag(:),binEdges,'Normalization','Probability')
xlabel(['Magnitude (' char(181) 'm/min)'])
ylabel('Probability')
box off
set(gca,'FontSize',16)

saveas(gcf,[saveDir filesep 'ExampleMagnitudeDistribution.png'])

%% Statistics: Theta Distribution
% Theta by default has a positive value when an object moves from the top of the image to the bottom.
% This is a common convention in image processing, but it can be confusing with our intuition that up is positive.
% By plotting negative theta, positive 90 degrees becomes upward motion.

figure(4)
polarhistogram(-theta(:),20,'Normalization','Probability')
set(gca,'FontSize',16)
title('XY Direction of Motion')
set(gca,'ThetaTick',0:45:360)


saveas(gcf,[saveDir filesep 'ExampleDirectionDistribution.png'])

%% MEASUREMENTS OVER TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This example will keep track of mean magnitude over time. This variable
% preallocates the storage for the mean measurement.
meanMag = NaN*ones(max(tRange)+3,1); % Exlcuded frames on the edges of timelapse will have value NaN.

%% Analyze each time point by looping through them
for exampleFrame = tRange

    disp([char(datetime('now')) ' - Processing frame ' num2str(exampleFrame) '...'])

    %% Reliability Thresholding
    
    % Get the appropriate file name for this time point
    relFile = [flowName '_rel_t' num2str(exampleFrame,'%04u') '.tiff'];
    % Load the reliability tif
    rel = TIFFvolume([flowDir filesep relFile],Nz);
    % Determine the threshold
    relThresh = prctile(rel(:),relPer);
    % Apply the threshold
    relMask = rel > relThresh;

    %% Load in velocities
    % Each file is first loaded in, then masked by the reliability
    % thresholding. Values are converted from pixels/frame (as they are
    % saved in the tif files) to real units. This conversion is important
    % to account for any anisotropy in z.
    
    vx = TIFFvolume([flowDir filesep flowName '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vx = vx.*relMask./relMask*xyscale/tscale;

    vy = TIFFvolume([flowDir filesep flowName '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vy = vy.*relMask./relMask*xyscale/tscale;

    vz = TIFFvolume([flowDir filesep flowName '_vz_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vz = vz.*relMask./relMask*zscale/tscale;

    %% Useful Quantities
    % Magnitude, theta, and phi are useful for characterizing the flow.
    % They could be plotted as quiver plots, saved as their own .tif files,
    % or can be used for statistical measurements.

    mag = sqrt(vx.^2 + vy.^2 + vz.^2);
    theta = atan2(vy,vx);
    phi = atan(vz./sqrt(vx.^2+vy.^2));

    % As an example, here we measure the mean magnitude over time.
    meanMag(exampleFrame) = mean(mag(~isnan(mag)));

    %% Example 2D Theta Quiver

    % Here the visualization will be based on the 2nd z-slice
    vD = theta(:,:,2);

    % Apply the optional gap parameter to clean up the visualization
    vxD = vx(1:gap:end);
    vyD = vy(1:gap:end);    
    vD = vD(1:gap:end);

    % Color is based on binned values of theta. Here 15 bins are used, and
    % a circular colormap is chosen from colorcet. See https://colorcet.com/gallery.html
    bins = linspace(-pi-eps,pi+eps,15);
    cmap = colorcet('C8');
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    figure(5)
    set(gcf,'Position',[100 100 560*2 420*2])
    % Loop through the bins to set the color appropriately
    for jj = 1:length(bins)-1     
        slice = (vD>=bins(jj)) & (vD<bins(jj+1));
        j = quiver(xD(slice),yD(slice),vxD(slice),vyD(slice),0,...
            'Color',cmap(jj,:),'LineWidth',0.25);
        hold on
        % Set vector scaling as necessary
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
    end
    hold off
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'ydir','reverse')
    set(gca,'Color','k')
    set(gca,'XTick',[],'YTick',[])
    % This sets the view to the full image; other axis limits could be
    % appropriate.
    xlim([0 (Nx-1)*xyscale])
    ylim([0 (Ny-1)*xyscale])

    % Draw the figure before saving
    drawnow;
    % Save the figure
    quivA = getframe(gca);
    imwrite(quivA.cdata,[saveDir filesep 'ThetaQuiver_Frame' num2str(exampleFrame,'%04u') '.png'])

    %% Optional Clean up

    % close all
    % disp([char(datetime('now')) ' - Frame ' num2str(exampleFrame) ' complete'])

end

%% Plot mean magnitude

tVec = 0:(max(tRange)+2);
tVec = tVec*tscale; % Put in minutes

figure(6)
plot(tVec,meanMag,'-o','LineWidth',2)
xlabel('Time (min)')
ylabel(['Mean Magnitude (' char(181) 'm/min)'])
set(gca,'FontSize',16)
ylim([0 max(meanMag)+0.05]) % Magnitude is bounded by zero
box off

saveas(gcf,[saveDir filesep 'ExampleMagnitudeVsTime.png'])

%% Clean up
disp('Processing Complete')
% close all
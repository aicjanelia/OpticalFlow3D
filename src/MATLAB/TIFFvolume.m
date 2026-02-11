function vol = TIFFvolume(filename,frames,startFrame)
% Use TIFF or fread to load in 3D volumes
%
% By default, TIFF files are read using the TIFF library.
% If a >4 GB ImageJ file is detected, the function switches to using fread
% instead as TIFF only sees the first xy plane of the file.
%
% vol = TIFFvolume(filename,Nx,Ny,frames)
%
% INPUTS:
% filename  = path for the saved file
% frames    = Number of frames from the volume to load. This can be time
%               and/or z depending on the input file.
% startFrame = First frame to load from. Allows for only loading part of
%               the image. Frames startFrame:startFrame+frames are loaded.
%               Default value = 1 (the first frame).
%
% OUTPUTS:
% vol       = Matrix of size [Ny,Nx,frames] and type double
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Use default inputs.
if ~exist('startFrame','var') || isempty(startFrame)
    startFrame = 1;
end
if startFrame < 1
    error('TIFFvolume requires the startFrame >=1 (this function uses 1-based indexing)')
end
type = 'tiff';

%%%% Check the metadata
meta = imfinfo(filename,'tif');
Nx = meta(1).Width;
Ny = meta(1).Height;
Nframes = length(meta);
bitD = meta(1).BitDepth;
if bitD == 64
    bitD = 'double';
elseif bitD == 32
    bitD = 'single';
elseif bitD == 16 || bitD == 8
    bitD = ['uint' num2str(bitD)];
else
    error('Unexpected bit depth')
end

if Nframes==1 && (frames+startFrame-1)>1
    imagej = isfield(meta,'ImageDescription');
    if imagej
        type = 'fread';
    else
        if startFrame+frames-1 > Nframes
            error(['startFrame is ' num2str(startFrame) ' and ' num2str(frames) ' frames were requested (requiring ' num2str(startFrame+frames) ' frames total) but only ' num2str(Nframes) ' frames were detected'])
        end
    end
end

%%%% Load in the files
if strcmp(type,'tiff')
    %  Open the file for reading
    intiff = Tiff(filename,"r");
    
    % If not starting from 1, progress to the correct directory
    if startFrame > 1
        for ii = 1:startFrame-1
            nextDirectory(intiff);
        end
    end
    
    % Load the number of frames requested
    vol = zeros(Ny,Nx,frames,bitD);
    for ii = 1:frames-1
        vol(:,:,ii) = read(intiff);
        nextDirectory(intiff);
    end
    vol(:,:,frames) = read(intiff);
    
    % Close the file
    intiff.close();

elseif strcmp(type,'fread') 
% This section is heavily inspired by https://www.mathworks.com/matlabcentral/fileexchange/61376-imread_big-read-in-tiff-stacks-larger-than-4gb

    offsets = meta(1).StripOffsets;
    byteCounts = meta(1).StripByteCounts;
    machinefmt = meta(1).ByteOrder;
    if strcmp(machinefmt,'big-endian')
        machinefmt = 'b';
    elseif strcmp(machinefmt,'little-endian')
        machinefmt ='l';
    else
        error('Unexpected order for reading byes. See https://www.mathworks.com/help/matlab/ref/fopen.html#btrnibn-1-machinefmt for options.')
    end

    if strcmp(bitD(1:4),'uint')
        uinttype = ['uint' bitD(5:end) '=>uint' bitD(5:end)];
    elseif strcmp(bitD,'double')
        uinttype = 'double=>double';
    elseif strcmp(bitD,'single')
        uinttype = 'single=>single';
    else
        error('Unexpected bit depth')
    end

    intiff = fopen(filename,'r');
    vol = zeros(Ny,Nx,frames,bitD);
    start_point = offsets(1) + (startFrame-1)*byteCounts(1); % offset for first requested frame
    fseek(intiff, start_point, 'bof'); % Move to the offset for the first requested frame
    for ii = 1:frames
        A = fread(intiff,[Nx Ny],uinttype,machinefmt);
        vol(:,:,ii) = A';
    end

    fclose(intiff);

else
    error('Unexpected read type')
end
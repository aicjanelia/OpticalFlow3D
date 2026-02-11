function TIFFwrite(filename,A)
% Use TIFF to save float tifs
%
% TIFF files are written using the bigtiff option to allow for any output
% file size. They are saved as 64-bit with LZW compression.
%
% USAGE: TIFFwrite(filename,A)
%
% INPUTS:
% filename  = path for the saved file
% A         = 2D or 3D matrix to save as a .tif
%
% OUTPUTS:
% A .tiff file is saved at the path specified by filename
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TIFF Tags
tagstruct.ImageLength = size(A,1);
tagstruct.ImageWidth = size(A,2);
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.BitsPerSample = 64;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.Compression = Tiff.Compression.LZW;

outtiff = Tiff(filename,'w8'); % Using bigtiff writing
for ii = 1:size(A,3)-1
    outtiff.setTag(tagstruct);
    outtiff.write(A(:,:,ii));
    outtiff.writeDirectory();
end
outtiff.setTag(tagstruct);
outtiff.write(A(:,:,size(A,3)));

% Close the file
outtiff.close();
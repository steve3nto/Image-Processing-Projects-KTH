function Y = yuv_import_y(filename,dims,numfrm,startfrm)

%Imports luminance component of YUV sequence into Matlab
%
%Y = yuv_import_y(filename,dims,numfrm,startfrm)
%
%Input:
% filename - YUV sequence file
% dims - dimensions of the frame [width height]
% numfrm - number of frames to read
% startfrm - [optional, default = 0] specifies from which frame to start reading
%            with the convention that the first frame of the sequence is denoted
%            with 0
%
%Output:
% Y - Cell with numfrm elements, where each element is the luminance
% component of a frame
%
%Example:
% Read 50 first frames from file
% Y = yuv_import_y('mother-daughter_qcif.yuv',[176 144],50);

fid=fopen(filename,'rb');
if (fid < 0) 
    error('File does not exist!');
end;

Y = cell(numfrm,1);
frelem = prod(dims)*3/2;

if (nargin == 4) %go to the starting frame
    fseek(fid, startfrm * frelem , 0);
else
    startfrm = 0;
end

for i=1:numfrm
    Y{i} = fread(fid,[dims(1) dims(2)],'uint8')';
    fseek( fid, (startfrm+i)*frelem, -1 );
end
fclose(fid);

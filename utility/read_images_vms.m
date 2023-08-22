function imagearray = read_images_vms(filename, rownum, columnnum, fileformat, imagenum)
%READ_IMAGES_VMS read a sequence of images from a data file, stores the images in a 3-D 
%	array rownum*columnnum*imagenum
%Input:
%	filename   -- string, name of the data file; e.g. 'myFile.rec'
%	rownum     -- integer, pixel number in the vertical direction
%	columnnum  -- integer, pixel number in the horizonal direction
%	fileformat -- format of the data file, e.g. 'int16', 'int32', 'double'
%	imagenum   -- number of images contained in the file
%Output:
%	imagearray -- rownum*columnnum*imagenum array, double type, stores the image gray scale
%Hanzhang Lu
%Date: 06/30/00

%open file
fileid = fopen (filename, 'r', 'ieee-le');
%read data into a dummy variable
dummy = fread( fileid, rownum*columnnum*imagenum, fileformat);
%copy the dummy variable into output array
imagearray = reshape (dummy, [rownum columnnum imagenum]);
%close the file
fclose (fileid);
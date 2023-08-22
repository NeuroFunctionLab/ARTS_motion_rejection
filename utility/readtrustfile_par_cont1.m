function varargout=readtrustfile_par_cont1(filename, dynamic,matrix)
% read in .rec images based on .par parameter
% No need specify the dynamic image number
% Input:
% filename: no extension, e.g.
% filename='C:\data\trustsampledata\3T0603110\3T603110_20_1';
% dynamic: number of dynamic scan ran in each eTe

% Output:
% varargout{1}=aslimg
% Feng Xu 
%Date: 10/29/2007


totaldyn=matrix(3);
rownum=matrix(1);columnnum=matrix(2);
tempstr=strcat(filename,'.REC');
aslimg=read_images_vms(tempstr,rownum,columnnum, 'int16',totaldyn);  

aslimg=aslimg;

nte=totaldyn/dynamic; % m=nte

aslimgfinal=aslimg;
aslimgodd=aslimgfinal(:,:,1:2:totaldyn);
aslimgeven=aslimgfinal(:,:,2:2:totaldyn);clear aslimgfinal;
asldiff=aslimgeven-aslimgodd;
varargout{1}=aslimg;
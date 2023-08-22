function varargout=readtrustfile_par_cont(filename, dynamic)
% read in .rec images based on .par parameter
% No need specify the dynamic image number
% Input:
% filename: no extension, e.g.
% filename='C:\data\trustsampledata\3T0603110\3T603110_20_1';
% dynamic: number of dynamic scan ran in each eTe

% Output:
% varargout{1}=con_all;
% varargout{2}=lab_all;
% varargout{3}=matrix;
% varargout{4}=dif_all;
% varargout{5}=con_mean;
% varargout{6}=lab_mean;
% varargout{7}=dif_mean;

% Feng Xu 
%Date: 10/29/2007

tempstr=strcat(filename,'.PAR');
[scale, matrix]=read_PAR(tempstr);
RI=scale(1);RS=scale(2);SS=scale(3);
totaldyn=matrix(3);
rownum=matrix(1);columnnum=matrix(2);
tempstr=strcat(filename,'.REC');
aslimg=read_images_vms(tempstr,rownum,columnnum, 'int16',totaldyn);   
aslimg=(aslimg.*RS+RI*ones(size(aslimg)))/RS/SS;

nte=totaldyn/dynamic;

aslimgfinal=aslimg; clear aslimg;
aslimgodd=aslimgfinal(:,:,1:2:totaldyn);
aslimgeven=aslimgfinal(:,:,2:2:totaldyn);clear aslimgfinal;
asldiff=aslimgeven-aslimgodd;

con_all=reshape(aslimgeven,rownum,columnnum,nte,dynamic/2);
lab_all=reshape(aslimgodd,rownum,columnnum,nte,dynamic/2);
dif_all=reshape(asldiff,rownum,columnnum,nte,dynamic/2);

con_mean=mean(con_all,4); 
lab_mean=mean(lab_all,4);
dif_mean=mean(dif_all,4);
% 
% for i=1:nte
%     t=dynamic/2*(i-1)+1;
% lab_mean(:,:,i)=mean(aslimgodd(:,:,t:t+dynamic/2-1),3);
% con_mean(:,:,i)=mean(aslimgeven(:,:,t:t+dynamic/2-1),3);
% dif_mean(:,:,i)=mean(asldiff(:,:,t:t+dynamic/2-1),3);
% con_all(:,:,i,:)=aslimgeven(:,:,t:t+dynamic/2-1);
% lab_all(:,:,i,:)=aslimgodd(:,:,t:t+dynamic/2-1);
% dif_all(:,:,i,:)=asldiff(:,:,t:t+dynamic/2-1);
% end

%output
varargout{1}=con_all;
varargout{2}=lab_all;
varargout{3}=matrix;
varargout{4}=dif_all;
varargout{5}=con_mean;
varargout{6}=lab_mean;
varargout{7}=dif_mean;
function [t2,ci,Yv]=trust5combine_SS_4te_6dyn_GUI_linux(filename1,matrix,seq,hct)
% Script edition with alignment
%***TRUST acquired based on the order: TE0,TE0,TE0, TE40,TE40,TE40,
%TE80,TE80,TE80, TE160,TE160,TE160
% Input:
% filename: .REC file directory
% nlist: pixel number of interest
% roiflag: dynamic drawing roi for different TE set '1' or static drawing roi from ...
%  first TE set '0'
% datastructure: data structure N,N,repeat*2 (include label and control, # of total dynamics each TE weighted scan has)
% nrep: n out of 4 repeat measurements want to use

% Output:
% T2 value in milisecond
% CI 95%
% std of repeated measurement
% SNR of difference image
% 
%%Feng Xu
%%09/21/2007 

%  filename='C:\data\Hypercapnia_timecourse\3T0819\3T0813_1'; % room air SS
% close all
% [fileNamePar,pathNamePar] = uigetfile('*.REC','Choose the REC image:','MultiSelect','off');
[fname_path fname_body fname_ext] = fileparts(filename1);
filename1=fullfile(fname_path,fname_body);
filename = [fname_path filesep fname_body];
etetemp=unique(seq);
etelen=length(etetemp);
ktemp=find(seq==1);
klen=length(ktemp);
ltemp=find(seq==40);
llen=length(ltemp);
mtemp=find(seq==80);
mlen=length(mtemp);
ntemp=find(seq==160);
nlen=length(ntemp);
roiflag=0;
nlist=4;
nrep=3;
bloodt1=1624;
ti=1200;
dynnum=matrix(3)/etelen;
if(matrix(3)~=length(seq))
    matrix(3)=length(seq);
end;
te=etetemp';
% matrix=zeros(1,2);
[noalign]=readtrustfile_par_cont1(filename,dynnum,matrix);
clear con_all
clear lab_all
nte=matrix(3)/dynnum;
te=te(1:nte)
% noalign=cat(4,con_all,lab_all);
columnnum=matrix(1);
rownum=matrix(2);

% % -------------do alignment-------------------------

rowspacing=3.438;
columnspacing=3.438;
slicespacing=5;
temp=filename;
%matrix
% g=exist([filename, '1'])
%seq=[1 1 1 1 1 1 40 40 40 40 40 40 80 80 80 80 80 80 160 160 160 160 160 160];
% seq=[1 1 40 40 80 80 160 160 1 1 40 40 80 80 160 160 1 1 40 40 80 80 160 160];
k=0;l=0;m=0;n=0;
if ~exist([filename, '_1'])
for i=1:length(seq)
    seq(i);
    switch seq(i)
        case 1
            tvar=1;
            k=k+1;
            if ~exist([filename, '_1'])
            tempstr=strcat(temp,'_',int2str(tvar));
            mkdir(tempstr);
            j=k;
            else
                tempstr=strcat(temp,'_',int2str(tvar));
                j=k;
            end;
         case 40
            tvar=2;
            l=l+1;
            if ~exist([filename, '_2'])
            tempstr=strcat(temp,'_',int2str(tvar));
            mkdir(tempstr);
            j=l;
            else
                tempstr=strcat(temp,'_',int2str(tvar));
                j=l;
            end;
         case 80
            tvar=3;
            m=m+1;
            if ~exist([filename, '_3'])
            tempstr=strcat(temp,'_',int2str(tvar));
            mkdir(tempstr);
            j=m;
            else
                tempstr=strcat(temp,'_',int2str(tvar));
                j=m;
            end;
        otherwise
            tvar=4;
            n=n+1;
            if ~exist([filename, '_4'])
            tempstr=strcat(temp,'_',int2str(tvar));
            mkdir(tempstr);
            j=n;
            else
                tempstr=strcat(temp,'_',int2str(tvar));
                j=n;
            end;
    end;

    if j<10
        a=strcat('00',int2str(j));
    elseif j<100
        a=strcat('0',int2str(j));
    else
        a=int2str(j);
    end
    tempstr1=strcat(tempstr,'\img_', a);
    write_ANALYZE(noalign(:,:,i),tempstr1,[64 64 1], [3.438 3.438 5], 1,16);

ct=[k l m n];
    if(k==klen||l==llen||m==mlen||n==nlen)
        bt=find(ct==6);
        switch bt %#ok<ALIGN>
            case 1
                k=0;
     
        case 2
            l=0;

    case 3
        m=0;

otherwise
    n=0;
end;
       tempstr=strcat(filename,'_',int2str(bt)); 
    %realign images using the first dynamic of first eTE as template
    %the images are realigned to the first dynamic in each scan
    target=strcat(strcat(tempstr,'\img_002.img'));
    
    for count=1:dynnum
        
    if count<10
        a=strcat('00',int2str(count));
    elseif count<100
        a=strcat('0',int2str(count));
    else
        a=int2str(count);
    end
    source=strcat(tempstr,'\img_', a, '.img');
    delete_if_exist(['img_' a '.mat']);
    other=source;
    spm_defaults;
    defs=defaults.realign;
    FlagsC = struct('quality',defs.estimate.quality,'fwhm',5,'rtm',0);
    spm_realign([target;source],FlagsC);
    FlagsR = struct('interp',defs.write.interp,...
		'wrap',defs.write.wrap,...
		'mask',defs.write.mask,...
		'which',2,'mean',1);
    FlagsR.which = 2; FlagsR.mean = 0;
    spm_reslice([target;source],FlagsR);
    warning off all
%     tvar
%     thisimg=read_images_vms(strcat(tempstr,'\rimg_', a, '.img'), rownum, columnnum, 'float',1);
    thisimg=loadimage(strcat(tempstr,'\rimg_', a, '.img'), 16);
    img_all(:,:,tvar,count)=thisimg(:,:);
    end
    end;
end 

else 
    for i=1:nte
        tempstr=strcat(filename,'_',int2str(i));
        for count=1:dynnum
            if count<10
                a=strcat('00',int2str(count));
            elseif count<100
                a=strcat('0',int2str(count));
            else
                a=int2str(count);
            end 
%             thisimg=read_images_vms(strcat(tempstr,'\rimg_', a, '.img'), rownum, columnnum, 'float',1);
            thisimg=loadimage(strcat(tempstr,'\rimg_', a, '.img'), 16);
            img_all(:,:,i,count)=thisimg(:,:);
        end
    end  
end;
con_all=img_all(:,:,:,2:2:dynnum);
lab_all=img_all(:,:,:,1:2:dynnum);
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dif_all=con_all-lab_all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
con_mean=mean(con_all,4);
lab_mean=mean(lab_all,4);
dif_mean=mean(dif_all,4);

%clear img_all thisimg noalign

mask=zeros(rownum,columnnum,nte);

for t=1:nte
    % draw ROI
    anat=dif_mean(:,:,t);  
    
    droi=zeros(rownum, columnnum);
    himage=imshow(anat,[min(anat(:)) max(anat(:))]);
    temp=roipoly();colorbar;
    droi=droi | temp;
    close;
    clear mylist;
    mycount=0;        
    for i=1:rownum
        for j=1:columnnum
            if droi(i,j)==1
                mycount=mycount+1;
                mylist(mycount,1)=i;
                mylist(mycount,2)=j;
                mylist(mycount,3)=anat(i,j,1);
            end
        end
    end
    [Y,I]=sort(mylist(:,3),'descend');
    if roiflag ==0
        mask=zeros(rownum,columnnum);
        for tt=1:nlist
            mask(mylist(I(tt),1),mylist(I(tt),2))=1;
        end
        break;
    else
        for tt=1:nlist
            mask(mylist(I(tt),1),mylist(I(tt),2),t)=1;
        end
    end
end

nn=ceil(sqrt(nte));
meanum=nrep*2;t2num=nte;
%fitting for all four repeated measurement
b_all=zeros(meanum/2,nte);
c_all=zeros(meanum/2,nte);
for i=1:nte
    temp3=dif_all(:,:,i,:);
    temp3=reshape(temp3,rownum,columnnum,meanum/2);
    if roiflag ==0
        temp4=average_actpixel1(1,mask,temp3);
%         clf;
        b_all(:,i)=temp4(:);
        temp3=con_all(:,:,i,:);
        temp3=reshape(temp3,rownum,columnnum,meanum/2);
        temp4=average_actpixel1(1,mask,temp3);
%         clf;
        c_all(:,i)=temp4(:);
    else
        temp4=average_actpixel1(1,mask(:,:,i),temp3);
%         clf; 
        b_all(:,i)=temp4(:);
    temp3=con_all(:,:,i,:);
    temp3=reshape(temp3,rownum,columnnum,meanum/2);
    temp4=average_actpixel1(1,mask(:,:,i),temp3);
%     clf;
    c_all(:,i)=temp4(:);
    end   
end
std_rep=std(b_all,0,1);
background=zeros(nte,1);
for i=1:nte
background(i)= mean(mean(abs(dif_mean(1:10,1:10,i))))+mean(mean(abs(dif_mean(rownum-9:rownum,1:10,i))))+...
    mean(mean(abs(dif_mean(1:10,columnnum-9:columnnum,i))))+mean(mean(abs(dif_mean(rownum-9:rownum,columnnum-9:columnnum,i))));
if background(i)==0
    background(i)=0.000001; % avoid the divisor equals zero
end
end
b=mean(b_all,1);b=reshape(b,nte,1);
snr_dif=b./background;

%%% T2 fitting by all b points

te_all=te(1:nte); temp=te(1:nte);
for i=1:meanum/2-1
    te_all=[te_all;temp];
end
b_all_vector=reshape(b_all',meanum/2*nte,1);
t2num=nte; 
c_all_vector=reshape(c_all',meanum/2*nte,1);
b_all_vector
c_all_vector
for i=1:meanum/2
    for j=1:nte-1
        b_all_vector((i-1)*t2num+j)=exp(-(ti-te(nte))/bloodt1)/exp(-(ti-te(j))/bloodt1)*b_all_vector((i-1)*t2num+j);
        c_all_vector((i-1)*t2num+j)=exp(-(ti-te(nte))/bloodt1)/exp(-(ti-te(j))/bloodt1)*c_all_vector((i-1)*t2num+j);
    end
end
b_all_vector
c_all_vector


[temp1,resid, jacob]=nlinfit_hlu(te_all/1000, b_all_vector,'monexp_model',[100,10]);
t2=1000./temp1(2);

conintval=nlparci(temp1,resid, jacob); %95% confidence interval for estimates
aa=1000./conintval;
ci=aa(2,:);

[temp1_c,resid_c, jacob_c]=nlinfit_hlu(te_all/1000, c_all_vector,'monexp_model',[100,10]);
conintval_c=nlparci(temp1_c,resid_c, jacob_c); %95% confidence interval for estimates
t2_c=1000./temp1_c(2);
aa_c=1000./conintval_c;
ci_c=aa_c(2,:);
Yv=prox_invitroblood_batch_R2vsyandhct_Golay2001_oneHct50_1(t2,hct);

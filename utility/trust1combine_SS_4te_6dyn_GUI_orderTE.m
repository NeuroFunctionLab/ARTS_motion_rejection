% Script edition with alignment
%***TRUST image acquired in the order TE0,TE40,TE80,TE160, TE0,TE40,...
% Input:
% filename: .REC file directory
% nlist: pixel number of interest
% foiflag: dynamic drawing roi for different TE set '1' or static drawing roi from ...
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

%  This code is used for hyperoxia study
% close all
% clear all

[fileNamePar,pathNamePar] = uigetfile('*.REC','Choose the REC image:','MultiSelect','off');
[fname_path fname_body fname_ext] = fileparts([pathNamePar fileNamePar]);

filename = [fname_path filesep fname_body];

roiflag=0;
nrep=3;
bloodt1=1624;
ti=1200;
nlist=4;
datastruct=6;
dynnum=6;

te=[1,40,80,160]'; 
matrix=zeros(1,2);
[con_all,lab_all,matrix]=readtrustfile_par_cont_orderTE(filename,dynnum);
nte=matrix(3)/datastruct;
te=te(1:nte);

noalign=cat(4,con_all,lab_all);
columnnum=matrix(1);
rownum=matrix(2);

% % -------------do alignment-------------------------

rowspacing=3.438;
columnspacing=3.438;
slicespacing=5;
temp=filename

if ~exist([filename, '1'])
for i=1:nte
    tempstr=strcat(temp,int2str(i));
    mkdir(tempstr);

    for j=1:dynnum
    if j<10
        a=strcat('00',int2str(j));
    elseif j<100
        a=strcat('0',int2str(j));
    else
        a=int2str(j);
    end
    fileid=fopen(strcat(tempstr,'\img_', a, '.img'),'w');
    
    fwrite(fileid,noalign(:,:,i,j),'int16');
    fclose(fileid);
    tempstr1=strcat(tempstr,'\img_', a);
    tempstr2=['U:\Data\TestData ', tempstr1,...
            ' ','2',' ',...
            num2str(rownum), ' ',...
            num2str(columnnum), ' ',...
            '1', ' ',...
            num2str(rowspacing), ' ',...
            num2str(columnspacing), ' ',...
            num2str(slicespacing)];
    dos(tempstr2);
    end
    %realign images using the first dynamic of first eTE as template
    %the images are realigned to the first dynamic in each scan
    target=strcat(strcat(tempstr,'\img_001.img'));
    for count=1:dynnum
    if count<10
        a=strcat('00',int2str(count));
    elseif count<100
        a=strcat('0',int2str(count));
    else
        a=int2str(count);
    end
    source=strcat(strcat(tempstr,'\img_', a, '.img'));
    other=source;
    spm_defaults;
    defs=defaults.realign
    FlagsC = struct('quality',defs.estimate.quality,'fwhm',5,'rtm',0);
    spm_realign([target;source],FlagsC);
    FlagsR = struct('interp',defs.write.interp,...
		'wrap',defs.write.wrap,...
		'mask',defs.write.mask,...
		'which',2,'mean',1);
    FlagsR.which = 2; FlagsR.mean = 0;
    spm_reslice([target;source],FlagsR);
    thisimg=read_images_vms(strcat(tempstr,'/rimg_', a, '.img'), rownum, columnnum, 'int16',1);
    img_all(:,:,i,count)=thisimg(:,:);
    end
end
else 
    for i=1:nte
        tempstr=strcat(filename,int2str(i));
        for count=1:dynnum
            if count<10
                a=strcat('00',int2str(count));
            elseif count<100
                a=strcat('0',int2str(count));
            else
                a=int2str(count);
            end 
            thisimg=read_images_vms(strcat(tempstr,'\rimg_', a, '.img'), rownum, columnnum, 'int16',1);
            img_all(:,:,i,count)=thisimg(:,:);
        end
    end  
end
con_all=img_all(:,:,:,1:dynnum/2);
lab_all=img_all(:,:,:,dynnum/2+1:dynnum);
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%
dif_all=con_all-lab_all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
con_mean=mean(con_all,4);
lab_mean=mean(lab_all,4);
dif_mean=mean(dif_all,4);

%clear img_all thisimg noalign

mask=zeros(rownum,columnnum,nte);
con_all=con_all(:,:,:,1:nrep);
lab_all=lab_all(:,:,:,1:nrep);
dif_all=dif_all(:,:,:,1:nrep);

show_imgs(dif_mean,1:nte,min(dif_mean(:)),max(dif_mean(:)));
show_imgs(cat(3,lab_mean,con_mean),1:2*nte,min(lab_mean(:)),max(con_mean(:)));

% get the top value mask
for t=1:nte
        % draw ROI
        anat=dif_mean(:,:,t); %susanxu 03/18
%         anat=con_mean(:,:,3); %susanxu 03/18
        droi=zeros(rownum, columnnum);
        reply_end='Y';
        while reply_end ~= 'N'            
            figure(100);
            imshow(anat,[min(anat(:)) max(anat(:))]);
            temp=roipoly();colorbar;
            droi=droi | temp;
            reply_end=input('more roi? Y/N\n','s');
        end
        
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

figure;
nn=ceil(sqrt(nte));
if roiflag ==0    
    for i= 1:nte subplot(nn,nn,i);hold on;display_fmri(dif_mean(:,:,i),mask);end
    % figure; display_fmri(anat,mask);
else
    for i=1:nte
        % fiugre; display_fmri(anat,mask(:,:,t));
        subplot(nn,nn,i);hold on;display_fmri(dif_mean(:,:,i),mask(:,:,i));
    end
end


% calculation for T2

meanum=nrep*2;t2num=nte;
%fitting for all four repeated measurement
b_all=zeros(meanum/2,nte);
c_all=zeros(meanum/2,nte);
for i=1:nte
    temp3=dif_all(:,:,i,:);
    temp3=reshape(temp3,rownum,columnnum,meanum/2);
    if roiflag ==0
        temp4=average_actpixel(1,mask,temp3);
        close;
        b_all(:,i)=temp4(:);
        temp3=con_all(:,:,i,:);
        temp3=reshape(temp3,rownum,columnnum,meanum/2);
        temp4=average_actpixel(1,mask,temp3);
        close;
        c_all(:,i)=temp4(:);
    else
        temp4=average_actpixel(1,mask(:,:,i),temp3);
        close; 
        b_all(:,i)=temp4(:);
    temp3=con_all(:,:,i,:);
    temp3=reshape(temp3,rownum,columnnum,meanum/2);
    temp4=average_actpixel(1,mask(:,:,i),temp3);
    close;
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
%  figure;;hold on;
%  for t=1:nte
%    plot(te,b_all(t,:),'*')
% end

%%% T2 fitting by all b points

te_all=te(1:nte); temp=te(1:nte);
for i=1:meanum/2-1
    te_all=[te_all;temp];
end
b_all_vector=reshape(b_all',meanum/2*nte,1);
t2num=nte; 
c_all_vector=reshape(c_all',meanum/2*nte,1);

for i=1:meanum/2
    for j=1:nte-1
        b_all_vector((i-1)*t2num+j)=exp(-(ti-te(nte))/bloodt1)/exp(-(ti-te(j))/bloodt1)*b_all_vector((i-1)*t2num+j);
        c_all_vector((i-1)*t2num+j)=exp(-(ti-te(nte))/bloodt1)/exp(-(ti-te(j))/bloodt1)*c_all_vector((i-1)*t2num+j);
    end
end

for tt=1:nrep
figure(1011);plot(te,log(b_all_vector(1+(tt-1)*nte : tt*nte)),'*');hold on;
figure(1021);plot(te,b_all_vector(1+(tt-1)*nte : tt*nte),'*');hold on;
end

[temp1,resid, jacob]=nlinfit(te_all/1000, b_all_vector,'monexp_model',[1000,10]);
conintval=nlparci(temp1,resid, jacob); %95% confidence interval for estimates
t2=1000./temp1(2)
aa=1000./conintval;
ci=aa(2,:)

figure(1011);hold on; cte=[1:180']';plot(cte,log(monexp_model(temp1,cte/1000)),'r-');hold off;
figure(1021);hold on;plot(cte,monexp_model(temp1,cte/1000),'r-');hold off;

[temp1_c,resid_c, jacob_c]=nlinfit(te_all/1000, c_all_vector,'monexp_model',[1000,10]);
conintval_c=nlparci(temp1_c,resid_c, jacob_c); %95% confidence interval for estimates
t2_c=1000./temp1_c(2);
aa_c=1000./conintval_c;
ci_c=aa_c(2,:);


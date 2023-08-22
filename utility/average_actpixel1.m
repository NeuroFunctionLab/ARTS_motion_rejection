function ave=average_actpixel(varargin)
%AVERAGE_ACTPIXEL display the averaged timecourse of the activated pixels.
%Input:
%	mode         --   only '1' is supported at this time
%	actpixel     --   N*N matrix, binary matrix, 1--activated;0--not activated
%	timecourse 1 --   N*N*#_of_times, time-course of the pixels
%	timecourse 2 --   ....
%	...
%	...               up to 4 timecourse arrays are allowed
%Output:
%	ave          --   #_of_times * arraynum; averaged timecourses of the actpixel
%Usage: e.g. average_actpixel(1,actpixel, img1);
%		 e.g. average_actpixel(1, actpixel, image1, image2, image3, image4);
%Note: the curve color was arranged such that image1--'blue', image2--'red', 
%image3--'green', image4--'black'.
%Hanzhang Lu
%Date: 01/03/01

%set variables
mode=varargin{1};
if (mode==1)
   actpixel = varargin{2};
   arraynum=nargin-2;
   temp=varargin{3};
   rownum=size(temp,1);
   columnnumber=size(temp, 2);
   timenum=size(temp,3);
   timecourses=zeros(rownum, columnnumber, timenum, arraynum);
   for i=1:arraynum
      temp=varargin{i+2};
      timecourses(:,:,:,i)=temp(:,:,:);
   end
end

%display
[row, column] = find(actpixel);
columnnum = ceil(sqrt(length(row)));
actnum=length(row);
ave=zeros(timenum, arraynum);

maximum=1;
minimum=1;
temp1=zeros(timenum, 1);
for i=1:arraynum
   temp=zeros(timenum, actnum);
   for j=1:actnum
      temp1(:)=timecourses(row(j), column(j), :, i);
      temp(:,j)=temp1(:);
   end
   dummy=mean(temp, 2);
   if i==1
      maximum=max(dummy);
      minimum=min(dummy);
   else
   	if (max(dummy)>maximum)
     		 maximum=max(dummy);
   	end
   	if (min(dummy)<minimum)
     		 minimum=min(dummy);
   	end
	end
   ave(:,i)=dummy(:);
end
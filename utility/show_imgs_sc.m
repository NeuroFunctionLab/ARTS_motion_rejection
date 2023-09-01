function output = show_imgs(varargin)
%SHOW_IMGS display the images stored in 'images'. The images
%to be displayed is indexed by 'index'. The output
%stores the data array corresponding to the displayed images
%Input:
%	varagin     --    The number of arguments is either 2 or 4. If 4,
%	                  they're 'images', 'index', 'low', 'high'. If 2,
%	                  only the first two arguments are used.
%	images     ---   image array to be displayed
%	index      ---   index of the images in 'images'
%	low        ---   the lower limit of the pixel values in the display
%	high       ---   the upper limit of the pixel values in the display
%Output:
%	output     ---   data of the images displayed
%Hanzhang Lu
%Date: 09/12/00
%Modified: 10/12/00

images=varargin{1};
index=varargin{2};
if (nargin==4)
   low=varargin{3};
   high=varargin{4};
end

%obtain the dimensionality of the images
row=size(images,1);
column=size(images,2);
numtotal=size(images,3);
numdisplay=max(size(index));
numrow=ceil(sqrt(numdisplay));
numcol = ceil(numdisplay/numrow);
%define temporary variables
temp=zeros(row,column);
output=zeros(row,column,numdisplay);

if nargin==2
   dummy=images(:,:, index);
   high=max(max(max(dummy)));
   low=min(min(min(dummy)));
end


figure;
%display the images specified in index
for i=1:numdisplay
   subplot(numrow,numcol,i);
   temp(:)=images(:,:,index(i));
   imagesc(temp, [low, high]);
   colorbar;  
%   colormap(jet(64));  % JET is a default which is equal to jet(64)
   colormap('gray');
   %hsv gray hot cool bone copper pink flag prism jet 
   %   colormap('gray');
   output (:,:,i)=temp(:,:);
   set(gca,'xtick',[]);
   set(gca,'ytick',[]);
   if numdisplay > 1
      title(sprintf('Img %d', i));
   end
end



function output = display_fmri (anat, actpixel)
%DISPLAY_FMRI combine the anatomical image with the activated
%pixels by setting a specific value for activated pixels.
%Then, a specifial colormap is used to display the image.
%Input:
%	anat      --  anatomical image to be modified; N*N matrix
%	actpixel  --  M*M matrix, '1' means active, '0' means inactive
%Output:
%	output     --  modified anatomical image
%Hanzhang Lu
%Date: 07/11/00

%normalize anatomical image and functional image
[rowa, columna]=size(anat);
rowf=size(actpixel, 1);
columnf=size(actpixel, 2);
scalerow=rowa/rowf;
scalecolumn=columna/columnf;

%modify anat
maximum = max(max(anat));
anat=anat/maximum;
maximum = max(max(anat));
%change the element values in the anatomical image according to activation map
for i=1:rowf
   for j=1:columnf
      if (actpixel(i, j)==1)
         anat(round((i-1)*scalerow+1):round(i*scalerow),round((j-1)*scalecolumn+1):round(j*scalecolumn)) = maximum*(1.0/60.0+1);
      end
   end
end

%display the modified image
figure;
%subplot(1,1,1);
imshow(anat, [0,maximum*(1.0/60.0+1)]);
%set new colormap so that the active pixel is pseudo-colored
c=gray(64);
c(64,2)=0;
c(64,3)=0;
colormap(c);

output = anat;
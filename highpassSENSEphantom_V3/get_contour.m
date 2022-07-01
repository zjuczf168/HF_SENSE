function [mask]=get_contour(I,Threshold)
% function [mask]=get_contour(I,Threshold)
%  I is the original image (3D or 2D). The salient feature 
%   should be positive.
%  Threshold is for filtering out noise (as 
%  normalized Gaussian). 
%  mask: a binary mask indicate the object/background region
% JIM Xiuquan JI, 2003, University of Illinois at Urbana-Champaign
% All Rights Reserved
%---------------------------------------------------
% Copyright (C) 2005 Jim Ji, Texas A&M University.
% All Rights Reserved.

% References:


% Version of 4-June-2005.

% Log:
% Created  2005 Jim Ji	Source code
% Updated
%-----------------------------------------------------

[Nx,Ny,total_slice]=size(I);

mask=zeros(Nx,Ny,total_slice);

for index=1:total_slice
a1=I(:,:,index);

%estimat the statistics of the background
	width=ceil(Nx/10);
	cornertl=a1(1:width, 1:width);
	cornertr=a1(1:width, Ny-width+1:Ny);
	cornerbl=a1(Nx-width+1:Nx, 1:width);
	cornerbr=a1(Nx-width+1:Nx, Ny-width+1:Ny);
	bkg=[cornertl(:); cornertr(:); cornerbl(:); cornerbr(:)];
% Extract the background
    maska = (a1-mean(bkg))/std(bkg)>Threshold;
    maska = imfill(maska,'holes');
    maska = bwselect(maska,round(Nx/2), round(Ny/2),8);
    
end

mask(:,:,index)=maska;

return





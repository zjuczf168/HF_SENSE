function [sos,coilimg]=recon_sumofsquares(kdata, IFT_off, R)
% This function takes multi coil k-space data and noise covariance matrix
% to provide a sum-of-quares image
% function [sos,coilimg]=recon_sumofsquares(kdata, IFT_off, R)
% Input:
%   kdata(Nfe,Npe,Ncoils): kspace phasedarray data
%   IFT_off: (optional), flag to turn off IFT if IFT_off==1. This means kdata represents images rather than Fourier data   
%   R: (optional) coil noise corvariance matrix, default R is identitity matrix (if
%    only kdata is provided)
%
% Output
%   sos: root sum of squares image
%    coilimg: the individual coil images
%   
%
% Jim Ji, Texas A&M University, 2005, All rights reserved
% jimji@tamu.edu

[Nx Ny Ncoil]=size(kdata);
sos = zeros(Nx,Ny);
coilimg=zeros(size(kdata));

IFT_flag = 1;
if nargin > 1
    if IFT_off==1
        IFT_flag = 0; % need IFT of the kdata
    end
end

if IFT_flag ==1
    for s=1:Ncoil
        coilimg(:,:,s) = ifftshift(ifft2(ifftshift(kdata(:,:,s))));
    end
else
    coilimg = kdata;
end


if nargin < 3
    sos = sqrt(sum([abs(coilimg)].^2,3));
elseif nargin == 3
    Rinv = inv(R);
    tmp = reshape(coilimg,Nx*Ny,Ncoil);
    sos = zeros(Nx*Ny,1);
    
    for index = 1:(Nx*Ny)
        vec = tmp(index,:); 
        sos(index) = vec * Rinv * vec';
    end    
    sos = sqrt(reshape(sos, Nx, Ny));
else
    display('Wrong number of input');
end

return

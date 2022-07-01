function [sensitivity_new]=sensitivity_hamming(roughsensitivity)
%function [sensitivity_new]=sensitivity_hamming(roughsensitivity);
% This function perform Hamming window filtering in the Fourier domain
% input: 
% rouoghsensitivity: raw,nisy sensitivity;
% output: 
% filtered sensitivity using Hamming window: 
% 
% Other low-pass windows such as Blackman, or Kaiser can be used by change 
% the the line "window_col=hamming(Npe);" to the desirable window
%
%----------------------------------------------
% Copyright (C) 2005 Jim Ji, Texas A&M University.
% All Rights Reserved.

% References:


% Version of 4-June-2005.

% Log:
% Created  2003 Jim Ji	Source code
% Updated
%-----------------------------------------------------

[Nfe,Npe,Ncoils]=size(roughsensitivity);

window_col=hamming(Npe);
hamm_window=repmat(window_col',[Nfe 1 Ncoils]);

sensitivity_inFourier = fftshift(fftshift(fft2(fftshift(fftshift(roughsensitivity,1),2)),2 ),1);
sensitivity_inFourier = sensitivity_inFourier.*hamm_window;

sensitivity_new = ifftshift(ifftshift(ifft2(ifftshift(ifftshift(sensitivity_inFourier,1),2)),2 ),1);

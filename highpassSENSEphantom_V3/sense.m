function [recon, sucessful_flag, gmap] = sense(channelsensitivity, reduced_kspace_data, FOV_reduction_factor)
%function [recon] = sense(channelsensitivity, reduced_kspace_data, FOV_reduction_factor)
%Input-
%   channelsensitivity: size (points_Freq Enc, lines_Phase Enc,
%       Number_channels), the spatial coil sensitivity function. It is
%       required that lines_Phase Enc must = Number_Phase Enc_Subsampled * FOV_reduction_factor; 
%       To enforce this condition, the channel sensitivity function may
%       need to be interpolated, i.e., by zero padding the k-space
%       calibration lines before inv-FT.
%   reduced_kspace_data: size (Number_Freq Enc, Number_Phase Enc_Subsampled,
%       Number_channels); As a convention, the center k-space line is in the 
%       middle. 
%   FOV_reduction_factor: It must be an integer, (FOV of sensitivity function)/(FOV of
%       reduced_kspace_data). If the data acquisition parameters results in
%       fractual reduction factor, then channelsensitivity functions must
%       be cropped and interpolated to make the reduction factor an integer 
%       before calling this function. 
%
%Output-
%   The sense reconstruction:  size (Number_Freq Enc,
%   Number_Phase Enc)
%   
%   sucessful_flag:  1 if reconstructin is successful and 0 if there is no
%   error in the process.
%
%   gmap: the g-factor map of the sense recon, size (Number_Freq Enc, Number_Phase Enc)
%
%--------------------------------------------------------
% The k-space subsampling is along the second dimension.; The leading dimension is
% the frequency encoding; The second dimension is phase encoding; The last
% dimension is the channel index;
%
% Jim Ji, 2002
% University of Illinois at Urbanan-Champaign
%----------------------------------------------------------
%----------------------------------------------
% Copyright (C) 2005 Jim Ji, Texas A&M University.
% All Rights Reserved.

% References:Sensitivity Encoding for fast MRI Scans (SENSE)
%                   K. Pruessmann, M. Weiger, P.Boesieger
%                   Magnetic Resonance in Medicine, 42:952-962 (1999)

% Version of 4-June-2005.

% Log:
% Created  2002 Jim Ji	Source code
% Updated  2005 Jim Ji  doc update, remove Receiver structure, add
%   FOV reduction_factor input
%-----------------------------------------------------


sucessful_flag = 0; % will return 1 if reconstructin is successful and there is no error in the process.

[Nfe,Npe_seg,Ncoils] = size(reduced_kspace_data);
Npe = size(channelsensitivity,2);
Rr = (Npe/Npe_seg); 

if FOV_reduction_factor - round(FOV_reduction_factor)~=0; 
    display('FOV_reduction_factor must be integer'); 
    return; 
end;

if Rr ~= FOV_reduction_factor; 
    display('Sensitivity function and k-space data does not agree with the FOV reduction factor!'); 
    display('You may need to interpolate the sensitivity functions. Type ''help sense'' to get more information.'); 
    return; 
end;

% generate the overlapped_img: (Nfe,Npe/R,Ncoils) 
% index as [-N/2:(N/2-1)].  Centered, unshifted 
overlapped_img = zeros(Nfe,Npe_seg,Ncoils);
for coil=1:Ncoils
 overlapped_img(:,:,coil)=ifftshift(ifft2(ifftshift(reduced_kspace_data(:,:,coil))));
end

% shift the overlapped image
% In the overlapped iamges above, the aliasing occurs as if (Npe/2+1) is
% the center of the unwrapped image, and the overlapped_img corresponds to
% an reduced FOV with the same center. Therefore, for R=2, the reduced
% FOV corresponds to (Npe/4 + 1): (Npe/2). But for R=3, the reduced FOV
% corresponds to (Npe/3+1);(Npe/3*2), exactly the same as (1:Npe/3). For
% the convenience of forming the SENSE equation, shift the first and sencod
% halves of the overlapped image if an even reduction factor is used.  
% For odd reduction factors, this shift is not necessary.
if mod(Rr,2)==0 %even R
    tmp_img  = overlapped_img;
	overlapped_img(:,1:(Npe_seg/2),:) = tmp_img(:,(Npe_seg/2+1):Npe_seg,:);
	overlapped_img(:,(Npe_seg/2+1):Npe_seg,:) = tmp_img(:,1:(Npe_seg/2),:);
end

% if mod(Rr,2)==0 %even R
%     tmp_img  = overlapped_img;
% 	overlapped_img(:,1:round(Npe_seg/2),:) = tmp_img(:,round(Npe_seg/2):Npe_seg,:);
% 	overlapped_img(:,round(Npe_seg/2):Npe_seg,:) = tmp_img(:,1:round(Npe_seg/2),:);
% end

%unwrap
recon=zeros(Nfe,Npe);
for pe=1:Npe_seg
   for i=1:Nfe
       if((pe+Npe_seg*(Rr-1))>Npe)
           delt=(pe+Npe_seg*(Rr-1))-Npe;
       else
           delt=0;
       end

      s1=squeeze(channelsensitivity(i,pe:Npe_seg:(pe-delt+Npe_seg*(Rr-1)),:));
      s=transpose(s1);
      I_mat=squeeze(overlapped_img(i,pe,:));
      recon(i,pe:Npe_seg:(pe-delt+Npe_seg*(Rr-1))) = transpose((inv(s'*s)*s'*I_mat));
      tmp = pinv(s'*s).*(s'*s); tmp = sqrt(abs(diag(tmp)));
      [l trash]=size(tmp);
      if (l ~= length(pe:Npe_seg:(pe+Npe_seg*(Rr-1))))
         tmp(l+1:length(pe:Npe_seg:(pe+Npe_seg*(Rr-1))),:)=0;
      end
      gmap(i,pe:Npe_seg:(pe+Npe_seg*(Rr-1)))= tmp';
  end 
end

sucessful_flag = 1;
return
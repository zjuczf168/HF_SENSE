close all; clear all; clc;
load('datao.mat');
% load('senseMap8ch.mat')
rawimage=datao(128:383,128:383,:);

for i=1:8
%     rawimage(:,:,i)=rawimage(:,:,i)+0.01*randi([-255,255],256,256);
    rawimage(:,:,i)=rawimage(:,:,i)+0.01*randi(max(max(max(round(abs(rawimage))))),256);

end

figure; imagesc(abs(rawimage(:,:,1))); colormap('gray');
figure; imagesc(abs(rawimage(:,:,4))); colormap('gray');
figure; imagesc(abs(rawimage(:,:,6))); colormap('gray');
figure; imagesc(abs(rawimage(:,:,8))); colormap('gray');

full_kspace_data = fftshift(fftshift(fft2(fftshift(fftshift(rawimage,1),2)),2 ),1);

[header.Nfe, header.Npe, header.num_coils] = size(full_kspace_data);
P=header.Nfe;
Q=header.Npe;
Rnoise=eye(header.num_coils);
FOV_reduction_factor=4;
c=24; w=8;
% c=18; w=6;


%%SENSE reconstruction
% k-space data subsampling
[reduced_kspace_data_SENSE, Subsampling_locations] = sample_kd(full_kspace_data,FOV_reduction_factor); 
coilimg=zeros(size(full_kspace_data));
 for s=1:header.num_coils
     coilimg(:,:,s) = ifftshift(ifft2(ifftshift(full_kspace_data(:,:,s))));
end
[recon,cmap]=adapt_array_2d(coilimg,Rnoise,0);  %Walsh MRM 2000
[recon_sense1, flag,gmap1] = sense(cmap, reduced_kspace_data_SENSE, FOV_reduction_factor); 
% [recon_sense1, flag,gmap1] = sense(senseMap, reduced_kspace_data_SENSE, FOV_reduction_factor); 

SOS = recon_sumofsquares(full_kspace_data,0, Rnoise);
error_image1 = abs(abs(recon_sense1) - SOS);
L2diff1=norm(error_image1(:));
% figure(1); imagesc(abs(recon_sense1)); colormap('gray'); axis square; colorbar;title('SENSE Recon with Walsh cmap'); 
% figure(2); imagesc(error_image1); colormap('gray'); axis square; colorbar; title(sprintf('L2 recon error: %s',L2diff1));
% figure(3); imagesc(abs(SOS)); colormap('gray'); axis square; colorbar; title('Sum of Squares Recon'); 
% figure(1); imagesc(abs(recon)); colormap('gray'); axis square; colorbar;title('Recon with Walsh cmap'); 
% figure(2); imagesc(abs(SOS)); colormap('gray'); axis square; colorbar; title('Sum of Squares Recon'); 
figure(3); imagesc(error_image1); colormap('gray'); axis square; colorbar; title(sprintf('L2 recon error: %s',L2diff1));
%% =====================================================================
% figure(1); imshow(abs(recon),[0,4.5e-6]); colormap('gray'); axis image; colorbar;title('Recon with Walsh cmap'); 
% figure(2); imshow(abs(SOS),[0,4.5e-6]); colormap('gray'); axis image; colorbar; title('Sum of Squares Recon'); 
% % figure(3); imagesc(error_image1); colormap('gray'); axis square; colorbar; title(sprintf('L2 recon error: %s',L2diff1));
% figure(3); imshow(error_image1,[0,2.2e-7]); colormap('gray'); axis image; colorbar; title(sprintf('L2 recon error: %s',L2diff1));
% 
figure(1); imshow(abs(recon_sense1),[0,240]); colormap('gray'); axis square; colorbar;title('Recon with Walsh cmap'); 
figure(2); imshow(abs(SOS),[0,240]); colormap('gray'); axis square; colorbar; title('Sum of Squares Recon'); 
% figure(3); imagesc(error_image1); colormap('gray'); axis square; colorbar; title(sprintf('L2 recon error: %s',L2diff1));
figure(3); imshow(error_image1,[0,173]); colormap('gray'); axis square; colorbar; title(sprintf('L2 recon error: %s',L2diff1));

%=================================
% hpSENSE
% ================================


reduced_kspace_data_hpSENSE=zeros(header.Nfe, header.Npe, header.num_coils); % 
highpassed_kspace_data=zeros(header.Nfe, header.Npe, header.num_coils); %
highpassed_reduced_kspace_data=zeros(header.Nfe, header.Npe/FOV_reduction_factor, header.num_coils); %
highpassed_kspace_data1=zeros(header.Nfe, header.Npe, header.num_coils);
highpassed_reduced_kspace_data2=zeros(header.Nfe, header.Npe/FOV_reduction_factor, header.num_coils); 
coilimg1=zeros(size(full_kspace_data));

H_1 = zeros(P,Q); %high pass filter, the same as hp-GRAPPA
H_1_inverse = zeros(P,Q);
for xx = (-P/2):1:(P/2)-1
    for yy = -Q/2:1:Q/2-1  
         H_1(xx+(P/2)+1,yy+Q/2+1) =1-(1+exp((sqrt(xx^2+yy^2)-c)/w))^-1+(1+exp((sqrt(xx^2+yy^2)+c)/w))^-1;
         H_1_inverse(xx+(P/2)+1,yy+Q/2+1) =1/(1-(1+exp((sqrt(xx^2+yy^2)-c)/w))^-1+(1+exp((sqrt(xx^2+yy^2)+c)/w))^-1);
    end
end

for coil = 1:header.num_coils
%     reduced_kspace_data_hpSENSE(:,Subsampling_locations,coil) = full_kspace_data(:,Subsampling_locations,coil);
%     highpassed_kspace_data(:, :, coil)=H_1.* reduced_kspace_data_hpSENSE(:, :, coil);
%     highpassed_reduced_kspace_data(:,:,coil) =  highpassed_kspace_data(:,Subsampling_locations,coil); 
    
    highpassed_kspace_data1(:, :, coil)=H_1.* full_kspace_data(:, :, coil);
%     coilimg1(:,:,coil) = ifftshift(ifft2(ifftshift(full_kspace_data(:,:,coil))));
    coilimg1(:,:,coil) = ifftshift(ifft2(ifftshift(highpassed_kspace_data1(:,:,coil))));
%     coilimg1(:,:,coil) = ifftshift(ifft2(ifftshift(highpassed_kspace_data(:,:,coil))));
    highpassed_reduced_kspace_data2(:,:,coil)=highpassed_kspace_data1(:,Subsampling_locations,coil); 

end

[recon1,cmap1]=adapt_array_2d(coilimg1,Rnoise,0);
% [recon_hpSENSE, flag,gmap1] = sense(cmap1, highpassed_reduced_kspace_data, FOV_reduction_factor); 
[recon_hpSENSE, flag,gmap1] = sense(cmap1, highpassed_reduced_kspace_data2, FOV_reduction_factor); 

recon_full_kspace_data = fftshift(fftshift(fft2(fftshift(fftshift(recon_hpSENSE,1),2)),2 ),1);
inverse_filtered_kdata=H_1_inverse.*recon_full_kspace_data;
recon_HIGHPASSSENSE = ifftshift(ifftshift(ifft2(ifftshift(ifftshift(inverse_filtered_kdata,1),2)),2 ),1);

% recon_full_kspace_data = fftshift(fftshift(fft2(fftshift(fftshift(recon1,1),2)),2 ),1);
% inverse_filtered_kdata=H_1_inverse.*recon_full_kspace_data;
% recon_HIGHPASSSENSE = ifftshift(ifftshift(ifft2(ifftshift(ifftshift(inverse_filtered_kdata,1),2)),2 ),1);

error_image_HIGHPASSSENSE=abs(abs(recon_HIGHPASSSENSE)-SOS);
L2diff_HIGHPASSSENSE=norm(error_image_HIGHPASSSENSE(:));

% error_image_HIGHPASSSENSEsubSENSE=abs(abs(recon_HIGHPASSSENSE)-recon_sense1);
% L2diff_error_image_HIGHPASSSENSEsubSENSE=norm(error_image_HIGHPASSSENSEsubSENSE(:));

% figure(4); imagesc(abs(recon_HIGHPASSSENSE)); colormap('gray'); axis square; colorbar;title('hpSENSE Recon'); 
figure(5); imagesc(error_image_HIGHPASSSENSE); colormap('gray'); axis square; colorbar; title(sprintf('high-pass SENSE, L2 recon error: %s',L2diff_HIGHPASSSENSE));
% figure(6); imagesc(abs(recon_hpSENSE)); colormap('gray'); axis square; colorbar;title('high pass filtered SENSE'); 
% figure(6); imagesc(abs(recon1)); colormap('gray'); axis image; colorbar;title('high pass filtered SENSE'); 

% figure(6); imagesc(abs(recon_hpSENSE)); colormap('gray'); axis square; colorbar;title('high pass filtered SENSE'); 
% figure(7); imagesc(error_image_HIGHPASSSENSEsubSENSE); colormap('gray'); axis square; colorbar; title(sprintf('HIGHPASSSENSEsubSENSE, L2 recon error: %s',L2diff_error_image_HIGHPASSSENSEsubSENSE));
% figure(9); imagesc(abs(coilimg(:,:,7))); colormap('gray'); axis image; colorbar;title('coil image'); 
figure(4); imshow(abs(recon_HIGHPASSSENSE),[0, 240]); colormap('gray'); axis square; colorbar;title('hpSENSE Recon'); 
figure(5); imshow(error_image_HIGHPASSSENSE,[0,173]); colormap('gray'); axis square; colorbar; title(sprintf('hpSENSE, L2 recon error: %s',L2diff_HIGHPASSSENSE));
% figure(6); imagesc(abs(recon_hpSENSE)); colormap('gray'); axis square; colorbar;title('high pass filtered SENSE'); 
figure(6); imshow(abs(recon_hpSENSE),[0, 240]); colormap('gray'); axis square; colorbar;title('high pass filtered SENSE'); 

NRMSE_SENSE=sqrt(sum(sum((SOS-abs(recon_sense1)).^2))/sum(sum(SOS.^2)));
NRMSE_HFSENSE=sqrt(sum(sum((SOS-abs(recon_HIGHPASSSENSE)).^2))/sum(sum(SOS.^2)));

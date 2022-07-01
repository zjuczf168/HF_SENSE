% clear all
% clc
% close all
% startup

%Step1, simulating the subsampled data and coil-calibration data
% =====================================================================
% 20170827 zjc revised
% apply the high pass filter after subsampling full kspace data
% then obtain the highpassed_reduced_kspace_data for SENSE reconstruction

[header.Nfe, header.Npe, header.num_coils] = size(full_kspace_data);
FOV_full = 1;
FOV_subsampled = 1/4;
FOV_reduction_factor = round(FOV_full/FOV_subsampled); 
% Rnoise = eye(8); %noise correlation matrix betweent channels
Rnoise = eye(16); %noise correlation matrix betweent channels
%----------------------------------------------
% the sampling location is the subsampling location; need to take care of
% the even/odd reduction factor; make sure that the central k-space
% (C/2+1) is sampled and the subsampled phase-encodings is an even number.
%----------------------------------------------
tmp1 = (header.Npe/2+1-FOV_reduction_factor):(-FOV_reduction_factor):1;
tmp2 = (header.Npe/2+1):FOV_reduction_factor:header.Npe;
Subsampling_locations = [flipud(tmp1(:)); tmp2(:)]; 
if length(tmp1)>length(tmp2)
    Subsampling_locations=Subsampling_locations(2:length(Subsampling_locations));
    else if length(tmp1)<length(tmp2)
        Subsampling_locations=Subsampling_locations(1:length(Subsampling_locations)-1);
    end
end

reduced_kspace_data=zeros(header.Nfe, header.Npe, header.num_coils);
for coil = 1:header.num_coils
    reduced_kspace_data(:,Subsampling_locations,coil) = full_kspace_data(:,Subsampling_locations,coil);
end

P=header.Nfe;
Q=header.Npe;
c=14; w=10;
% c=14; w=12;
% c=14; w=2;
H_1 = zeros(P,Q);
H_1_inverse = zeros(P,Q);
highpassed_kspace_data=zeros(header.Nfe, header.Npe, header.num_coils);
highpassed_reduced_kspace_data=zeros(header.Nfe, header.Npe/FOV_reduction_factor, header.num_coils);
% -------------------------------------------------------
% test the effect of high pass filter

for xx = (-P/2):1:(P/2)-1
    for yy = -Q/2:1:Q/2-1  
H_1(xx+(P/2)+1,yy+Q/2+1) =1-(1+exp((sqrt(xx^2+yy^2)-c)/w))^-1+(1+exp((sqrt(xx^2+yy^2)+c)/w))^-1;
H_1_inverse(xx+(P/2)+1,yy+Q/2+1) =1/(1-(1+exp((sqrt(xx^2+yy^2)-c)/w))^-1+(1+exp((sqrt(xx^2+yy^2)+c)/w))^-1);
    end
end

for coil=1:header.num_coils
    highpassed_kspace_data(:, :, coil)=H_1.* reduced_kspace_data(:, :, coil);
    highpassed_reduced_kspace_data(:,:,coil) =  highpassed_kspace_data(:,Subsampling_locations,coil); 
end

% =========================================================================


% k-space data subsampling
% [reduced_kspace_data, Subsampling_locations] = sample_kd(highpassed_kspace_data,FOV_reduction_factor);  %high-pass SENSE


[reduced_kspace_data_sense, Subsampling_locations] = sample_kd(full_kspace_data,FOV_reduction_factor); %SENSE

% save('subsampled_kdata','reduced_kspace_data', 'Subsampling_locations','FOV_reduction_factor','Rnoise');

% coil sensitivity calibration data from the central k-space
Num_centrallines = 32;
[Nfe,trash,Ncoil] = size(reduced_kspace_data_sense);

for ch=1:Ncoil
    central_kdata(:,:,ch) = datacrop2d(full_kspace_data(:,:,ch),Nfe,Num_centrallines);
end
save('16ch_centralkdata','central_kdata');


% clear 
%Recon processing
% load('subsampled_kdata');

%sensitivity estimation
Npe_tobe = size(reduced_kspace_data_sense,2)*FOV_reduction_factor;

% sensitivity estimate

%filename='coilimages.mat'; datatype =1;
filename='16ch_centralkdata.mat'; datatype =2;


% Filters.type ='nofilter'; 
% Filters.type ='polynomial'; Filters.Poly_Order = 1; Filters.poly_Wc=5; Filters.poly_Wr=5; %for poly 
% Filters.type ='walsh';   Filters.walsh_width = 3; % for walsh method, averaging window size
% Filters.type ='median';   Filters.median_width = 3; % median filter window size
% Filters.type ='wavelet';  Filters.wavelet_name = 'bior5.5'; Filters.wavelet_levels=3; % for wavelet
 Filters.type ='hamming'; %
 [channelsensitivity, ROI_mask]= sensitivity_estimation(filename, datatype, Npe_tobe,Filters,Rnoise); %sense sensitivity

%show_group(abs(channelsensitivity),3,4); colorbar; title('filtered map: magnitude'); colorbar
%show_group(angle(channelsensitivity),4,4); colorbar; title('filtered map: phase'); colorbar

%sense reconstruction
[recon_sense, flag,gmap] = sense(channelsensitivity, reduced_kspace_data_sense, FOV_reduction_factor); %SENSE reconstruction
% [recon_sense_highpassfiltered, flag,gmap] = sense(channelsensitivity, reduced_kspace_data, FOV_reduction_factor); %high-pass SENSE
[recon_sense_highpassfiltered, flag,gmap] = sense(channelsensitivity, highpassed_reduced_kspace_data, FOV_reduction_factor); %high-pass SENSE

%error map

%sum of square recon as "standard"
currentdir = pwd;
% load(sprintf('%s\\..\\..\\data\\data_brain',currentdir))
% load(sprintf('%s\\data_brain',currentdir))
load('16ch_full_kspace_data');

sos = recon_sumofsquares(full_kspace_data,0, Rnoise);

error_image = abs(abs(recon_sense) - sos);

L2diff=norm(error_image(:));
% ---------------------------------------------------------
% 20170826 zjc added
% inverse 2DFT the recon_sense to obtain the recon_full_kspace_data
% inverse high pass filtering the reconstructed kspace data
% then inverse 2DFT to obtain the final high-pass SENSE image
% 
recon_full_kspace_data=fftshift(fftshift(fft2(fftshift(fftshift(recon_sense_highpassfiltered,1),2)),2 ),1);
inverse_filtered_kdata=H_1_inverse.*recon_full_kspace_data;
recon_HIGHPASSSENSE= ifftshift(ifftshift(ifft2(ifftshift(ifftshift(inverse_filtered_kdata,1),2)),2 ),1);

error_image_HIGHPASSSENSE=abs(abs(recon_HIGHPASSSENSE)-sos);
L2diff_HIGHPASSSENSE=norm(error_image_HIGHPASSSENSE(:));
% --------------------------------------------------------------

figure(1); imagesc(abs(recon_sense)); colormap('gray'); axis square; colorbar;title('sense recon'); 
figure(2); imagesc(abs(sos)); colormap('gray'); axis square; colorbar;title('sum of squares recon'); 
figure(3); imagesc(error_image); colormap('gray'); axis square; colorbar; title(sprintf('L2 recon error: %s',L2diff));
figure(4); imagesc(gmap); colorbar; axis square; title('g-factor map');
figure(5); imagesc(abs(recon_HIGHPASSSENSE)); colormap('gray'); axis square; colorbar;title('high-pass SENSE recon'); 
figure(6); imagesc(error_image_HIGHPASSSENSE); colormap('gray'); axis square; colorbar; title(sprintf('high-pass SENSE, L2 recon error: %s',L2diff_HIGHPASSSENSE));


% figure(1); imagesc(abs(recon_sense)); colormap('gray');  colorbar;title('sense recon'); 
% figure(2); imagesc(abs(sos)); colormap('gray');  colorbar;title('sum of squares recon'); 
% figure(3); imagesc(error_image); colormap('gray'); colorbar; title(sprintf('L2 recon error: %s',L2diff));
% figure(4); imagesc(gmap); colorbar;  title('g-factor map');

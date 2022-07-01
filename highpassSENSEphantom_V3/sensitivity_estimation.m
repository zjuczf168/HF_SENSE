function   [channelsensitivity, ROI_mask]= sensitivity_estimation(filename, datatype, Npe_tobe,Filters,Rnoise)
%function   [channelsensitivity, ROI_mask]= channelsensitivityimation(filename, datatype, Npe_tobe,Filters,Rnoise);
% Input
% filename: the file which contains the coil calibration data
%           This should be a .mat
%           file with only one variable with the same name as the filename. The
%           dimension is (Frequency_enc, Phase_enc, coil)
% datatype: 
%   1) file contains central k-space data; the data should be
%   collecged with a FOV large enough to avoid aliasing; It should comes
%   from the central portion of the k-space symmetrically from
%   -L/2:1:(L/2-1) assuming total L lines.
%   2) file contains low-resolution coil images, the FOV must be equal to
%   reduced FOV times FOV_reduction_factor; but it can have arbitrary
%   number of lines --- the program will resample the image so it will have
%   the correct # of lines (i.e., Npe_tobe = Npe_of_the_reduced_k-space_data *
%   FOV_reduction_factor)
%   3) file contains the sensitivity functions alreay estimated by measured elsewhere. 
%
% Npe_tobe: the desirable lines along Phase Encoding direction of the
% sensitivity function.
%
% For correct unaliasing operation, the size the
% sensitivity functions must be correctly set. For example, if the reduced
% k-space data were acquired with reduced_FOV, with Npe_of_the_reduced_k-space_data total lines, then
% the FOV of the sensitivity calibration data must be
% reduced_FOV*FOV_reduction_factor; correspondingly, the total lines of the
% sensitivity function must be Npe_tobe = Npe_of_the_reduced_k-space_data * FOV_reduction_factor.
% For this purpose, it needs to know the number of lines in
% the subsampled k-space data; the FOV of the coil calibration images/data,
% and the reduced FOV of the subsampled data
%
% Filters: the control for post filtering
%    Filters.type:
%       'nofilter': no filter applied
%       'polynomial': filter the rough sensitivity using
%                   polynomial fitting as discussed in the SENSE method
%       'median': apply medial filter
%       'wavelet': apply wavelet filter
%       'hamming': applyy hamming filter
%       'walsh': using singular value decomposition to obtain sensitivity;
%   Filters.Poly_Order = 1; Filters.poly_Wc=5; Filters.poly_Wr=5; %for poly
%   Filters.walsh_width = 3; % for walsh method, averaging window size
%   Filters.median_width = 3; % median filter window size
%   Filters.wavelet_name = 'bior5.5'; filters.wavelet_levels=3; % for wavelet

% Rnoise: the noise corvariance matrix of the phased-array coils.
% 
% OUTPUT:
%    channelsensitivity: estimated sensitivity function of all coils with the
%        desirable size and FOV
%    ROI_mask: a binary mask indicating region of the interest, i.e., valid
%         object region. can be rough estimate
%
%
% Copyright (C) 2005 Jim Ji, Texas A&M University.
%  All Rights Reserved.

 if exist(filename)==0
    display(sprintf('Coil sensitivity calibration file %s does not exist',filename));
    return
end
     
% input or produce the coil images
switch datatype
    case 1  %read in the coil calibration images functions
        coilimages = double(importdata(filename));
        [Nfe,Npe,Ncoil] = size(coilimages);
        
        %if the coil image is too big or small, interpolate/decimate it so that the FOV matches the reduced_FOV*FOV_reduction_factor; 
        if Npe ~= Npe_tobe 
            tmpcoilimages = zeros(Nfe,Npe_tobe,Ncoil);  
            for ch =1:Ncoil
                tmp =  fftshift(fft2(fftshift(coilimages(:,:,ch))));    
                if Npe < Npe_tobe 
                    tmp = zerofill2d(tmp,Nfe,Npe_tobe);     
                else % image bigger than Npe_tobe
                    tmp = datacrop2d(tmp, Nfe, Npe_tobe);
                end
                tmpcoilimages(:,:,ch) =  ifftshift(ifft2(ifftshift(tmp)));    
            end
            coilimages=tmpcoilimages;
        end
        clear tmpcoilimages
        
    case 2 % the central k-space calibration data
        calibration_kspacedata = double(importdata(filename));
        [Nfe,Npe,Ncoil] = size(calibration_kspacedata);
            
        if (Npe ~= Npe_tobe)
            tmpdata = zeros(Nfe,Npe_tobe,Ncoil);  
            for ch =1:Ncoil
                if (Npe < Npe_tobe)
                    tmpdata(:,:,ch) = zerofill2d(calibration_kspacedata(:,:,ch),Nfe,Npe_tobe);
                else
                    tmpdata(:,:,ch) = datacrop2d(calibration_kspacedata(:,:,ch),Nfe,Npe_tobe);
                end
            end
            calibration_kspacedata = tmpdata;
        end
        
        for ch =1:Ncoil
            coilimages(:,:,ch) =  ifftshift(ifft2(ifftshift(calibration_kspacedata(:,:,ch))));
        end
        clear tmpdata
        
    case 3 %matlab file contains the sensitivity file roughsensitivity(Nfe,Npe,Ncoil)
        roughsensitivity =double(importdata(filename));
        [Nfe,Npe,Ncoil] = size(roughsensitivity);
        
        %if the coil image is too big or small, interpolate/decimate it so that the FOV matches the reduced_FOV*FOV_reduction_factor; 
        if Npe ~= Npe_tobe 
            coilimages = zeros(Nfe,Npe_tobe,Ncoil);  
            for ch =1:Ncoil
                tmp =  fftshift(fft2(fftshift(roughsensitivity(:,:,ch))));    
                if Npe < Npe_tobe 
                    tmp = zerofill2d(tmp,Nfe,Npe_tobe);     
                else % image bigger than Npe_tobe
                    tmp = datacrop2d(tmp, Nfe, Npe_tobe);
                end
                coilimages(:,:,ch) =  ifftshift(ifft2(ifftshift(tmp)));    
            end
            roughsensitivity=coilimages;
        end
        
    otherwise
        display('Data type for coil calibration must be 1, 2, or 3');
        
end


    coilimages(find(abs(coilimages)==0)) = 1e-18; %1e-15 is a small regularization
    sum_of_squares = recon_sumofsquares(coilimages,1, Rnoise); 
    ROI_mask = get_contour(sum_of_squares,2);

if datatype ==1 | datatype ==2
    %--------sum of squares image--------------------   
    for ch=1:Ncoil
        roughsensitivity(:,:,ch)=coilimages(:,:,ch)./(sum_of_squares); 
    end
end


%filtering options
switch Filters.type
    case 'nofilter'
        channelsensitivity = roughsensitivity;
        
    case 'polynomial'
        [channelsensitivity]=sensitivity_polyfilt(roughsensitivity,Filters.Poly_Order,Filters.poly_Wc,Filters.poly_Wr);

    case 'walsh'
        if datatype ==3
            display('For Walsh method, datatype must be 1 or 2, i.e., input file must contain either coil calibration images or k-space data');
            return
        end
        channelsensitivity = sensitivity_walsh(coilimages, Filters.walsh_width);
        
    case 'median'
        [channelsensitivity]=sensitivity_median(roughsensitivity,Filters.median_width);
        
    case 'wavelet'
        [channelsensitivity]=sensitivity_wavelet(roughsensitivity, Filters.wavelet_name,Filters.wavelet_levels);
        
    case 'hamming'
        [channelsensitivity]=sensitivity_hamming(roughsensitivity);
        
    otherwise
                channelsensitivity = roughsensitivity;
end



%show_group(abs(roughsensitivity),1,4); colorbar; title('rough map');
%show_group(angle(roughsensitivity),2,4); colorbar; title('rough map');
%show_group(abs(channelsensitivity),3,4); colorbar; title('filtered map');
%show_group(angle(channelsensitivity),4,4); colorbar; title('filtered map');

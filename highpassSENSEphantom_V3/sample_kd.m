function [reduced_kspace_data, Subsampling_locations] = sample_kd(full_kspace_data,FOV_reduction_factor)
%function [reduced_kspace_data, Subsampling_locations] = sample_kd(full_kspace_data,FOV_reduction_factor)
% simulate subsampled k-space data from each coil; reduction in the second
% direction, i.e., phase encoding direction;
%----------------------------------------------
%input:
% full kspace data : the kspace data without subsampling, i.e., without
% aliasing
% FOV_reduction_factor:  data subsampling_factor, must be integer
%-------------------------------------------
% output: 
% reduced/subsampled kspace data: the data has an apparent reduced FOV.  
% Subsampling locations: the corresponding line numbers along the phase encoding that were kept from the full kspace data;
%-------------------------------------------
% Copyright (C) 2005 Jim Ji, Texas A&M University.
% All Rights Reserved.
%
% References:
%
% Version of 31-July-2005.

% Log:
% Modified 
% 2004 from kspacedata_sim.m Swati Rane
% 2005/7/31: JJ cleanup the structure Sampling, change to generic parameters. Edit function header.
% Updated
%-----------------------------------------------------

[Ro,C, num_coils]=size(full_kspace_data);
% [header.Nfe, header.Npe, header.num_coils] = size(full_kspace_data);

%----------------------------------------------
% the sampling location is the subsampling location; need to take care of
% the even/odd reduction factor; make sure that the central k-space
% (C/2+1) is sampled and the subsampled phase-encodings is an even number.
%----------------------------------------------
% 20171128revised
% zjc
tmp1 = (C/2+1-FOV_reduction_factor):(-FOV_reduction_factor):1;
tmp2 = (C/2+1):FOV_reduction_factor:C;
Subsampling_locations = [flipud(tmp1(:)); tmp2(:)]; 
% if length(tmp1)>length(tmp2)
%     Subsampling_locations=Subsampling_locations(2:length(Subsampling_locations));
%     else if length(tmp1)<length(tmp2)
%         Subsampling_locations=Subsampling_locations(1:length(Subsampling_locations)-1);
%     end
% end
% Subsampling_locations=1:FOV_reduction_factor:header.Npe;
for coil = 1:num_coils
%   for coil = 1:header.num_coils

    reduced_kspace_data(:,:,coil) = full_kspace_data(:,Subsampling_locations,coil);
end

return
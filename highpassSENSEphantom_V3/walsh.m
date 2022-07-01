function walsh(MR)
%Routine to use Walsh et al. his method to compute the coil sensitivity
% maps for 2D/3D scans. This is by far the fastest way to generate the coil
% maps. 
%
% 20170717 - T.Bruijnen

%% Logic
if ~strcmpi(MR.UMCParameters.AdjointReconstruction.CoilSensitivityMaps,'walsh') 
    return;end

%% walsh
% Check whether its multi 2D or 3D data
if (strcmpi(MR.Parameter.Scan.ScanMode,'2D')) || (strcmpi(MR.UMCParameters.AdjointReconstruction.NufftType,'2D') && strcmpi(MR.Parameter.Scan.AcqMode,'Radial'))
    
    % Track progress
    parfor_progress(MR.UMCParameters.AdjointReconstruction.KspaceSize{1}(3));

    % (m)2D cases
    for z=1:MR.UMCParameters.AdjointReconstruction.KspaceSize{1}(3)
        
        % Estimate csm
        [~,MR.Parameter.Recon.Sensitivities(:,:,:,z)]=openadapt(permute(MR.Data{1}(:,:,z,:),[4 1 2 3]));
        
        % Track progress
        parfor_progress;
        
    end
    
    % Reset parfor
    parfor_progress(0);
    
elseif strcmpi(MR.Parameter.Scan.ScanMode,'3D')
    
    % 3D case
    [~,MR.Parameter.Recon.Sensitivities]=openadapt(permute(MR.Data{1},[4 1 2 3]));
        
end

% Reshape to [x,y,z,coil] dimension
MR.Parameter.Recon.Sensitivities=permute(MR.Parameter.Recon.Sensitivities,[2 3 4 1]);

% END
end
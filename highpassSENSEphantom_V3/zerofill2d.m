function [outdata]=zerofill2d(data,Nx_tobe,Ny_tobe)
%function [outdata]=zerofill2d(data,Nx_tobe,Ny_tobe)
% Pad zeros to the matrix data so it have the overall size of [Nx_tobe, Ny_tobe)
% Assume data is centered at zero frequency, -Nx_orig/2: Nx_orig/2-1,
%     -Ny_orig/2: Ny_orig/2-1 
% fill the other part with zeros
%----------------------------------------------
%
% Copyright (C) 2005 Jim Ji, Texas A&M University.
% All Rights Reserved.

[Nx_orig,Ny_orig]=size(data);

if Nx_orig>Nx_tobe  | Ny_orig>Ny_tobe
   'Wrong dimension parameter Nx_tobe or Ny_tobe';
   return
end

outdata=zeros(Nx_tobe,Ny_tobe);
outdata( (-Nx_orig/2+Nx_tobe/2+1):(Nx_orig/2+Nx_tobe/2),...
         (-Ny_orig/2+Ny_tobe/2+1):(Ny_orig/2+Ny_tobe/2) ) =data;
return
      



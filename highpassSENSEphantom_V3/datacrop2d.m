function [outdata] = datacrop2d(data, Wr_tobe, Wc_tobe)
% This function is to crop the central portion of a matrix data 
% Input: 
%   data, 2D matrix
%   Wc_tobe, Wr_tobe: the desirable columns and rows after cropping.
% Output:
%   outdata: Wr_tobe * Wc_tobe matrix
%
%  It assumes the center of the matrix is (Nc/2+1, Nr/2+1) for even lines,
%  and ((Nc+1)/2, (Nr+1)/2) for odd lines
%----------------------------------------------
%
% Copyright (C) 2005 Jim Ji, Texas A&M University.
% All Rights Reserved.
% -----------------------------------------------------
% 20170825 zjc revised
% refer to the PULSAR user manual for SENSE, 
% the output should be Wr_tobe * Wc_tobe*Ncoil,
% thus revise the function input to "data, Wr_tobe, Wc_tobe",
% and revise Nc and Nr for the 2nd and 1st dimension
% of input data respectively
% this has no effect on the entire code, just for clarity
% -----------------------------------------------------
outdata = zeros(Wr_tobe,Wc_tobe);

%center of the data
Nc = size(data,2);
if mod(Nc,2) == 0; 
    center_col = Nc/2+1;
else
    center_col = (Nc+1)/2;
end

Nr = size(data,1);
if mod(Nr,2) == 0; 
    center_row = Nr/2+1;
else
    center_row = (Nr+1)/2;
end

% start index
if mod(Wc_tobe,2) == 0; %even 
    start_col = -Wc_tobe/2;
else
    start_col = -(Wc_tobe-1)/2;
end

if mod(Wr_tobe,2) == 0; %even 
    start_row = -Wr_tobe/2;
else
    start_row = -(Wr_tobe-1)/2;
end

outdata = data(center_row + (start_row:(start_row+Wr_tobe-1)), center_col + (start_col:(start_col+Wc_tobe-1)));

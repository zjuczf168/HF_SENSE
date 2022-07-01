function [outdata] = datacrop2d(data, Wc_tobe, Wr_tobe)
% This function is to crop the central portion of a matrix data 
% function [outdata] = datacrop2d(data, Wc,Wr)
% Input: 
%   data, 2D matrix
%   Wc_tobe, Wr_tobe: the desirable columns and rows after cropping.
% Output:
%   outdata: Wc_tobe * Wr_tobe matrix
%
%  It assumes the center of the matrix is (Nc/2+1, Nr/2+1) for even lines,
%  and ((Nc+1)/2, (Nr+1)/2) for odd lines
%----------------------------------------------
%
% Copyright (C) 2005 Jim Ji, Texas A&M University.
% All Rights Reserved.

outdata = zeros(Wc_tobe,Wr_tobe);

%center of the data
Nc = size(data,1);
if mod(Nc,2) == 0; 
    center_col = Nc/2+1;
else
    center_col = (Nc+1)/2;
end

Nr = size(data,2);
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

outdata = data(center_col + (start_col:(start_col+Wc_tobe-1)), center_row + (start_row:(start_row+Wr_tobe-1)));

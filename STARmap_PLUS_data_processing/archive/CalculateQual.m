function [quality_matrix] = CalculateQual(cc_img, Nchannel)
%CalculateQual
%   test

% Calculate comp
sq = cc_img .^ 2;
comp_denominator = sum(sq, 4) .^ 0.5; % dim4 - channel
comp_denominator = repmat(comp_denominator, [1 1 1 Nchannel]);
comp_matrix = cc_img ./ comp_denominator; %4D

% Get max comp and barcode 
[max_matrix, barcode_matrix] = max(comp_matrix, [], 4);

% Get mask
max_matrix = eq(comp_matrix, repmat(max_matrix, [1 1 1 Nchannel]));
mask_matrix = sum(double(max_matrix), 4) ~= 1;

% Change pixel's barcode to 0 if pixel has multiple max
barcode_matrix(mask_matrix) = 0;

% Calculate quality score 
quality_sub = sum((comp_matrix - double(max_matrix)).^2, 4) .^ 0.5;
quality_matrix = ones(size(barcode_matrix)) - quality_sub;
quality_matrix(mask_matrix) = 0; % Change pixel's quality score to 0 if pixel has multiple max

 % Change pixel's comp to 0 if pixel has multiple max
comp_matrix(repmat(mask_matrix, [1 1 1 Nchannel])) = 0;


end


function [output_matrix] = im_cell2mat( input_cell )
%UNTITLED This functions is used to convert a cell containing multiple 3d
%images to a corresponding matrix 
%   -----IO-----
%   input_cell: a cell contains multiple image stacks 
%   output_matrix: a multi-dimensional matrix with all input images 

    N = numel(input_cell);
    
    output_matrix = uint8(zeros([size(input_cell{1}), N]));
    for i = 1:N
        output_matrix(:,:,:,:,i) = input_cell{i};
    end


end


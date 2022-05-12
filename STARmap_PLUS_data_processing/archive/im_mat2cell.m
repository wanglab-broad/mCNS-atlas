function [output_cell] = im_mat2cell( input_matrix )
%UNTITLED This functions is used to convert a 5D matrix containing multiple 3d
%images to a corresponding cell 
%   -----IO-----
%   input_matrix: a 5D matrix contains multiple image stacks 
%   output_cell: a cell with all input images 
    
    N = size(input_matrix, 5);
    output_cell = {N};

    for i = 1:N
        output_cell{i} = input_matrix(:,:,:,:,i);
    end


end


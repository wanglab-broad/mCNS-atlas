function [thres_stack, bw_stack] = lowerboundarythresh_each_channel(input_stack, lower_bound_factor)
%lowerboundarythresh_each_channel is used to apply global lower boundary thresholding to all
%channels, and output both a thresholded version and binary mask of the
%input stack
%   -----IO-----
%   input_stack: 5D image stack
%   lower_bound_factor: fraction of the max intensity used as lower
%   boundary
%   thres_stack: Threshold applied input_stack
%   bw_stack: Binary mask of the thresholding

    Nround = size(input_stack, 5);
    Nchannel = size(input_stack, 4);
    
    thres_stack = uint8(zeros(size(input_stack)));
    bw_stack = zeros(size(input_stack));
    
    for r = 1:Nround
        curr_round = input_stack(:,:,:,:,r);
        for c = 1:Nchannel
            curr_channel = curr_round(:,:,:,c);
            curr_max = max(curr_channel, [], 'all');
            curr_boundary = curr_max * lower_bound_factor;
            curr_bw = curr_channel > curr_boundary;
            curr_thres = curr_channel;
            curr_thres(~curr_bw) = 0;
            
            bw_stack(:,:,:,c,r) = curr_bw;
            thres_stack(:,:,:,c,r) = uint8(curr_thres);
        end
    end
    
    bw_stack = logical(bw_stack);
    
end


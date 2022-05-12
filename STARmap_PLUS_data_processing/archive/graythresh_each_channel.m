function [thres_stack, bw_stack] = graythresh_each_channel(input_stack)
%graythres_each_channel is used to apply global graythresholding to all
%channels, and output both a thresholded version and binary mask of the
%input stack
%   -----IO-----
%   input_stack: 5D image stack
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
            T = graythresh(curr_channel);
            curr_bw = imbinarize(curr_channel, T);
            curr_thres = curr_channel;
            curr_thres(~curr_bw) = 0;
            
            bw_stack(:,:,:,c,r) = curr_bw;
            thres_stack(:,:,:,c,r) = uint8(curr_thres);
        end
    end
    
    bw_stack = logical(bw_stack);
    
end


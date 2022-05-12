function [img_stack] = smooth_channels(img_stack, sigma, f_size)
%smooth_channels is used to equalize or normalize intensity profile
%among all channels in each round
%   -----IO-----
%   img_stack: 5D array XxYxZxChxRd
%   img_stack: Histogram matched img_stack
    if nargin < 2
        sigma = 1;
    end 
    
    if nargin < 3
        f_size = [5 5 3];
    end 
    
    Nround = size(img_stack, 5);
    Nchannel = size(img_stack, 4);
    for r=1:Nround 
        fprintf('Smooth each channel in round %d\n', r);
        currStack = uint8(img_stack(:,:,:,:,r));        
        for c=1:Nchannel              
            currStack(:,:,:,c) = imgaussfilt3(uint8(currStack(:,:,:,c)), sigma, 'FilterSize', f_size);
        end    
        img_stack(:,:,:,:,r) = currStack;
    end


end


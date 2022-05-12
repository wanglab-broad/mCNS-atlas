function [ output_img ] = MinMaxNorm( input_img )
%MinMaxNorm is used to normalize the intensity profile of each channel
%   -----IO-----
%   input_img: cell with input images 
%   output_img: cell with normalized images


    Nround = numel(input_img);
    output_img = cell(Nround, 1);
    
    for r=1:Nround
        fprintf(sprintf("Normalizing Round %d...\n", r))
        curr_round = input_img{r};
        
        Nchannel = size(curr_round, 4);
        for c=1:Nchannel 
            curr_channel = curr_round(:,:,:,c);
            
            curr_channel = im2double(curr_channel);
            curr_min = min(curr_channel, [] ,'all');
            curr_max = max(curr_channel, [] ,'all');
            
            
            curr_channel = (curr_channel - curr_min) ./ (curr_max - curr_min);
            curr_channel = uint8(curr_channel .* 255); 
            
            
            curr_round(:,:,:,c) = curr_channel;
        end
        output_img{r} = curr_round;
    end

end


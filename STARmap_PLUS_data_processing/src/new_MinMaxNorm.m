function input_img = new_MinMaxNorm( input_img )
%MinMaxNorm is used to normalize the intensity profile of each channel
%   -----IO-----
%   input_img: mat with input images 
%   output_img: mat with normalized images


    Nround = size(input_img, 5);
    Nchannel = size(input_img, 4);
    
    fprintf("====Min-Max intensity normalization====\n");
    
    for r=1:Nround
        tic
        fprintf(sprintf("Normalizing Round %d...", r))
        
        for c=1:Nchannel 
            curr_channel = input_img(:,:,:,c,r);
            
            curr_channel = im2double(curr_channel);
            curr_min = min(curr_channel, [] ,'all');
            curr_max = max(curr_channel, [] ,'all');
            
            
            curr_channel = (curr_channel - curr_min) ./ (curr_max - curr_min);
            curr_channel = uint8(curr_channel .* 255); %% WARNING
            
            
            input_img(:,:,:,c,r) = curr_channel;
        end
        fprintf(sprintf('[time = %.2f s]\n', toc));
    end

end


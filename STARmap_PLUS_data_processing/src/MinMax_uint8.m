function img = MinMax_uint8(img)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    img = im2double(img);
    curr_min = min(img, [] ,'all');
    curr_max = max(img, [] ,'all');
    
    
    img = (img - curr_min) ./ (curr_max - curr_min);
    img = uint8(img .* 255); %% WARNING

end


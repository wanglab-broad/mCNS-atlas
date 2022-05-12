function [pixelQual,baseCall] = RefineBaseCall(registed_img,cc_stats)
%UNTITLED RefineBaseCall
%   test

Nround = size(registed_img, 1);

pixelQual = zeros([size(registed_img{1}), Nround]);
baseCall = {};


for r=1:Nround
    
    curr_round = registed_img{r};
    curr_cc_img = zeros(size(curr_round));
    Nchannel = size(curr_cc_img, 4);
    
    for c=1:Nchannel
        curr_ch = curr_round(:,:,:,c);
        curr_cc_ch = curr_cc_img(:,:,:,c);
        curr_idx = cell2mat(cc_stats.VoxelIdxList);
        curr_cc_ch(curr_idx) = curr_ch(curr_idx);
        curr_cc_img(:,:,:,c) = curr_cc_ch;
    end
    
    curr_qual_matrix = CalculateQual(curr_cc_img, Nchannel);
    
    curr_qual_smooth = imgaussfilt3(curr_qual_matrix, 1);
    curr_max_qual = imregionalmax(curr_qual_smooth,26);

    curr_qual_stats = regionprops3(curr_max_qual, 'Centroid');
    % qual_centroid = qual_stats(qual_stats.Volume > 0, "Centroid"); % volume filtration
    curr_qual_centroid = int16(table2array(curr_qual_stats));
    
    baseCall{r} = curr_qual_centroid;
    pixelQual(:,:,:,r) = curr_qual_matrix;
end



end


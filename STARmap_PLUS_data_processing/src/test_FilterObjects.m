function output_img = test_FilterObjects( input_img, low_v, high_v, smooth)


bw_img = zeros(size(input_img), 'logical');
output_img = zeros(size(input_img), 'uint8');
[~, ~, Nslice, Nchannel, Nround] = size(input_img);

% Binarize
for r=1:Nround
    for c=1:Nchannel
        for z=1:Nslice
            if smooth
                curr_slice = imgaussfilt(input_img(:,:,z,c,r), 0.5);
            else
                curr_slice = input_img(:,:,z,c,r);
            end
            curr_bw = imbinarize(curr_slice);
            bw_img(:,:,z,c,r) = curr_bw;
        end
    end
end



for r=1:Nround
    for c=1:Nchannel
        curr_channel = bw_img(:,:,:,c,r);
        gray_img = input_img(:,:,:,c,r);
        CC = bwconncomp(curr_channel, 26);
        S = regionprops3(CC, 'Volume');
        L = labelmatrix(CC);
        curr_channel_out = ismember(L, find(([S.Volume] >= low_v) & ([S.Volume] <= high_v)));
        gray_img(~curr_channel_out) = 0;
        output_img(:,:,:,c,r) = gray_img;
    end
end


end


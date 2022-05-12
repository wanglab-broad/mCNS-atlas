function locations = SpotFindingMax3D( input_img, ref_index )
%SpotFindingMax3D 

    locations = [];

    ref_round = input_img(:,:,:,:,ref_index);
    Nchannel = size(ref_round, 4);
    
    for c=1:Nchannel
        curr_channel = ref_round(:,:,:,c);
        curr_max = imregionalmax(curr_channel);

        max_intensity = max(curr_channel, [], 'all');
        curr_out = curr_max & curr_channel > 0.2 * max_intensity;

        curr_centroid = regionprops3(curr_out, "Centroid");
        curr_centroid = int16(curr_centroid.Centroid);
        locations = [locations; curr_centroid]; %% WARNING
    end




end


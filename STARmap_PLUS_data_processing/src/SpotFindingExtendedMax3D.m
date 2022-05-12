function locations = SpotFindingExtendedMax3D( input_img )
%SpotFindingMax3D 

    locations = [];

    curr_round = input_img(:,:,:,:,1);
    Nchannel = size(curr_round, 4);
    
    for c=1:Nchannel
        curr_channel = curr_round(:,:,:,c);
        curr_max = imextendedmax(curr_channel, 100);

        max_intensity = max(curr_channel, [], 'all');
        curr_out = curr_max & curr_channel > 0.2 * max_intensity;

        curr_centroid = regionprops3(curr_out, "Centroid");
        curr_centroid = int16(curr_centroid.Centroid);
        locations = [locations; curr_centroid]; %% WARNING
    end




end


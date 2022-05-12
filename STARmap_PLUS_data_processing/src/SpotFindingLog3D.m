function locations = SpotFindingLog3D( input_img, fsize, fsigma )
%SpotFindingLog3D 

    locations = [];

    curr_round = input_img(:,:,:,:,1);
    Nchannel = size(curr_round, 4);
    
    log3d = fspecial3('log', fsize, fsigma);
    
    for c=1:Nchannel
        curr_channel = curr_round(:,:,:,c);
        curr_log = imfilter(curr_channel, log3d);
        curr_min = imregionalmin(curr_log);
        max_intensity = max(curr_channel, [], 'all');
        curr_out = curr_min & curr_channel > 0.2 * max_intensity;
        curr_centroid = regionprops3(curr_out, "Centroid");
        curr_centroid = int16(curr_centroid.Centroid);
        locations = [locations; curr_centroid]; %% WARNING
    end



end


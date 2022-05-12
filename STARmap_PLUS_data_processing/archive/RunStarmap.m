function [ good_points, read_counts ] = RunStarmap( input_path )
%URunStarmap

    Nround = 6;
    prefix = '';
    doNormalize = false;

    test = STARMapDataset(input_path, Nround, prefix);
    test = test.LoadRawImages(prefix, doNormalize, 4, [1:36]);
    test = test.MinMaxNormalize();
    test = test.MorphoRecon();
    test = test.SwapChannels2And3();
    Ntile = [1 1];

    test = test.RegisterImages(Ntile);
    test = test.LocalRegisterImages('method', "sum");

    max_centroid = [];

    curr_round = test.registeredImages{1};
    for c=1:4
        curr_channel = curr_round(:,:,:,c);
        curr_max = imregionalmax(curr_channel);

        max_intensity = max(curr_channel, [], 'all');
        curr_out = curr_max & curr_channel > 0.2 * max_intensity;

        curr_centroid = regionprops3(curr_out, "Centroid");
        curr_centroid = int16(curr_centroid.Centroid);
        max_centroid = [max_centroid; curr_centroid];
    end
    
    
    %thresholdRel = 0.2;
    scoreThres = 0.5;
    regionSize = [3 3 1];

    test = test.changePoints(max_centroid, regionSize, scoreThres);
    
    test = test.LoadCodebook();
    %test = test.LoadCodebook(false);

    % Filter
    endBases = ['C'; 'C']; % used for 16genes
    test = test.FilterReads(endBases);
    
    good_points = test.goodPoints;
    read_counts = test.SaveTopGenes;
    
end


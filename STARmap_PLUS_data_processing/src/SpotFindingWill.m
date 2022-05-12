function locations = SpotFindingWill( input_img, varargin )
%SpotFindingWill

    p = inputParser();
    p.addParameter('thresholdRel', 0.2, @isnumeric);
    p.addParameter('ksize', 3, @isnumeric);
    p.addParameter('kdims', [5 5], @isnumeric);
    p.parse(varargin{:});

    thresholdRel = p.Results.thresholdRel;
    ksize = p.Results.ksize;
    kdims = p.Results.kdims;
    
    locations = [];

    curr_round = input_img(:,:,:,:,1);
    Nchannel = size(curr_round, 4);
    
    for c=1:Nchannel
        currColor = curr_round(:,:,:,c);
        currCentroids = new_FindSpots2D(currColor, thresholdRel, kdims, ksize);
        locations = [locations; currCentroids];
    end       
    
    
end


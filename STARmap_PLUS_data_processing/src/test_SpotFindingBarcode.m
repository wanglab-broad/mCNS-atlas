function [allReads, allSpots, allScores, basecsMat] = test_SpotFindingBarcode( input_img, seqTogene, quality_factor, volume_threshold, showPlots, useGPU)
%SpotFindingBarcode

    Nround = size(input_img, 5);
    goodCodes = str2double(keys(seqTogene));
    allSpots = [];
    allReads = {};
    basecsMat = [];
    allScores = [];
    
    local_threshold = quality_factor;
    
    % get barcode for each pixel
    [base_call_stack, quality_matrix, round_quality_matrix, barcodes] = GetBarcodes(input_img, local_threshold, false); % 
    mean_quality_matrix = mean(round_quality_matrix, 4);
    
    % get global quailty
    base_call_qual = quality_matrix(quality_matrix ~= 0);
    base_call_qual = reshape(base_call_qual, 1, numel(base_call_qual));
    
    % get global quailty threshold
    quality_threshold = quality_factor * max(quality_matrix, [], 'all');
    fprintf('The global quality score threshold for barcode calling is: %4.2f \n', quality_threshold)
    
    if showPlots
        figure
        hold on
        plot(sort(base_call_qual, 'descend'))
        plot(1:numel(base_call_qual), ones(numel(base_call_qual),1) * mean(base_call_qual))
        title('Global barcode quality')
        hold off
    end
    
            

    % iterate through correct barcode clusters
    props = {'Volume'; 'Centroid'; 'VoxelIdxList'};
    % props = {'Volume'; 'Centroid'; 'EquivDiameter'; 'VoxelIdxList'};

    all_volume = []; 
    fprintf('Iterate through clusters...')
    tic;

    iter_codes = goodCodes;
%     iter_codes = barcodes;

    
    parfor c = 1:numel(iter_codes)

        curr_barcode = iter_codes(c); 

        curr_img = zeros(size(base_call_stack)); 
        curr_img(base_call_stack == curr_barcode & quality_matrix > quality_threshold) = 1; % global quality filtration 
        curr_img = logical(curr_img);


        curr_stats = regionprops3(curr_img, props);
        all_volume = [all_volume; curr_stats.Volume]; %% WARNING
        curr_valid_stats = curr_stats(curr_stats.Volume > volume_threshold, ["Centroid", "VoxelIdxList"]); % volume filtration
        curr_valid_cc = int16(curr_valid_stats.Centroid);
        if ~isempty(curr_valid_cc)    
            Nspots = size(curr_valid_cc, 1);
            allSpots = [allSpots; curr_valid_cc]; %% WARNING
            curr_barcode = Colorseq2Str(curr_barcode);
            curr_read = repelem({curr_barcode}, Nspots, 1);
            allReads = [allReads; curr_read];
            curr_mat = reshape(curr_barcode, [Nround, 1]);
            curr_mat = arrayfun(@str2double, curr_mat);
            curr_mat = repmat(curr_mat', [Nspots, 1]);
            basecsMat = [basecsMat; curr_mat];
            
            for i=1:Nspots
                curr_spot = curr_valid_stats.VoxelIdxList{i};
                curr_score = mean_quality_matrix(curr_spot);
                curr_score = mean(curr_score);
                allScores = [allScores; curr_score];
            end
        end

    end
    

    fprintf(sprintf('[time = %.2f s]\n', toc));


    if showPlots
        figure
        histogram(all_volume)
        title("All CC volume");

        figure
        histogram(all_volume(all_volume > volume_threshold))
        title(sprintf("CC volume largerd than %f", volume_threshold));
    end
    
            
    
end


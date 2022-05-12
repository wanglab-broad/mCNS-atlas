function locations = SpotFindingBarcode( input_img, seqTogene, quality_factor, volume_threshold, method, showPlots, useGPU)
%SpotFindingBarcode

    goodCodes = str2double(keys(seqTogene));
    locations = [];
    
    local_threshold = quality_factor;
    % get barcode for each pixel
    [base_call_stack, quality_matrix, ~, barcodes] = GetBarcodes(input_img, local_threshold, false); % 

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
    
    switch method
        case "iteration"
            
%             % fill missing value
%             for i =1:1
%             base_call_stack(base_call_stack == 0) = NaN;
%             base_call_stack = fillmissing(base_call_stack, 'movmedian', [1 1]);
%             base_call_stack = fillmissing(base_call_stack, 'movmedian', [1 1], 2);
%             base_call_stack = fillmissing(base_call_stack, 'constant', 0);
% 
%             quality_matrix(quality_matrix == 0) = NaN;
%             quality_matrix = fillmissing(quality_matrix, 'movmedian', [1 1]);
%             quality_matrix = fillmissing(quality_matrix, 'movmedian', [1 1], 2);
%             quality_matrix = fillmissing(quality_matrix, 'constant', 0);
%             end

            % iterate through all barcode clusters
            %Ncluster = numel(barcodes);
            props = {'Volume'; 'Centroid'};
            % props = {'Volume'; 'Centroid'; 'EquivDiameter'; 'VoxelIdxList'};

            all_volume = []; 
            fprintf('Iterate through clusters...')
            tic;
            
            iter_codes = barcodes;
            %iter_codes = goodCodes;
             
            parfor c = 1:numel(iter_codes)

                curr_barcode = iter_codes(c); 
               
                curr_img = zeros(size(base_call_stack)); 
                curr_img(base_call_stack == curr_barcode & quality_matrix > quality_threshold) = 1; % global quality filtration 
                curr_img = logical(curr_img);


                curr_stats = regionprops3(curr_img, props);
                all_volume = [all_volume; curr_stats.Volume]; %% WARNING
                curr_valid_cc = curr_stats(curr_stats.Volume > volume_threshold, "Centroid"); % volume filtration
                curr_valid_cc = int16(table2array(curr_valid_cc));
                if ~isempty(curr_valid_cc)       
                    locations = [locations; curr_valid_cc]; %% WARNING
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

        case "image"
            % NEW
            if useGPU
                base_call_stack = gather(base_call_stack);
            end

            out_stack = zeros(size(base_call_stack));
            se = ones(3);
            
            base_call_stack(quality_matrix <= quality_threshold) = 0; % global quality filtration
            for z=1:size(base_call_stack,3)
                curr_slice = base_call_stack(:,:,z);
                good_pixels = ismember(curr_slice, goodCodes);
                curr_slice(~good_pixels) = 0;
                
                %bmask = boundarymask(base_call_stack(:,:,z));
                bmask = imdilate(curr_slice, se) > curr_slice;
                
                %base_bw = imbinarize(base_call_stack(:,:,z));
                base_bw = imbinarize(curr_slice);
                base_bw = imfill(base_bw, 'holes');
                bmask = imcomplement(bmask);
                out_bw = base_bw .* bmask;
                out_stack(:,:,z) = out_bw;
                
%                 figure             
%                 imshow(out_bw)
            end

            props = {'Volume'; 'Centroid'};
            % out_stack(quality_matrix <= quality_threshold) = 0; % global quality filtration
            %out_stack( (0.9 < quality_matrix) & ( quality_matrix <= 1)) = 0; % test 
            curr_stats = regionprops3(logical(out_stack), props);
            curr_valid_cc = curr_stats(curr_stats.Volume > volume_threshold, "Centroid"); % volume filtration
            curr_valid_cc = int16(table2array(curr_valid_cc));
            if ~isempty(curr_valid_cc)       
                locations = [locations; curr_valid_cc]; %% WARNING
            end

            if showPlots
                figure
                histogram(curr_stats.Volume)
                title("All CC volume");

                figure
                histogram(curr_stats.Volume(curr_stats.Volume > volume_threshold))
                title(sprintf("CC volume largerd than %f", volume_threshold));
            end
            
            
            
        case "iteration_test"
            
            
            
            
            
    end
    
end


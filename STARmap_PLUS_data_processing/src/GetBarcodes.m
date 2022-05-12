function [base_call_stack, all_quality_stack, round_quality_stack, barcodes] = GetBarcodes( input_img, local_threshold, useGPU )
% GetBarcodes

    % Convert to double for calculation
    input_img = im2double(input_img);
    [dim_x, dim_y, dim_z, Nchannel, Nround] = size(input_img);
    
    if useGPU
        base_call_stack = zeros(dim_x, dim_y, dim_z, Nround, 'gpuArray');
        all_quality_stack = zeros(dim_x, dim_y, dim_z, Nround, 'gpuArray'); 
    else
        base_call_stack = zeros(dim_x, dim_y, dim_z, Nround);
        all_quality_stack = zeros(dim_x, dim_y, dim_z, Nround);
    end

    
    for r = 1:Nround
        tic;
        fprintf('Get barcodes for Round %d....', r)
        curr_round = input_img(:, :, :, :, r);
        [curr_barcode_matrix, curr_quality_matrix] = GetBarcodesEachRound(curr_round, Nchannel);
        
%         % fill missing value
%         for i =1:1
%         curr_barcode_matrix(curr_barcode_matrix == 0) = NaN;
%         curr_barcode_matrix = fillmissing(curr_barcode_matrix, 'movmedian', 1);
%         curr_barcode_matrix = fillmissing(curr_barcode_matrix, 'movmedian', 1, 2);
%         curr_barcode_matrix = fillmissing(curr_barcode_matrix, 'constant', 0);
% 
%         curr_quality_matrix(curr_quality_matrix == 0) = NaN;
%         curr_quality_matrix = fillmissing(curr_quality_matrix, 'movmedian', 1);
%         curr_quality_matrix = fillmissing(curr_quality_matrix, 'movmedian', 1, 2);
%         curr_quality_matrix = fillmissing(curr_quality_matrix, 'constant', 0);
%         end
        
        % local quality filtration
        idx = curr_quality_matrix >= local_threshold;

        curr_barcode_matrix(curr_quality_matrix < local_threshold) = 0;
        curr_quality_matrix(curr_quality_matrix < local_threshold) = 0;
        
        count = sum(idx, 'all');
        base_call_stack(:,:,:,r) = curr_barcode_matrix;
        all_quality_stack(:,:,:,r) = curr_quality_matrix;
        fprintf('Q: %.2f, C: %d', mean(curr_quality_matrix(curr_quality_matrix ~= 0), 'all'), count);
        fprintf(sprintf('[time = %.2f s]\n', toc));
    end
    
    filtration_mask = base_call_stack(:,:,:,1);
    for r=2:Nround
        filtration_mask = filtration_mask .* base_call_stack(:,:,:,r);
    end
    filtration_mask(filtration_mask ~= 0) = 1;
    
    % get all barcode permutations 
    barcodes = permn(1:Nchannel, Nround);
    
    % assemble barcode of each pixel
    [base_call_stack, barcodes] = AssembleBarcode(base_call_stack, barcodes, Nround);
    
    % get quality score for each pixel
    round_quality_stack = all_quality_stack;
    all_quality_stack = sum(all_quality_stack, 4);
    
    % filtration 
    base_call_stack = base_call_stack .* filtration_mask;
    all_quality_stack = all_quality_stack .* filtration_mask;
    
end
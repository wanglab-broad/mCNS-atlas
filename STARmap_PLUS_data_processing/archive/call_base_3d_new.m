function [base_call_stack, all_quality_stack, barcodes] = call_base_3d_new(input_stack)
    input_stack = im2double(input_stack);
    [dim_x, dim_y, dim_z, Nchannel, Nround] = size(input_stack);
    
    
    base_call_stack = zeros(dim_x, dim_y, dim_z, Nround);
    all_quality_stack = zeros(dim_x, dim_y, dim_z, Nround);
    
    fprintf('Start base calling....\n')
    for n = 1:Nround
        fprintf('Round %d....\n', n)
        curr_round = input_stack(:, :, :, :, n);
        [curr_barcode_matrix, curr_quality_matrix] = call_base_each_round(curr_round, Nchannel);
        
        %sum(curr_quality_matrix < 0.7, 'all')
        curr_barcode_matrix(curr_quality_matrix < 0.7) = 0;
        curr_quality_matrix(curr_quality_matrix < 0.7) = 0;
        
        base_call_stack(:,:,:,n) = curr_barcode_matrix;
        all_quality_stack(:,:,:,n) = curr_quality_matrix;
    end
    
    barcodes = permn(1:Nchannel, Nround);
    
    [base_call_stack, barcodes] = get_barcode(base_call_stack, barcodes);
    all_quality_stack = sum(all_quality_stack, 4);
    
end
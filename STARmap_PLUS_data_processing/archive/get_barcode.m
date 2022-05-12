function [output_stack, output_barcodes] = get_barcode(input_stack, input_barcodes)

    input_dim = size(input_stack);
    if numel(input_dim) < 4
        Nround = 1;
    else
         Nround = input_dim(4);
    end

    barcode_scalar = 10.^(0 : Nround - 1);
    
    output_stack = zeros(input_dim);
    output_barcodes = zeros(size(input_barcodes));
    
    for i = 1:Nround
        output_stack(:,:,:,i) = input_stack(:,:,:,i) * barcode_scalar(i);
        output_barcodes(:, i) = input_barcodes(:, i) * barcode_scalar(i);
    end
    
    output_stack = sum(output_stack, 4);
    output_barcodes = sum(output_barcodes, 2);
end
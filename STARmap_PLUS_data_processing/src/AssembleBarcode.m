function [input_stack, input_barcodes] = AssembleBarcode( input_stack, input_barcodes, Nround )

    %barcode_scalar = 10.^(0 : Nround - 1);
    barcode_scalar = 10.^(Nround - 1 : -1 : 0);
    
    for i = 1:Nround
        input_stack(:,:,:,i) = input_stack(:,:,:,i) * barcode_scalar(i);
        input_barcodes(:, i) = input_barcodes(:, i) * barcode_scalar(i);
    end
    
    input_stack = sum(input_stack, 4);
    input_barcodes = sum(input_barcodes, 2);
end
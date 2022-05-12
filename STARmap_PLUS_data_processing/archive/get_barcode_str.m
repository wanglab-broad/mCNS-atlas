function barcode_str = get_barcode_str(input_vector)
% Convert multiple num codes to one str
    barcode_str = '';
    for i = 1:numel(input_vector)
        barcode_str = [barcode_str num2str(input_vector(i))];
    end

end
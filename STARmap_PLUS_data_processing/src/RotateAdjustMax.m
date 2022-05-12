function output_image = RotateAdjustMax(input_image, rotation_degree, lower_pct, upper_pct, output_path)
%RotateAdjustMax 

    if nargin < 5
        output_path = [];
    end
    
    input_image = imrotate(input_image, rotation_degree);
    output_image = max(input_image, [], 3);
    output_image = MinMax_uint8(output_image);
    output_image = imadjust(output_image, stretchlim(output_image, [lower_pct upper_pct]), []);
    
    if ~isempty(output_path)
        imwrite(output_image, output_path);
    end
    
end


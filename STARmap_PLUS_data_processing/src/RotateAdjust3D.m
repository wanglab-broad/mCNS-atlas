function output_image = RotateAdjust3D(input_image, rotation_degree, lower_pct, upper_pct, output_path)
%RotateAdjustMax 

    if nargin < 5
        output_path = [];
    end
    
    output_image = imrotate(input_image, rotation_degree); % rotate the image
    current_limit = stretchlim(output_image(:), [lower_pct upper_pct]);
    output_image = imadjustn(output_image, current_limit);
  
    
    if ~isempty(output_path)
        if exist(output_path, 'file') == 2
            delete(output_path);
        end
        for j=1:size(output_image, 3)
            imwrite(output_image(:,:,j), output_path, 'writemode', 'append');        
        end
    end
    
end


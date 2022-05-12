function SaveSpotImages(obj, file_name, esize)
%

    spot_stack = get_spot_stack(obj.dims, obj.goodSpots, esize);

    output_path = obj.outputPath;
    file_path = fullfile(output_path, file_name);

    if exist(file_path, 'file') == 2
       delete(file_path);
    end
    
    for z=1:size(spot_stack, 3)        
        imwrite(spot_stack(:,:,z), file_path, 'writemode', 'append');        
    end

end


function SaveSingleTiff(input_img, filename)

    if exist(filename, 'file') == 2
        delete(filename);
    end

    for j=1:size(input_img, 3)
        imwrite(input_img(:,:,j), filename, 'writemode', 'append');        
    end
        
end


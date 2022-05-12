function SaveImages(img_stack, output_path, dir_name)
%save_thres_imgs is used to save images
%   -----IO-----
%   output_path: output location
%   img_stack: the thresholded image stack

    Nround = size(img_stack, 5);
    Nchannel = size(img_stack, 4);
    
    thres3DPath = fullfile(output_path, dir_name);
    if ~exist(thres3DPath, 'dir')
        mkdir(thres3DPath);
    end


    for r = 1:Nround
        curr_dir = fullfile(thres3DPath, sprintf('round%d',r));
        if ~exist(curr_dir, 'dir')
           mkdir(curr_dir);
        end
        thres_temp = uint8(img_stack(:,:,:,:,r));
        for c = 1:Nchannel
            fname = fullfile(curr_dir, sprintf('round%02d_ch%02d.tif',r,c));
            if exist(fname, 'file') == 2
                delete(fname);
            end
            for z=1:size(thres_temp,3)        
                imwrite(thres_temp(:,:,z,c), fname, 'writemode', 'append');        
            end
        end
    end

    
end


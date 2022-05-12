function save_thres_imgs(output_path, thres_stack)
%save_thres_imgs is used to save thresholded images
%   -----IO-----
%   output_path: output location
%   thres_stack: the thresholded image stack

    Nround = size(thres_stack, 5);
    Nchannel = size(thres_stack, 4);
    
    thres3DPath = fullfile(output_path, 'thres');
    if ~exist(thres3DPath, 'dir')
        mkdir(thres3DPath);
    end


    for r = 1:Nround
        curr_dir = fullfile(thres3DPath, sprintf('round%d',r));
        if ~exist(curr_dir, 'dir')
           mkdir(curr_dir);
        end
        thres_temp = uint8(thres_stack(:,:,:,:,r));
        for c = 1:Nchannel
            fname = fullfile(curr_dir, sprintf('thres_round%02d_ch%02d.tif',r,c));
            if exist(fname, 'file') == 2
                delete(fname);
            end
            for z=1:size(thres_temp,3)        
                imwrite(thres_temp(:,:,z,c), fname, 'writemode', 'append');        
            end
        end
    end

    
end


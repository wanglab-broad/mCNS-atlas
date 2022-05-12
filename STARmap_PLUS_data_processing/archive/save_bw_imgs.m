function save_bw_imgs(output_path, bw_stack)
%save_bw_imgs is used to save binary thresholded mask 
%   -----IO-----
%   output_path: the output location
%   bw_stack: the binary image stack 

    Nround = size(bw_stack, 5);
    Nchannel = size(bw_stack, 4);
    
    bw3DPath = fullfile(output_path, 'bw');
    if ~exist(bw3DPath, 'dir')
        mkdir(bw3DPath);
    end


    for r = 1:Nround
        curr_dir = fullfile(bw3DPath, sprintf('round%d',r));
        if ~exist(curr_dir, 'dir')
           mkdir(curr_dir);
        end
        % bw_temp = uint8(bw_stack(:,:,:,:,r));
        bw_temp = bw_stack(:,:,:,:,r);
        for c = 1:Nchannel
            fname = fullfile(curr_dir, sprintf('bw_round%02d_ch%02d.tif',r,c));
            if exist(fname, 'file') == 2
                delete(fname);
            end
            for z=1:size(bw_temp,3)        
                imwrite(bw_temp(:,:,z,c), fname, 'writemode', 'append');        
            end
        end
    end

    
end


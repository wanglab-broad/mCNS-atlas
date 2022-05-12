function new_SaveImg( dirName, imgs )

[~, ~, z, Nchannel, Nround] = size(imgs);

fprintf(sprintf('Saving images to %s\n', dirName));
tic;
for r=1:Nround
    curr_dir = fullfile(dirName, sprintf('round%d', r));
    if ~exist(curr_dir, 'dir')
       mkdir(curr_dir);
    end
    temp = uint8(imgs(:,:,:,:,r));
    for ch=1:Nchannel
        fname = fullfile(curr_dir, sprintf('reg_round%02d_ch%02d.tif', r, ch));
        if exist(fname, 'file') == 2
            delete(fname);
        end
        for j=1:z        
            imwrite(temp(:,:,j,ch), fname, 'writemode', 'append');        
        end
    end
    fprintf(sprintf('Round %d finished...[time=%02f]\n', r, toc));
end


function SaveCellImg( dirName, imgs, curr_data_dir, sub_dirs )

Nchannel = numel(imgs);
z = size(imgs{1}, 3);


fprintf(sprintf('Saving cell images to %s\n', dirName));
tic;

%sub_dirs = ["plaque", "tau", "pi", "merged_dots"];

for ch=1:Nchannel
    ch_dir = fullfile(dirName, sub_dirs(ch));
    if ~exist(ch_dir, 'dir')
       mkdir(ch_dir);
    end
    fname = fullfile(ch_dir, sprintf('%s.tif', curr_data_dir));
    if exist(fname, 'file') == 2
        delete(fname);
    end
    for j=1:z        
        imwrite(imgs{ch}(:,:,j), fname, 'writemode', 'append');        
    end
end
fprintf(sprintf('Finished!...[time=%02f]\n', toc));



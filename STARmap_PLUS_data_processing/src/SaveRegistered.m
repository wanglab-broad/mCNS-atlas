function  SaveRegistered(dirName, imgs )

%SAVECOLORMAX Summary of this function goes here
%   Detailed explanation goes here
for i=1:numel(imgs)
    %i
    curr_dir = fullfile(dirName, sprintf('round%d',i));
    if ~exist(curr_dir, 'dir')
       mkdir(curr_dir);
    end
    temp = uint8(imgs{i});
    for ch=1:4
        fname = fullfile(curr_dir, sprintf('reg_round%02d_ch%02d.tif',i,ch));
        if exist(fname, 'file') == 2
            delete(fname);
        end
        for j=1:size(temp,3)        
            imwrite(temp(:,:,j,ch), fname, 'writemode', 'append');        
        end
    end
end


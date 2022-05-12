function images = LoadImageStacks( input_dir, prefix_name, Nround, Nchannel, zrange)
%LOADIMAGESTACKS Load image stacks for each round
warning('off','all');
if nargin < 4
    Nchannel = 4;
end

if nargin < 5
    zrange = [];
end

images = cell(Nround,1);
for i=1:Nround
    
    fprintf('Loading round %d...\n', i);
    
    %dirname = fullfile(input_dir);
    %dirname = fullfile(input_dir, ['round' num2str(i) '/'], 'uint8/');
    dirname = fullfile(input_dir, ['round' num2str(i) '/']);
    
    fnames = cell(5,1);
    
    finfo = dir(strcat(dirname, '*.tif'));
    
    for j=1:Nchannel % 4 colors without DAPI
        finfo(j).name
        %fnames{j} = fullfile(dirname, sprintf('%sseq%d_decon_ch0%d.tif', prefix_name, i, j-1));
        fnames{j} = fullfile(dirname, finfo(j).name);
        %fnames{j} = fullfile(dirname, sprintf('ch0%d.tif', j-1));
        %fnames{j} = fullfile(dirname, sprintf('%sseq%d_decon_crop_ch0%d.tif', prefix_name, i, j-1));
        %fnames{j} = fullfile(dirname, sprintf('%sseq%d_decon_decon_crop_ch0%d.tif', prefix_name, i, j-1));
    end
    
    % Load channel 1 
    currImg = uint8(LoadMultipageTiff(fnames{1}));
    
    if ~isempty(zrange)
        currImg = currImg(:,:,zrange);
    end
    
    currImgs = zeros(size(currImg,1), size(currImg,2), size(currImg,3), 4,'uint8');
    currImgs(:,:,:,1) = currImg;
    
    % Load channel 2-4
    for j=2:Nchannel 
        currFrame = uint8(LoadMultipageTiff(fnames{j}));
        if isempty(zrange)
            currImgs(:,:,:,j) = currFrame;
        else
            currImgs(:,:,:,j) = currFrame(:,:,zrange);
        end
    end
    
    images{i} = currImgs;
end

% collapse to common sized array
maxX = 1E10; maxY = 1E10;
for i=1:Nround
   [currX,currY,~] = size(images{i}); 
   if currX < maxX
       maxX = currX;
   end
   if currY < maxY
       maxY = currY;
   end
end



fprintf('Collapsed to size %d by %d by %d\n', maxX, maxY, size(images{1},3));
for i=1:Nround
    fprintf('Collapsing round %d\n', i);
    temp = images{i};
    images{i} = temp(1:maxX, 1:maxY, :,:);
end


end


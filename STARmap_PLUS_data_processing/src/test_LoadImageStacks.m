function [output_imgs, dims] = test_LoadImageStacks( inputPath, sub_dir, input_dim, input_format, zrange, output_class, useGPU )
%LOADIMAGESTACKS Load image stacks for each round
    
        
    % Suppress all warnings 
    warning('off','all');
    
    % Get directories containing all images 
    dirs = dir(strcat(inputPath, 'round*'));
    Nround = numel(dirs);
    Nchannel = 4;
    
    
    switch output_class
        case "mat"
            
            if ~isempty(input_dim)
                output_imgs = zeros(input_dim, 'uint8');
            else
                output_imgs = [];
            end
                
            for r=1:Nround
                tic
                fprintf('Loading round %d...', r);

                curr_dir = dirs(r).name;

                curr_files = dir(fullfile(inputPath, curr_dir, sub_dir, '*.tif'));

%                 if ~isempty(input_dim)
%                     round_img = zeros(input_dim(1:4), 'uint8');
%                 else
%                     round_img = [];
%                 end

                % Load all channels
                for c=1:Nchannel 
                    curr_path = strcat(curr_files(c).folder, '/', curr_files(c).name);

                    curr_img = new_LoadMultipageTiff(curr_path, input_format, 'uint8', useGPU);
                    
                    if ~isempty(zrange)
                        curr_img = curr_img(:,:,zrange(1):zrange(2));
                    end
                    
%                     % append blank layer if z doesn't match 
%                     while size(curr_img, 3) ~= input_dim(3)
%                         curr_img(:,:,end+1) = zeros(input_dim(1:2), 'uint8');
%                     end
                    
                    output_imgs(:,:,:,c,r) = curr_img;
                end

                fprintf(sprintf('[time = %.2f s]\n', toc));
            end

%             % Collapse to common sized array
%             maxX = 1E10; maxY = 1E10;
%             for r=1:Nround
%                 curr_round = output_imgs(:,:,:,:,r);
%                [currX, currY, ~] = size(curr_round); 
%                if currX < maxX
%                    maxX = currX;
%                end
%                if currY < maxY
%                    maxY = currY;
%                end
%             end
% 
%             if strcmp(zrange, '')
%                zrange = 1:size(curr_round, 3);
%             end
% 
%             dims = [maxX maxY numel(zrange) size(curr_round, 4) Nround];
              dims = size(output_imgs);
% 
% 
%             % Show message for re-sizing
%             fprintf('Collapsed to size %d by %d by %d\n', maxX, maxY, numel(zrange));
% 
%             for r=1:Nround
%                 fprintf('Collapsing round %d\n', r);
%                 curr_round = output_imgs(:,:,:,:,r);
%                 if isempty(zrange)
%                     output_imgs(:,:,:,:,r) = curr_round(1:maxX, 1:maxY, :, :);
%                 else
%                     output_imgs(:,:,:,:,r) = curr_round(1:maxX, 1:maxY, zrange, :);
%                 end
%             end
            
            fprintf('Raw image as 5-D array\n');
            
            
        case "cell"
            % Set output_imgs
            output_imgs = cell(Nround, 1);
            
            for r=1:Nround
                tic
                fprintf('Loading round %d...', r);

                curr_dir = dirs(r).name;

                curr_files = dir(fullfile(inputPath, curr_dir, sub_dir, '*.tif'));

                if ~isempty(input_dim)
                    round_img = zeros(input_dim(1:4), 'uint8');
                else
                    round_img = [];
                end

                % Load all channels
                for c=1:Nchannel 
                    curr_path = strcat(curr_files(c).folder, '/', curr_files(c).name);

                    curr_img = new_LoadMultipageTiff(curr_path, 'uint8', 'uint8', useGPU);
                    
                    round_img(:,:,:,c) = curr_img;
                end

                output_imgs{r} = round_img;
                fprintf(sprintf('[time = %.2f s]\n', toc));
            end

            % Collapse to common sized array
            maxX = 1E10; maxY = 1E10;
            for r=1:Nround
                curr_round = output_imgs{r};
               [currX, currY, ~] = size(curr_round); 
               if currX < maxX
                   maxX = currX;
               end
               if currY < maxY
                   maxY = currY;
               end
            end

            if strcmp(zrange, '')
               zrange = 1:size(curr_round, 3);
            end

            dims = [maxX maxY numel(zrange) size(curr_round, 4) Nround];


            % Show message for re-sizing
            fprintf('Collapsed to size %d by %d by %d\n', maxX, maxY, numel(zrange));

            for r=1:Nround
                fprintf('Collapsing round %d\n', r);
                curr_round = output_imgs{r};
                if isempty(zrange)
                    output_imgs{r} = curr_round(1:maxX, 1:maxY, :, :);
                else
                    output_imgs{r} = curr_round(1:maxX, 1:maxY, zrange, :);
                end
            end
            
            fprintf('Raw image as cell');
    
    
    
    
    end

    
end


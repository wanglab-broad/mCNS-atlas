function output_stack = load_raw_imgs(input_path)
%LOADIMAGESTACKS Load image stacks for each round
%   -----IO-----
%   input_path: obj.path
%   output_stack: 5D image stack

    % Suppress all warnings 
    warning('off','all');
    
    dirs = dir(strcat(input_path, 'round*'));
    Nround = numel(dirs);
    
    
    files = dir(strcat(input_path, dirs(1).name,'/*.tif'));
    Nchannel = numel(files);
    
    temp_cell = {Nround};
    
    for r=1:Nround
        
        fprintf('Loading round %d...\n', r);

        curr_dir = dirs(r).name;

        curr_files = dir(strcat(input_path, curr_dir,'/*.tif'));

        % Load channel 1 
        ch1_path = strcat(curr_files(1).folder, curr_files(1).name);
        
        curr_img = uint8(LoadMultipageTiff(ch1_path));
        

        if ~isempty(zrange)
            curr_img = curr_img(:,:,zrange);
        end

        curr_round = zeros(size(curr_img,1), size(curr_img,2), size(curr_img,3), Nchannel, 'uint8');
        curr_round(:,:,:,1) = curr_img;

        % Load rest channels
        for c=2:Nchannel 
            curr_path = strcat(curr_files(c).folder, curr_files(c).name);
            curr_img = uint8(LoadMultipageTiff(curr_path));
            
            if isempty(zrange)
                curr_round(:,:,:,c) = curr_img;
            else
                curr_round(:,:,:,c) = curr_img(:,:,zrange);
            end
        end

        temp_cell{r} = curr_round;
    end

    % Collapse to common sized array
    maxX = 1E10; maxY = 1E10;
    for r=1:Nround
       [currX,currY,~] = size(temp_cell{r}); 
       if currX < maxX
           maxX = currX;
       end
       if currY < maxY
           maxY = currY;
       end
    end


    fprintf('Collapsed to size %d by %d by %d\n', maxX, maxY, size(temp_cell{1},3));
    for r=1:Nround
        fprintf('Collapsing round %d\n', r);
        temp = temp_cell{r};
        temp_cell{r} = temp(1:maxX, 1:maxY, :,:);
    end

    output_stack = im_cell2mat(temp_cell);
end


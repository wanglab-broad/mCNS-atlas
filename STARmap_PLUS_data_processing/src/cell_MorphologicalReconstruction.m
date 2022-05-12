function input_img = cell_MorphologicalReconstruction( input_img, method, varargin )
%MorphologicalReconstruction

    % Input parser
    p = inputParser;
    
    % Defaults
    defaultRadius = 2; % org = 6
    defaultHeight = 3;


    addRequired(p,'input_img');
    addRequired(p,'method');
    addOptional(p,'radius',defaultRadius);
    addOptional(p,'height',defaultHeight);


    parse(p, input_img, method, varargin{:});

    % p.Results
    
    % get dims 
    Nround = numel(input_img);
    [~, ~, Nslice, Nchannel] = size(input_img{1});
    
    % p.Results.radius
    % switch case 2D/3D
    
    switch method
        case "2d"
            % setup structure element
            se = strel('disk', p.Results.radius);

            for r=1:Nround
                tic
                fprintf(sprintf("Processing Round %d...", r));

                for c=1:Nchannel
                    curr_channel = input_img{r}(:,:,:,c);
                    % curr_channel = input_img(:,:,:,c,r);
                    for z=1:Nslice
                        curr_slice = curr_channel(:,:,z);
                        marker = imerode(curr_slice, se); % Morphological opening is useful for removing small objects from an image while preserving the shape and size of larger objects in the image
                        obr = imreconstruct(marker, curr_slice);
                        curr_out = curr_slice - obr;
                        %mask = imbinarize(curr_out);
%                         mask = curr_out > 0;
%                         mask = bwareaopen(mask, 1);
%                         curr_out(~mask) = 0;

                        %bw = imopen(mask, se_2);
                        %curr_out(~bw) = 0;

                        curr_out = imsubtract(imadd(curr_out, imtophat(curr_out, se)), imbothat(curr_out, se));
                        
                        
                        curr_channel(:,:,z) = curr_out;
                    end
                    input_img{r}(:,:,:,c) = uint8(curr_channel);
                    % input_img(:,:,:,c,r) = uint8(curr_channel);
                end
                fprintf(sprintf('[time = %.2f s]\n', toc));
            end 
   
        case "2d_thres"
            % setup structure element
            se = strel('disk', p.Results.radius);

            for r=1:Nround
                tic
                fprintf(sprintf("Processing Round %d...", r));

                for c=1:Nchannel
                    curr_channel = input_img{r}(:,:,:,c);
                    % curr_channel = input_img(:,:,:,c,r);
                    for z=1:Nslice
                        curr_slice = curr_channel(:,:,z);
                        marker = imerode(curr_slice, se); % Morphological opening is useful for removing small objects from an image while preserving the shape and size of larger objects in the image
                        obr = imreconstruct(marker, curr_slice);
                        curr_out = curr_slice - obr;
                        
                        curr_bw = curr_out > 0 & curr_out < 10; % Threshold
                        curr_out(curr_bw) = 0;
                        
                        curr_out = im2double(curr_out);
                        curr_max = max(curr_out, [], 'all');
                        curr_min = min(curr_out, [], 'all');
                        curr_out = (curr_out - curr_min) ./ (curr_max - curr_min);
                        curr_out = uint8(curr_out .* 255); %% WARNING
    
                        curr_channel(:,:,z) = curr_out;
                    end
                    input_img{r}(:,:,:,c) = uint8(curr_channel);
                    % input_img(:,:,:,c,r) = uint8(curr_channel);
                end
                fprintf(sprintf('[time = %.2f s]\n', toc));
            end 
            
    end
    
end


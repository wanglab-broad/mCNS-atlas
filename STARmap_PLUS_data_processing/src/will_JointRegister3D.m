function input_img = will_JointRegister3D( input_img, ref_idx, nblocks, useOverlay, log_file )
% 

    if nargin < 2    
        nblocks = [1 1];    
    end

    if nargin < 3
        useOverlay = false;
    end

      
    % start time counting
    tic;

    % get dims 
    % Nround = size(input_img, 5);
    % Nround = Nround(1);
    [m, n, ~, Nchannel, Nround] = size(input_img);
    
    % blocks
    nblocksX = nblocks(1);
    nblocksY = nblocks(2);
    
    ol = 0.1; % 10% overlap
    stepsize_m = floor(m ./ nblocksX);
    mod_m = mod(m, nblocksX);
    stepsize_n = floor(n ./ nblocksY);
    mod_n = mod(n, nblocksY);
    fprintf('Using block size of %f x %f \n', stepsize_m, stepsize_n);
    % Compute overlap number of pixels
    overlap_m = floor(stepsize_m .* ol);
    overlap_n = floor(stepsize_n .* ol);

%     if useOverlay
%         ref_round = sum(input_img{ref_idx} / 4, 4);
%     else
%         ref_round = max(input_img{ref_idx}, [], 4);
%     end
    
    if useOverlay
        ref_round = sum(input_img(:,:,:,:,ref_idx) / 4, 4);
    else
        ref_round = max(input_img(:,:,:,:,ref_idx), [], 4);
    end
    
    % reg image
    regImgs = cell(Nround, 1);
    regImgs{1} = uint8(ref_round);
    
    % compute transform
    % use overlay
    if useOverlay
        params = cell(Nround-1, 1);
        for r=2:Nround
            % curr_overlay = sum(input_img{r} / 4, 4);
            curr_overlay = sum(input_img(:,:,:,:,r) / 4, 4);
            msg = sprintf('Registering Round %d vs. Round %d...[time=%02f]\n', r, ref_idx, toc);
            fprintf(msg);
            fprintf(log_file, msg);
            [params{r-1}, currReg] = DFTRegister3D(ref_round, curr_overlay, false);
%             disp(['Shifting by ' num2str(params{r-1}.shifts)]);
            fprintf(sprintf('Shifted by %s\n', num2str(params{r-1}.shifts)));
            fprintf(log_file, sprintf('Shifted by %s\n', num2str(params{r-1}.shifts)));
            regImgs{r} = uint8(currReg);
        end

    else % use max 
        params = cell(Nround-1, 1);
        for r=2:Nround
            % curr_max = max(input_img{r}, [], 4);
            curr_max = max(input_img(:,:,:,:,r), [], 4);
            msg = sprintf('Registering Round %d vs. Round %d...[time=%02f]\n', r, ref_idx, toc);
            fprintf(msg);      
            fprintf(log_file, msg);
            [params{r-1}, currReg] = DFTRegister3D(ref_round, curr_max, false);
%             disp(['Shifting by ' num2str(params{r-1}.shifts)]);
            fprintf(sprintf('Shifted by %s\n', num2str(params{r-1}.shifts)));
            fprintf(log_file, sprintf('Shifted by %s\n', num2str(params{r-1}.shifts)));
            regImgs{r} = uint8(currReg);
        end
    end

    % compute transform for blocks 
    if nblocksX > 1 || nblocksY > 1
        if useOverlay
            blockParams = cell(Nround, nblocksX, nblocksY);    
            for t = 2:Nround
                fprintf('Registering Round %d...\n', t);
                % Loop through each block generating the correct subset and registration
                for i = 1:nblocksX
                    for j = 1:nblocksY
                        % Store m block range
                        [bin_m, bin_n] = GetBlocks(i,j);
                        % Compute fourier transforms
                        block1 = ref_round(bin_m, bin_n, :);
                        block2 = regImgs{t}(bin_m, bin_n, :);
                        % Compute offset for block
                        [blockParams{t-1,i,j},~] = DFTRegister3D(block1, block2, false);
                        % Store output
                        disp(['Shifting (' num2str(i) ' ' num2str(j) ') by ' num2str(blockParams{t-1,i,j}.shifts)]);
                    end
                end
            end        
        else
            blockParams = cell(Nround-1, nblocksX, nblocksY);    

            for t = 2:Nround
                fprintf('Registering Round %d...\n', t);
                % Loop through each block generating the correct subset and registration
                for i = 1:nblocksX
                    for j = 1:nblocksY
                        % Store m block range
                        [bin_m, bin_n] = GetBlocks(i,j);
                        % Compute fourier transforms
                        block1 = regImgs{1}(bin_m, bin_n, :);
                        block2 = regImgs{t}(bin_m, bin_n, :);
                        % Compute offset for block
                        [blockParams{t-1, i, j}, ~] = DFTRegister3D(block1, block2, false);
                        % Store output
                        disp(['Shifting (' num2str(i) ' ' num2str(j) ') by ' num2str(blockParams{t-1,i,j}.shifts)]);
                    end
                end
            end
        end
    end

    % apply the transformation
    % use max 
    for r=2:Nround
        msg = sprintf('Applying registration to Round %d...[time=%02f]\n', r, toc);
        fprintf(msg);
        for c=1:Nchannel    
            %msg = sprintf('Applying registration to round [%d / %d], channel %d [time=%02f]\n', r, Nround, c, toc);
            %fprintf(msg);        
            %input_img{r}(:,:,:,c) = uint8(DFTApply3D(input_img{r}(:,:,:,c), params{r-1}, false));
            input_img(:,:,:,c,r) = uint8(DFTApply3D(input_img(:,:,:,c,r), params{r-1}, false));
        end        
    end
    
    % apply the transformation to blocks
    if nblocksX > 1 || nblocksY > 1
        for t = 2:Nround
            msg = sprintf('Applying registration to Round %d... [time=%02f]\n', t, toc);
            fprintf(msg);
            % Loop through each block generating the correct subset and registration
            for kk=1:Nchannel
                for i = 1:nblocksX
                    for j = 1:nblocksY
                        [bin_m, bin_n] = GetBlocks(i, j);
                        % Compute fourier transforms
                        % Compute offset for block
                        % input_img{t}(bin_m, bin_n, :, kk) = uint8(DFTApply3D(input_img{t}(bin_m, bin_n, :, kk), blockParams{t-1, i, j})); 
                        input_img(bin_m, bin_n, :, kk, t) = uint8(DFTApply3D(input_img(bin_m, bin_n, :, kk, t), blockParams{t-1, i, j}));  
                        %disp(['Shifting (' num2str(i) ' ' num2str(j) ')...']);
                    end
                end
            end
        end    
    end
    
    
function [bin_m, bin_n] = GetBlocks(i, j)
% Store m block range
if nblocksX > 1
    if i == 1
        bin_m = 1:(stepsize_m + overlap_m);
    elseif i == nblocksX
        bin_m = (stepsize_m.*(i-1) - overlap_m):(stepsize_m.*i + mod_m);
    else
        bin_m = (stepsize_m.*(i-1) - overlap_m):(stepsize_m.*i + overlap_m);
    end
elseif nblocksX == 1
    bin_m = 1:stepsize_m;
end

if nblocksY > 1
    % Store n block range
    if j == 1
        bin_n = 1:(stepsize_n + overlap_n);
    elseif j == nblocksY
        bin_n = (stepsize_n.*(j-1) - overlap_n):(stepsize_n.*j + mod_n);
    else
        bin_n = (stepsize_n.*(j-1) - overlap_n):(stepsize_n.*j + overlap_n);
    end
elseif nblocksY == 1
    bin_n = 1:stepsize_n;
end
    
end

end



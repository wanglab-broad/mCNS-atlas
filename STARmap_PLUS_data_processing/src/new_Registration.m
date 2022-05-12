function obj = new_Registration( obj, ref_idx, nblocks, useOverlay, Iterations, AFS )
%

    obj.registeredImages = obj.rawImages;
    
    if obj.useGPU
        obj.registeredImages{ref_idx} = gpuArray(obj.registeredImages{ref_idx});
    end
    
    if useOverlay
        method = "overlay";
    else
        method = "max";
    end
    
    ref_round = GetMovingImg(obj.registeredImages{ref_idx}, method);
    
    % start time counting
    tic;

    % get dims 
    Nround = obj.Nround;
    rounds = 1:Nround;
    rounds = rounds(rounds ~= ref_idx);
    [m, n, z, Nchannel] = size(obj.rawImages{1});
    
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

    
    % reg image
    regImgs = cell(Nround, 1);
    regImgs{ref_idx} = uint8(ref_round);
    
    global_params = {};
    blockParams = cell(Nround, nblocksX, nblocksY);
    
    pyd_level = floor(log2(z)); 
    
    for r=rounds
        
        % GPU
        if obj.useGPU
           obj.registeredImages{r} = gpuArray(obj.registeredImages{r});      
        end
       
        curr_round = GetMovingImg(obj.registeredImages{r}, method);
       
        % Global registration
        % compute transform
        msg = sprintf('(G) Registering Round %d vs. Round %d...[time=%02f]\n', r, ref_idx, toc);
        fprintf(msg);
        [global_params{r}, currReg] = DFTRegister3D(ref_round, curr_round, false);
        % global_params{r} = DFTRegister3D(ref_round, curr_round, false);
        disp(['Shifting by ' num2str(global_params{r}.shifts)]);
        regImgs{r} = uint8(currReg);
        
        if nblocksX > 1 || nblocksY > 1
            fprintf('(GB) Registering Round %d...\n', r);
            % Loop through each block generating the correct subset and registration
            for i = 1:nblocksX
                for j = 1:nblocksY
                    % Store m block range
                    [bin_m, bin_n] = GetBlocks(i,j);
                    % Compute fourier transforms
                    block1 = ref_round(bin_m, bin_n, :);
                    block2 = regImgs{r}(bin_m, bin_n, :);
                    % Compute offset for block
                    [blockParams{r,i,j},~] = DFTRegister3D(block1, block2, false);
                    % Store output
                    disp(['Shifting (' num2str(i) ' ' num2str(j) ') by ' num2str(blockParams{r,i,j}.shifts)]);
                end
            end
        end
        
        % apply transform
        msg = sprintf('(G) Applying registration to Round %d...[time=%02f]\n', r, toc);
        fprintf(msg);
        for c=1:Nchannel       
            obj.registeredImages{r}(:,:,:,c) = uint8(DFTApply3D(obj.registeredImages{r}(:,:,:,c), global_params{r}, false));
        end  
        
        % apply the transformation to blocks
        if nblocksX > 1 || nblocksY > 1
            msg = sprintf('(GB) Applying registration to Round %d... [time=%02f]\n', r, toc);
            fprintf(msg);
            % Loop through each block generating the correct subset and registration
            for kk=1:Nchannel
                for i = 1:nblocksX
                    for j = 1:nblocksY
                        [bin_m, bin_n] = GetBlocks(i, j);
                        % Compute fourier transforms
                        % Compute offset for block
                        obj.registeredImages{r}(bin_m, bin_n, :, kk) = uint8(DFTApply3D(obj.registeredImages{r}(bin_m, bin_n, :, kk), blockParams{r, i, j}));    
                        %disp(['Shifting (' num2str(i) ' ' num2str(j) ')...']);
                    end
                end
            end   
        end
        
        
%         % Non-rigid registration
%         currReg = GetMovingImg(obj.registeredImages{r}, method);
%         
%         % Block
%         [D, ~] = imregdemons(currReg, ref_round, Iterations, ...
%             'PyramidLevels', pyd_level, ...
%             'AccumulatedFieldSmoothing', AFS, ...
%             'DisplayWaitbar', false);
% 
%         if obj.useGPU
%             d = gather(D);
%             obj.registeredImages{r} = gather(obj.registeredImages{r});
%         else
%             d = D;
%         end
% 
%         % Apply displacement field on each channel
%         for c=1:Nchannel
%             obj.registeredImages{r}(:,:,:,c) = imwarp(obj.registeredImages{r}(:,:,:,c), d);
%         end
%             
% 
%         fprintf(sprintf('[time = %.2f s]\n', toc));
        
          
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


function MovingImg = GetMovingImg(InputImg, method)
    switch method
        case "max"
            MovingImg = max(InputImg, [], 4);
        case "overlay"
            MovingImg = sum(InputImg / 4, 4);
    end
end






end


function output_cell = JointRegister2D( input_cell, nblocks )
% given cell array of images, register multiple images to reference (last
% image in array)

if nargin < 2    
    nblocks = [1 1];    
end

nblocksX = nblocks(1);
nblocksY = nblocks(2);

output_cell = cell(size(input_cell));
Nround = numel(input_cell);

tic;

ref_round = input_cell{1};
ref_img = ref_round{1};
output_cell{1} = ref_round;

for r=2:Nround
    curr_round = input_cell{r};
    curr_dapi = curr_round{1};
    
    ol = 0.1;
    
    regImgs = cell(numel(curr_round),1);
    [m,n] = size(curr_dapi);
    Nimg = numel(curr_round);


    stepsize_m = floor(m./nblocksX);
    mod_m = mod(m,nblocksX);
    stepsize_n = floor(n./nblocksY);
    mod_n = mod(n,nblocksY);
    fprintf('Using block size of %f x %f \n', stepsize_m, stepsize_n);
    % Compute overlap number of pixels
    overlap_m = floor(stepsize_m.*ol);
    overlap_n = floor(stepsize_n.*ol);

    
    fprintf('%s\n', repmat('-',1,20));

    
 
    msg = sprintf('Registering round [%d / %d] [time=%02f]\n', r, Nround, toc);
    fprintf(msg);          
    [curr_param, currReg] = DFTRegister2D(ref_img, curr_dapi, false);

    disp(['Shifting by ' num2str(curr_param.shifts)]);

    regImgs{1} = uint8(currReg);

    if nblocksX > 1 || nblocksY > 1

        blockParams = cell(nblocksX, nblocksY);    

        % Loop through each block generating the correct subset and registration
        for i = 1:nblocksX
            for j = 1:nblocksY
                % Store m block range
                [bin_m, bin_n] = GetBlocks(i,j);
                % Compute fourier transforms
                block1 = ref_img(bin_m, bin_n,:);
                block2 = curr_dapi(bin_m, bin_n,:);
                % Compute offset for block
                [blockParams{i,j},~] = DFTRegister2D(block1, block2, false);

                % Store output
                disp(['Shifting (' num2str(i) ' ' num2str(j) ') by ' num2str(blockParams{i,j}.shifts)]);
            end
        end

    end
    
    % Apply to other images 
    for i=2:Nimg
        temp = curr_round{i};    
   
        msg = sprintf('Applying registration to image %d in current round [time=%02f]\n', i, toc);
        fprintf(msg);        
        temp = uint8(DFTApply2D(temp, curr_param, false));
        regImgs{i} = temp;
    end
    
    if nblocksX > 1 || nblocksY > 1
        for z = 2:Nimg
            for i = 1:nblocksX
                for j = 1:nblocksY
                    [bin_m, bin_n] = GetBlocks(i,j);
                    % Compute fourier transforms
                    % Compute offset for block
                    temp(bin_m, bin_n) = uint8(DFTApply2D(temp(bin_m, bin_n), blockParams{i,j}));                  
                end
            end
            regImgs{z} = temp;
        end
    end
    
    output_cell{r} = regImgs;
end



function [bin_m, bin_n] = GetBlocks(i,j)
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




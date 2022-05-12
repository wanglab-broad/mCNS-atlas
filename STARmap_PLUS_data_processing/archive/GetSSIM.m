function [ ssim_array_2d, ssim_array_3d, maps ] = GetSSIM( input_images, ref_round)
%GetSSIM
%   -----IO-----
%   input_images = input cell with images 
%   ref_round = the reference round 

    ref = input_images{ref_round};
    ref_3d = max(ref, [], 4);
    ref_2d = max(ref_3d, [], 3);
    
    Nround = numel(input_images);
    ssim_array_2d = [];
    ssim_array_3d = [];
    maps = cell(Nround);
    
    fprintf("-----SSIM--2D--3D-----\n");
    for r=1:Nround
        curr_round = input_images{r};
        curr_round_3d = max(curr_round, [], 4);
        curr_round_2d = max(curr_round_3d, [], 3);
        
        [curr_ssim_2d, map_2d] = ssim(curr_round_2d, ref_2d, 'Exponents', [0 0 1]);
        curr_ssim_3d = ssim(curr_round_3d, ref_3d, 'Exponents', [0 0 1]);
        
        ssim_array_2d = [ssim_array_2d; curr_ssim_2d];
        ssim_array_3d = [ssim_array_3d; curr_ssim_3d];
        maps{r} = map_2d;
        
        fprintf(sprintf("Round %d vs. Round %d: %.3f -- %.3f\n", r, ref_round, curr_ssim_2d, curr_ssim_3d));
    end
    
end


function ShowRegisteredTogether( input_images, ref_round )
%ShowRegistered is used to plot registration results 
%   -----IO-----
%   input_images: cell with images 
%   ref_round: reference round

    ref = input_images{ref_round};
    ref_3d = max(ref, [], 4);
    ref_2d = max(ref_3d, [], 3);
    Nround = numel(input_images);
    
    figure
    for r=2:Nround
        curr_round = input_images{r};
        curr_round_3d = max(curr_round, [], 4);
        curr_round_2d = max(curr_round_3d, [], 3);
        
        curr_ssim_3d = ssim(curr_round_3d, ref_3d);     
        curr_ssim_2d = ssim(curr_round_2d, ref_2d);
        
        subplot(1, Nround-1, r-1), imshowpair(ref_2d, curr_round_2d), title(sprintf("Round %d vs. Round %d", ref_round, r)), xlabel(sprintf("3D SSIM: %.3f\n2D SSIM: %.3f", curr_ssim_3d, curr_ssim_2d))
    end

end


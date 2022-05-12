function ShowRegistered( input_images, ref_round, prefix, saveAS )
%ShowRegistered is used to plot registration results 
%   -----IO-----
%   input_images: cell with images 
%   ref_round: reference round
    if nargin < 3
        prefix = '';
        saveAS = false;
    end
    
    ref = input_images{ref_round};
    ref_3d = max(ref, [], 4);
    ref_2d = max(ref_3d, [], 3);
    Nround = numel(input_images);
    
    if saveAS
        output_dir = sprintf("%s/", prefix);

        % curr_dir = fullfile(dirName, sprintf('round%d',i));
        if ~exist(output_dir, 'dir')
           mkdir(output_dir);
        end
    end
    
    rounds = 1:Nround;
    rounds = rounds(rounds ~= ref_round);
            
    for r=rounds
        curr_round = input_images{r};
        curr_round_3d = max(curr_round, [], 4);
        curr_round_2d = max(curr_round_3d, [], 3);
        
        curr_ssim_3d = ssim(curr_round_3d, ref_3d);     
        curr_ssim_2d = ssim(curr_round_2d, ref_2d);
        
        figure
        imshowpair(ref_2d, curr_round_2d), title(sprintf("Round %d vs. Round %d", ref_round, r)), xlabel(sprintf("3D SSIM: %.3f\n2D SSIM: %.3f", curr_ssim_3d, curr_ssim_2d))
        
        if saveAS
            fig_name = sprintf("%sim_%d%d.png", output_dir, ref_round, r);
            saveas(gcf, fig_name)
        end
        
        drawnow;
    end

end


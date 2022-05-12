function cell_imgs = RegisterCellImage( obj, input_path, cell_folder, sub_dir, amplicon_channel, useGPU )


fprintf("====Register cell images====\n");
fprintf(obj.log, "====Register cell images====\n");

tic;
cell_path = fullfile(input_path, cell_folder, sub_dir);
cell_files = dir(fullfile(cell_path, '*.tif'));
nfiles = numel(cell_files);
cell_imgs = cell(nfiles, 1);

% Load all channels
for c=1:nfiles 
    curr_path = strcat(cell_files(c).folder, '/', cell_files(c).name);

    curr_img = new_LoadMultipageTiff(curr_path, 'uint8', 'uint8', false);
    
    cell_imgs{c} = curr_img;
end

if useGPU
    ref_round = gpuArray(obj.rawImages(:,:,:,:,1));
else
    ref_round = obj.rawImages(:,:,:,:,1);
end
fix = max(ref_round, [], 4);

params = DFTRegister3D(fix, cell_imgs{amplicon_channel}, false);
fprintf(sprintf('Shifted by %s\n', num2str(params.shifts)));
fprintf(obj.log, sprintf('Shifted by %s\n', num2str(params.shifts)));

for c=1:nfiles
    if useGPU
        curr_reg = DFTApply3D(gpuArray(cell_imgs{c}), params, false);
        cell_imgs{c} = uint8(gather(curr_reg));
    else
        curr_reg = DFTApply3D(cell_imgs{c}, params, false);
        cell_imgs{c} = uint8(curr_reg);
    end
end

fprintf(sprintf('[time = %.2f s]\n', toc));
reset(gpuDevice)

end


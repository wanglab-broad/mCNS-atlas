%% set batch id
% label_file=1;
label_file=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];

% label_file=[13];
all_file_name=dir('J:\ClusterMap2\ccf\AP_histology-master_sagittal\slices\slices\*.tif');
% label_id=num2str(label_file);
filename={};
for lengi=1:length(label_file)
    filename{lengi}=all_file_name(label_file(lengi)).name;
    filename{lengi}= erase(filename{lengi},'.tif');
end

%% 1) Load CCF and set paths for slide and slice images

% Load CCF atlas
allen_atlas_path = 'J:\ClusterMap2\ccf\AP_histology-master_sagittal';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% Set paths for histology images and directory to save slice/alignment
im_path = 'J:\ClusterMap2\ccf\AP_histology-master_sagittal\slices';
slice_path = [im_path filesep 'slices'];

%% 2) 
% read our data in slices
% im_path=slice_path;
slice_im_path=slice_path;

% Initialize guidata
gui_data = struct;
gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;

% Load in slice images
gui_data.slice_im_path = slice_im_path;
for lengi=1:length(label_file)
    slice_im_dir = dir([im_path filesep 'slices\' filename{lengi} '.tif']);
    slice_im_fn{lengi} = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
        {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
end

%% pad and rotate
% Pad all slices centrally to the largest slice and make matrix
slice_im = cell(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
   slice_im{curr_slice} = imread(string(slice_im_fn{curr_slice}));  
end
slice_size_max = [5000,5000,3];%max(cell2mat(cellfun(@size,slice_im,'uni',false)),[],1);
slice_im_pad = ...
    cell2mat(cellfun(@(x) x(1:slice_size_max(1),1:slice_size_max(2),:), ...
    reshape(cellfun(@(im) padarray(im, ...
    [ceil((slice_size_max(1) - size(im,1))./2), ...
    ceil((slice_size_max(2) - size(im,2))./2)],0,'both'), ...
    slice_im,'uni',false),1,1,1,[]),'uni',false));

%% rotate for each dataset
im_aligned = zeros([5000,5000,3,length(label_file)],class(slice_im_pad));
% num2str(label_file(lengi))
for lengi=1:length(label_file)
    align_axis=load(['align_axis',filename{lengi},'.mat']);
    align_axis=align_axis.align_axis;
    align_angle = squeeze(atan2d(diff(align_axis(:,1,:),[],1),diff(align_axis(:,2,:),[],1)));
    align_center = squeeze(nanmean(align_axis,1));
    % Set target angle as the nearest multiple of 90
    target_angle = round(nanmean(align_angle)/90)*90;
    % Set target position as the average center of the reference lines
    if size(align_center,1)>1
        target_position = nanmean(align_center,2);
    else
        target_position = align_center;
    end
    
    curr_im =1;% 1:length(slice_im)

    angle_diff = target_angle - align_angle(curr_im);
    if size(align_center,1)>1
        x_diff = target_position(2) - align_center(2,curr_im);
        y_diff = target_position(1) - align_center(1,curr_im);
    else
        x_diff=0;
        y_diff=0;
    end
    im_aligned(:,:,:,lengi) = ...
        imrotate(imtranslate(slice_im_pad(:,:,:,lengi), ...
        [x_diff,y_diff]),angle_diff,'bilinear','crop');    
    
end

% replace with new cell type map
gui_data.slice_im = cell(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
    gui_data.slice_im{curr_slice} = im_aligned(:,:,:,curr_slice);
end

%% Load corresponding CCF slices
histology_ccf_all=struct([]);
for lengi=1:length(label_file)
    
    ccf_slice_fn_part = [slice_im_path filesep 'histology_ccf',filename{lengi},'.mat'];
    load(ccf_slice_fn_part);
    histology_ccf_all= cat(1,histology_ccf_all,histology_ccf); 
    
end
histology_ccf=histology_ccf_all;
gui_data.histology_ccf = histology_ccf;

% Load histology/CCF alignment
atlas2histology_tform_all=cell(1,length(label_file));
for lengi=1:length(label_file)
    
    ccf_alignment_fn = [slice_im_path filesep 'atlas2histology_tform',filename{lengi},'.mat'];
    load(ccf_alignment_fn);
    atlas2histology_tform_all(1,lengi) = atlas2histology_tform; 
end

atlas2histology_tform=atlas2histology_tform_all;
gui_data.histology_ccf_alignment = atlas2histology_tform;

% Warp area labels by histology alignment
gui_data.histology_aligned_av_slices = cell(length(gui_data.slice_im),1);
for curr_slice = 1:length(gui_data.slice_im)
    curr_av_slice = gui_data.histology_ccf(curr_slice).av_slices;
    curr_av_slice(isnan(curr_av_slice)) = 1;

    %our data
    curr_slice_im = gui_data.slice_im{curr_slice};

    % transform
    tform = affine2d;
    tform.T = gui_data.histology_ccf_alignment{curr_slice};   
    tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);

    gui_data.histology_aligned_av_slices{curr_slice} = ...
        imwarp(curr_av_slice,tform,'nearest','OutputView',tform_size);
end



% ccf_reftest=gui_data.histology_ccf(1).av_slices;
% [Gmag, ~] = imgradient(ccf_reftest,'prewitt');
%     Gmag=uint8(Gmag);
%     point_at_boundary=find(Gmag>0);
%     sz=size(Gmag);
%     [row,col] = ind2sub(sz,point_at_boundary);
%     figure;imshow(ccf_reftest);
    


%% Create figure, set button functions
for ij=1:length(gui_data.slice_im)
    
    gui_data.curr_slice = ij;
    gui_fig = figure('KeyPressFcn',@keypress,'WindowButtonMotionFcn',@mousehover);
    % Set up axis for histology image
    gui_data.histology_ax = axes('YDir','reverse'); 
    hold on; colormap(gray); axis image off;
    gui_data.histology_im_h = image(gui_data.slice_im{gui_data.curr_slice}, ...
        'Parent',gui_data.histology_ax);

    % Create title to write area in
    gui_data.histology_ax_title = title(gui_data.histology_ax,'','FontSize',10);

    % Set up histology-aligned atlas overlay
    % (and make it invisible to mouse clicks)
    histology_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.slice_im{gui_data.curr_slice},1),size(gui_data.slice_im{gui_data.curr_slice},2));
    %
       % 
    gui_data.histology_aligned_atlas_boundaries = ...
        imagesc(histology_aligned_atlas_boundaries_init,'Parent',gui_data.histology_ax, ...
        'AlphaData',histology_aligned_atlas_boundaries_init,'PickableParts','none');

    % plot boundary
    ourimage=gui_data.slice_im{gui_data.curr_slice};
    ccf_ref=gui_data.histology_aligned_av_slices{gui_data.curr_slice};
    [Gmag, ~] = imgradient(ccf_ref,'prewitt');
    Gmag=uint8(Gmag);
    point_at_boundary=find(Gmag>0);
    sz=size(Gmag);
    [row,col] = ind2sub(sz,point_at_boundary);
    hold on;
    plot(col,row,'r.','markersize',3);

    % Upload gui data
    guidata(gui_fig,gui_data);

    % Update the slice
    update_slice(gui_fig);
 end

%% final) Display histology within 3D CCF
slice_im_path=slice_path;
gui_data.tv = tv;
gui_data.av = av;

% Warp histology to CCF
gui_data.atlas_aligned_histology = cell(length(gui_data.slice_im),1);
for curr_slice = 1:length(gui_data.slice_im)
    curr_av_slice = gui_data.histology_ccf(curr_slice).av_slices;
    curr_av_slice(isnan(curr_av_slice)) = 1;
    curr_slice_im = gui_data.slice_im{curr_slice};
    
    tform = affine2d;
    tform.T = gui_data.histology_ccf_alignment{curr_slice};
    % (transform is CCF -> histology, invert for other direction)
    tform = invert(tform);

    tform_size = imref2d([size(gui_data.histology_ccf(curr_slice).av_slices,1), ...
        size(gui_data.histology_ccf(curr_slice).av_slices,2)]);
    
    gui_data.atlas_aligned_histology{curr_slice} = ...
        imwarp(curr_slice_im,tform,'nearest','OutputView',tform_size);
end

% Create figure
gui_fig = figure;

% Set up 3D plot for volume viewing
axes_atlas = axes;
[~, brain_outline] = plotBrainGrid([],axes_atlas);
set(axes_atlas,'ZDir','reverse');
hold(axes_atlas,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])


% Turn on rotation by default
h = rotate3d(axes_atlas);
h.Enable = 'on';

% Draw all aligned slices
histology_surf = gobjects(length(gui_data.slice_im),1);
% Rcolor=rand(3,3,'double');
for curr_slice = 1:length(gui_data.slice_im)
    
    test= gui_data.atlas_aligned_histology{curr_slice};        
    color_valu=zeros(800*1140,3);
    t=1;
    for i=1:1140
        for j=1:800
            color_valu(t,:)=test(j,i,:);
            t=t+1;
        end
    end

    cell_ind=find(sum(color_valu,2)>100 & max(color_valu,[],2)>30 & sum(color_valu,2)<735);

    X=gui_data.histology_ccf(curr_slice).plane_ap(:);
    Y=gui_data.histology_ccf(curr_slice).plane_ml(:);
    Z=gui_data.histology_ccf(curr_slice).plane_dv(:);
    % C =  gui_data.atlas_aligned_histology{curr_slice};
    XX=X(cell_ind);
    YY=Y(cell_ind);
    ZZ=Z(cell_ind);
    CC=color_valu(cell_ind,:)/255;
    hold on;
    scatter3(XX,YY,ZZ,3,CC,'filled');
end


























function keypress_rotate(gui_fig,eventdata)

shift_on = any(strcmp(eventdata.Modifier,'shift'));

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % 1/2: switch slice
    % Shift + 1/2: move slice in stack

    case '1'
        if ~shift_on
            gui_data.curr_slice = max(gui_data.curr_slice - 1,1);
            set(gui_data.histology_im_h,'CData',gui_data.im_aligned(:,:,:,gui_data.curr_slice))
            guidata(gui_fig,gui_data);
        elseif shift_on && gui_data.curr_slice ~= 1
            slice_flip = [gui_data.curr_slice-1,gui_data.curr_slice];
            gui_data.im_aligned(:,:,:,slice_flip) = flip(gui_data.im_aligned(:,:,:,slice_flip),4);
            gui_data.curr_slice = slice_flip(1);
            guidata(gui_fig,gui_data);
        end
        
    case '2'
        if ~shift_on
            gui_data.curr_slice = ...
                min(gui_data.curr_slice + 1,size(gui_data.im_aligned,4));
            set(gui_data.histology_im_h,'CData',gui_data.im_aligned(:,:,:,gui_data.curr_slice))
            guidata(gui_fig,gui_data);
        elseif shift_on && gui_data.curr_slice ~= size(gui_data.im_aligned,4)
            slice_flip = [gui_data.curr_slice,gui_data.curr_slice+1];
            gui_data.im_aligned(:,:,:,slice_flip) = flip(gui_data.im_aligned(:,:,:,slice_flip),4);
            gui_data.curr_slice = slice_flip(2);
            guidata(gui_fig,gui_data);
        end
        
    % Arrow keys: flip slice
    case {'leftarrow','rightarrow'}
        gui_data.im_aligned(:,:,:,gui_data.curr_slice) = ...
            fliplr(gui_data.im_aligned(:,:,:,gui_data.curr_slice));
        set(gui_data.histology_im_h,'CData',gui_data.im_aligned(:,:,:,gui_data.curr_slice))
        guidata(gui_fig,gui_data);
               
    case {'uparrow','downarrow'}
        gui_data.im_aligned(:,:,:,gui_data.curr_slice) = ...
            flipud(gui_data.im_aligned(:,:,:,gui_data.curr_slice));
        set(gui_data.histology_im_h,'CData',gui_data.im_aligned(:,:,:,gui_data.curr_slice))
        guidata(gui_fig,gui_data);
        
    % Escape: save and close
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?','Confirm exit',opts);
        if strcmp(user_confirm,'Yes')
            % Overwrite old images with new ones
            for curr_im = 1:size(gui_data.im_aligned,4)
                imwrite(gui_data.im_aligned(:,:,:,curr_im),gui_data.slice_fn{curr_im},'tif');
            end
            disp(['Saved padded/centered/rotated slices']);
            close(gui_fig)
        end
        
end

end

function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % Arrow keys: rotate atlas slice
    case 'leftarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [1,0]);
        update_atlas_slice(gui_fig)
    case 'rightarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [-1,0]);
        update_atlas_slice(gui_fig)
    case 'uparrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [0,-1]);
        update_atlas_slice(gui_fig)
    case 'downarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [0,1]);
        update_atlas_slice(gui_fig)
    
    % 1/2 keys: cycle through histology slices
    % (if there's a saved plane point, move atlas to that position)
    case '1'
        gui_data.curr_histology_slice = max(gui_data.curr_histology_slice - 1,1);            
        guidata(gui_fig,gui_data);
        update_histology_slice(gui_fig);
        
    case '2'
        gui_data.curr_histology_slice = ...
            min(gui_data.curr_histology_slice + 1,length(gui_data.slice_im));
        guidata(gui_fig,gui_data);
        update_histology_slice(gui_fig);
        
    % Enter: save slice coordinates
    case 'return'        
        % Store camera vector and point
        % (Note: only one camera vector used for all slices, overwrites)
        gui_data.slice_vector = get_camera_vector(gui_data);
        gui_data.slice_points(gui_data.curr_histology_slice,:) = ...
            gui_data.atlas_slice_point;
        guidata(gui_fig,gui_data);
                
        update_histology_slice(gui_fig);
        title(gui_data.histology_ax,'New saved atlas position');
        
    % Escape: save and exit
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?','Confirm exit',opts);
        if strcmp(user_confirm,'Yes')
            
            % Check that a CCF slice point exists for each histology slice
            if any(isnan(gui_data.slice_points(:)))
                createmode = struct;
                createmode.Interpreter = 'tex';
                createmode.WindowStyle = 'modal';
                msgbox('\fontsize{12} Some histology slice(s) not assigned CCF slice', ...
                    'Not saving','error',createmode);
                return
            end
            
            % Go through each slice, pull full-resolution atlas slice and
            % corrsponding coordinates       
            histology_ccf_init = cell(length(gui_data.slice_im),1);
            histology_ccf = struct( ...
                'tv_slices',histology_ccf_init, ...
                'av_slices',histology_ccf_init, ...
                'plane_ap',histology_ccf_init, ...
                'plane_ml',histology_ccf_init, ...
                'plane_dv',histology_ccf_init);
            
            h = waitbar(0,'Saving atlas slices...');
            for curr_slice = 1:length(gui_data.slice_im)
                gui_data.atlas_slice_point = gui_data.slice_points(curr_slice,:);
                [histology_ccf(curr_slice).tv_slices, ...
                    histology_ccf(curr_slice).av_slices, ...
                    histology_ccf(curr_slice).plane_ap, ...
                    histology_ccf(curr_slice).plane_ml, ...
                    histology_ccf(curr_slice).plane_dv] = ...
                    grab_atlas_slice(gui_data,1);
                waitbar(curr_slice/length(gui_data.slice_im),h, ...
                    ['Saving atlas slices (' num2str(curr_slice) '/' num2str(length(gui_data.slice_im)) ')...']);
            end                     
            close(h);
            
            save_fn = [gui_data.slice_im_path filesep 'histology_ccf.mat'];
            save(save_fn,'histology_ccf','-v7.3');
            close(gui_fig);            
        end
end

end

function update_histology_slice(gui_fig)
% Draw histology slice (and move atlas if saved position)

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h,'CData',gui_data.slice_im{gui_data.curr_histology_slice})

% If there's a saved atlas position, move atlas to there
if all(~isnan(gui_data.slice_points(gui_data.curr_histology_slice,:)))
    gui_data.atlas_slice_point = ...
        gui_data.slice_points(gui_data.curr_histology_slice,:);
    title(gui_data.histology_ax,'Saved atlas position')
    guidata(gui_fig,gui_data);
    update_atlas_slice(gui_fig);
else
    title(gui_data.histology_ax,'No saved atlas position')
end

% Upload gui data
guidata(gui_fig, gui_data);

end

function cam_vector = get_camera_vector(gui_data)
% Get the camera viewing vector to define atlas slice plane

% Grab current camera angle

% (Old way: more confusing, easily messed up by axes directions)
% [cam_az,cam_el] = view(gui_data.atlas_ax);
% 
% % Camera azimuth is 90 degrees offset from spherical standard (?!)
% cam_az_sphere = cam_az - 90;
% % Camera elevation is reversed (because of CCF orientation)
% cam_el_sphere = -cam_el;
% 
% [cam_vector_x,cam_vector_y,cam_vector_z] = ...
%     sph2cart(deg2rad(cam_az_sphere),deg2rad(cam_el_sphere),1);
% cam_vector = [cam_vector_x,cam_vector_y,cam_vector_z];

% (New way: just a normalized line from the camera to the center)
curr_campos = campos(gui_data.atlas_ax);
curr_camtarget = camtarget(gui_data.atlas_ax);
cam_vector = (curr_camtarget - curr_campos)./norm(curr_camtarget - curr_campos);

end

function scroll_atlas_slice(gui_fig,eventdata)
% Move point to draw atlas slice perpendicular to the camera

% Get guidata
gui_data = guidata(gui_fig);

% Move slice point along camera -> center axis
cam_vector = get_camera_vector(gui_data);

% Move slice point
gui_data.atlas_slice_point = gui_data.atlas_slice_point + ...
    eventdata.VerticalScrollCount*cam_vector;

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_atlas_slice(gui_fig)

end

function update_atlas_slice(gui_fig)
% Draw atlas slice through plane perpendicular to camera through set point

% Get guidata
gui_data = guidata(gui_fig);

% Get slice (larger spacing for faster pulling)
[tv_slice,av_slice,plane_ap,plane_ml,plane_dv] = grab_atlas_slice(gui_data,3);

% Update the slice display
set(gui_data.atlas_slice_plot,'XData',plane_ap,'YData',plane_ml,'ZData',plane_dv,'CData',tv_slice);

% Upload gui_data
guidata(gui_fig, gui_data);

end

function [tv_slice,av_slice,plane_ap,plane_ml,plane_dv] = grab_atlas_slice(gui_data,slice_px_space)
% Grab anatomical and labelled atlas within slice

% Get plane normal to the camera -> center axis, grab voxels on plane
cam_vector = get_camera_vector(gui_data);
plane_offset = -(cam_vector*gui_data.atlas_slice_point');

% Define a plane of points to index
% (the plane grid is defined based on the which cardinal plan is most
% orthogonal to the plotted plane. this is janky but it works)

[~,cam_plane] = max(abs(cam_vector./norm(cam_vector)));

switch cam_plane
    
    % Note: ML and DV directions are flipped to match 2D histology and 3D
    % atlas axes, so make ML and DV coordinates go backwards for true CCF
    % coordinates
    
    case 1
        [plane_ml,plane_dv] = ...
            meshgrid(1:slice_px_space:size(gui_data.tv,3), ...
            1:slice_px_space:size(gui_data.tv,2));
        plane_ap = ...
            (cam_vector(2)*plane_ml+cam_vector(3)*plane_dv + plane_offset)/ ...
            -cam_vector(1);
        
    case 2
        [plane_ap,plane_dv] = ...
            meshgrid(1:slice_px_space:size(gui_data.tv,1), ...
            1:slice_px_space:size(gui_data.tv,2));
        plane_ml = ...
            (cam_vector(1)*plane_ap+cam_vector(3)*plane_dv + plane_offset)/ ...
            -cam_vector(2);
        
    case 3
        [plane_ap,plane_ml] = ...
            meshgrid(size(gui_data.tv,3):-slice_px_space:1, ...
            1:slice_px_space:size(gui_data.tv,3));
        plane_dv = ...
            (cam_vector(1)*plane_ap+cam_vector(2)*plane_ml + plane_offset)/ ...
            -cam_vector(3);
        
end

% Get the coordiates on the plane
ap_idx = round(plane_ap);
ml_idx = round(plane_ml);
dv_idx = round(plane_dv);

% Find plane coordinates in bounds with the volume
% (CCF coordinates: [AP,DV,ML])
use_ap = ap_idx > 0 & ap_idx < size(gui_data.tv,1);
use_dv = dv_idx > 0 & dv_idx < size(gui_data.tv,2);
use_ml = ml_idx > 0 & ml_idx < size(gui_data.tv,3);
use_idx = use_ap & use_ml & use_dv;

curr_slice_idx = sub2ind(size(gui_data.tv),ap_idx(use_idx),dv_idx(use_idx),ml_idx(use_idx));

% Find plane coordinates that contain brain
curr_slice_isbrain = false(size(use_idx));
curr_slice_isbrain(use_idx) = gui_data.av(curr_slice_idx) > 0;

% Index coordinates in bounds + with brain
grab_pix_idx = sub2ind(size(gui_data.tv),ap_idx(curr_slice_isbrain),dv_idx(curr_slice_isbrain),ml_idx(curr_slice_isbrain));

% Grab pixels from (selected) volume
tv_slice = nan(size(use_idx));
tv_slice(curr_slice_isbrain) = gui_data.tv(grab_pix_idx);

av_slice = nan(size(use_idx));
av_slice(curr_slice_isbrain) = gui_data.av(grab_pix_idx);

end

function mousehover(gui_fig,eventdata)
% Display area of atlas on mouse hover

% Get guidata
gui_data = guidata(gui_fig);

% Get mouse position
mouse_position = get(gui_data.histology_ax,'CurrentPoint');
mouse_x = round(mouse_position(1,1));
mouse_y = round(mouse_position(1,2));

curr_av_slice_warp = gui_data.histology_aligned_av_slices{gui_data.curr_slice};

% Don't use if mouse out of bounds
if ~ismember(mouse_x,1:size(curr_av_slice_warp,2)) || ...
        ~ismember(mouse_y,1:size(curr_av_slice_warp,1))
    return
end
    
curr_av = curr_av_slice_warp(mouse_y,mouse_x);

% Don't use if AV = 0
if curr_av == 0
    return
end

% Grab area name and set title
curr_area_name = gui_data.st.safe_name(curr_av);
set(gui_data.histology_ax_title,'String',curr_area_name);

end

function align_ccf_to_histology(gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

curr_av_slice_warp = gui_data.histology_aligned_av_slices{gui_data.curr_slice};
av_warp_boundaries = round(conv2(curr_av_slice_warp,ones(3)./9,'same')) ~= curr_av_slice_warp;

set(gui_data.histology_aligned_atlas_boundaries, ...
    'CData',av_warp_boundaries, ...
    'AlphaData',av_warp_boundaries*0.3);

% Upload gui data
guidata(gui_fig, gui_data);

end

function update_slice(gui_fig)
% Draw histology and CCF slice

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h,'CData',gui_data.slice_im{gui_data.curr_slice})

% Upload gui data
guidata(gui_fig, gui_data);

% Update atlas boundaries
align_ccf_to_histology(gui_fig)

end

function keypress_manual(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % 1/2: move slice
    case '1'
        gui_data.curr_slice = max(gui_data.curr_slice - 1,1);
        guidata(gui_fig,gui_data);
        update_slice_align(gui_fig);
        
    case '2'
        gui_data.curr_slice = ...
            min(gui_data.curr_slice + 1,length(gui_data.slice_im));
        guidata(gui_fig,gui_data);
        update_slice_align(gui_fig);
        
    % O: toggle overlay visibility
    case 'space'
        curr_visibility = ...
            get(gui_data.histology_aligned_atlas_boundaries,'Visible');
        set(gui_data.histology_aligned_atlas_boundaries,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility)))
        
    % C: clear current points
    case 'c'
        gui_data.histology_control_points{gui_data.curr_slice} = zeros(0,2);
        gui_data.atlas_control_points{gui_data.curr_slice} = zeros(0,2);
        
        guidata(gui_fig,gui_data);
        update_slice_align(gui_fig);
        
    % S: save
    case 's'
        atlas2histology_tform = ...
            gui_data.histology_ccf_manual_alignment;
        save_fn = [gui_data.slice_im_path filesep 'atlas2histology_tform.mat'];
        save(save_fn,'atlas2histology_tform');
        disp(['Saved ' save_fn]);
        
    % Escape: save and exit
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?','Confirm exit',opts);
        if strcmp(user_confirm,'Yes')            
            atlas2histology_tform = ...
                gui_data.histology_ccf_manual_alignment;
            save_fn = [gui_data.slice_im_path filesep 'atlas2histology_tform.mat'];
            save(save_fn,'atlas2histology_tform');
            disp(['Saved ' save_fn]);
            close(gui_fig);            
        end
        
end

end


function mouseclick_histology_manual(gui_fig,eventdata)
% Draw new point for alignment

% Get guidata
gui_data = guidata(gui_fig);

% Add clicked location to control points
gui_data.histology_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.histology_control_points{gui_data.curr_slice}, ...
    eventdata.IntersectionPoint(1:2));

set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,1), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,2));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) || ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 3 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 3)
    align_ccf_to_histology_manual(gui_fig)
end

end


function mouseclick_atlas_manual(gui_fig,eventdata)
% Draw new point for alignment

% Get guidata
gui_data = guidata(gui_fig);

% Add clicked location to control points
gui_data.atlas_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.atlas_control_points{gui_data.curr_slice}, ...
    eventdata.IntersectionPoint(1:2));

set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,1), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,2));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) || ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 3 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 3)
    align_ccf_to_histology_manual(gui_fig)
end

end


function align_ccf_to_histology_manual(gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) && ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) >= 3 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) >= 3)    
    % If same number of >= 3 control points, use control point alignment
    tform = fitgeotrans(gui_data.atlas_control_points{gui_data.curr_slice}, ...
        gui_data.histology_control_points{gui_data.curr_slice},'affine');
    title(gui_data.histology_ax,'New alignment');
    
      
elseif size(gui_data.histology_control_points{gui_data.curr_slice},1) >= 1 ||  ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) >= 1
    % If less than 3 or nonmatching points, use auto but don't draw
    title(gui_data.histology_ax,'New alignment');
    
    % Upload gui data
    guidata(gui_fig, gui_data);
    return
    
else
    % If no points, use automated outline
    if isfield(gui_data,'histology_ccf_auto_alignment')
        tform = affine2d;
        tform.T = gui_data.histology_ccf_auto_alignment{gui_data.curr_slice};
        title(gui_data.histology_ax,'Previous alignment');
    end
end

curr_av_slice = gui_data.histology_ccf(gui_data.curr_slice).av_slices;
curr_av_slice(isnan(curr_av_slice)) = 1;
curr_slice_im = gui_data.slice_im{gui_data.curr_slice};

tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
curr_av_slice_warp = imwarp(curr_av_slice, tform, 'OutputView',tform_size);

av_warp_boundaries = round(conv2(curr_av_slice_warp,ones(3)./9,'same')) ~= curr_av_slice_warp;

set(gui_data.histology_aligned_atlas_boundaries, ...
    'CData',av_warp_boundaries, ...
    'AlphaData',av_warp_boundaries*0.3);

% Update transform matrix
gui_data.histology_ccf_manual_alignment{gui_data.curr_slice} = tform.T;

% Upload gui data
guidata(gui_fig, gui_data);

end


function update_slice_align(gui_fig)
% Draw histology and CCF slice

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h,'CData',gui_data.slice_im{gui_data.curr_slice})
set(gui_data.atlas_im_h,'CData',gui_data.histology_ccf(gui_data.curr_slice).tv_slices);

% Plot control points for slice
set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,1), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,2));
set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,1), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,2));

% Reset histology-aligned atlas boundaries if not
histology_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.slice_im{1},1),size(gui_data.slice_im{1},2));
set(gui_data.histology_aligned_atlas_boundaries, ...
    'CData',histology_aligned_atlas_boundaries_init, ...
    'AlphaData',histology_aligned_atlas_boundaries_init);

% Upload gui data
guidata(gui_fig, gui_data);

% Update atlas boundaries
align_ccf_to_histology_manual(gui_fig)

end







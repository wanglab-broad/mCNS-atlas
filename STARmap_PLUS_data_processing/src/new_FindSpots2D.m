function centroids = new_FindSpots2D( imgstack, thresholdRel, kdims, ksize )
%new_FindSpots2D Find spots in single plane

    % cc_props = {'Volume'; 'Centroid'; 'EquivDiameter'; 'VoxelIdxList'};
    cc_props = {'Centroid'};


    %imgstack = single(imgstack);
    thresholdRel = thresholdRel * max(imgstack(:));

    Nslice = size(imgstack, 3);

    h = fspecial('log', kdims, ksize);

    filteredImg = zeros(size(imgstack));
    for i=1:Nslice
        % apply LoG filtering to each plane
        filteredImg(:,:,i) = imfilter(single(imgstack(:,:,i)), h);
    end

    rmax = boolean(zeros(size(imgstack)));

    for i=1:Nslice
        rmax(:,:,i) = boolean(gather(imregionalmin(filteredImg(:,:,i))));
    end

    for i=1:Nslice
        rmax(:,:,i) = rmax(:,:,i) & (imgstack(:,:,i) > thresholdRel);
    end

    CC = bwconncomp(rmax);
    CC_stats = regionprops3(CC, cc_props);

    centroidProps = CC_stats.Centroid;
    centroids = int16(centroidProps);

end


function plot_centroids(input_centroid, input_img, msize, color, input_title)
% Plot centroids on a gray-scale and a bw image 
if nargin < 4
    color = 'red';
end

if nargin < 5
    input_title = '';
end

    figure
    imshow(input_img)
    hold on
    plot(input_centroid(:,1), input_centroid(:,2), '.', "Color", color, "MarkerSize", msize)
    title(input_title)
    hold off
    
end
function output_dist = GetDist( input_points, dist )
%GETDIST

    dist_mat = squareform(pdist(input_points, 'euclidean'));
    upper_dist = triu(dist_mat, 1);
    close_points = upper_dist > 0 & upper_dist < dist;
    num_close_points = sum(close_points, 'all');
    
    fprintf(sprintf("Number of paired points with distance less than %d: %d\n", dist, num_close_points));


end


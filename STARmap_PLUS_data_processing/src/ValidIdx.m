function output_idx = ValidIdx( input_idx, size )
%


    valid_x = logical([]);
    valid_y = logical([]);
    
    for i=1:numel(input_idx)
        curr_idx = input_idx{i};
        curr_x = curr_idx(1,2) - curr_idx(1,1) + 1;
        curr_y = curr_idx(2,2) - curr_idx(2,1) + 1;

        curr_x = curr_x >= size;
        curr_y = curr_y >= size;

        valid_x = [valid_x; curr_x];
        valid_y = [valid_y; curr_y];
    end

    output_idx = input_idx(valid_y & valid_x);



end


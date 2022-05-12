function output = BlockTest( input_stack, block_size, overlap_percent )
%BlockTest

    [x, y, ~, ~, ~] = size(input_stack);
    block_x = block_size(1);
    block_y = block_size(2);
    
    offset_x = block_x * overlap_percent;
    offset_y = block_y * overlap_percent;
%     block_z = block_size(3);

    % block limits

    % ... starting points
    r0 = 1:block_x-offset_x:x;
    c0 = 1:block_y-offset_y:y;
%     s0 = 1:block_z:z;

    % ... end points
    rx = r0 + block_x - 1;
    cx = c0 + block_y - 1;
%     sx = s0 + block_z - 1;
    rx = min(rx, x);
    cx = min(cx, y);
%     sx = min(sx, z);

    NR = length(r0);
    NC = length(c0);
%     NS = length(s0);
    
%     numblocks = NR * NC * NS;
    
    output = {};
    
    for i=1:NR
        for j=1:NC
            output{end+1} = [r0(i) rx(i); c0(j) cx(j)];
            
%             for k=1:NS
%                 curr_block = input_stack(...
%                     r0(i):rx(i),...
%                     c0(j):cx(j),...
%                     s0(k):sx(k),...
%                     :,:);

%                 output{end+1} = [r0; rx; c0; cx];
               
                
                % test 
%                 curr_img = curr_block(:,:,1,1,1);              
%                 output{end+1} = curr_img;
%                 figure
%                 imshow(curr_img)
%                 drawnow;
                 
%             end
        end
    end
    
    output = cellfun(@int16, output, 'UniformOutput', false);
end


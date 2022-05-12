function [img_stack] = hist_match_of_channels(img_stack)
%hist_match_of_channels is used to equalize or normalize intensity profile
%among all channels in each round
%   -----IO-----
%   img_stack: 5D array XxYxZxChxRd
%   img_stack: Histogram matched img_stack
    Nround = size(img_stack, 5);
    Nchannel = size(img_stack, 4);
    
    for r=1:Nround 
        fprintf('Equalizing intensity histogram between each channel in round %d\n', r);
        currStack = uint8(img_stack(:,:,:,:,r));
        signal_min = 0;
        signal = 0;
        for c=1:Nchannel
            curr_signal = currStack(:,:,:,c) ~= 0;
            curr_signal = sum(curr_signal, 'all');
            if signal_min == 0
                signal_min = c;
                signal = curr_signal;
            elseif curr_signal < signal
                signal_min = c;
                signal = curr_signal;
            end
        end    
        signal_min
        for c=1:Nchannel
            currStack(:,:,:,c) = uint8(imhistmatchn(uint8(currStack(:,:,:,c)), uint8(currStack(:,:,:,signal_min)), 200 ));
            %currStack(:,:,:,c) = gather(histeq(gpuArray(currStack(:,:,:,c))));
        end  
        
        img_stack(:,:,:,:,r) = currStack;
    end


end


function output_img = ArbitraryThresholding( input_img, varargin )
%ArbitraryThresholding

    % Input parser
    p = inputParser;
    % Defaults
    defaultfactor = 0.2;
    defaultDoSM = true;
    defaultSigma = 1;
    defaultFilter = [5 5 3];
    
    addRequired(p,'input_img');
    addOptional(p,'factor',defaultfactor);
    addOptional(p,'DoSmoothing',defaultDoSM);
    addParameter(p,'Sigma',defaultSigma);
    addParameter(p,'FilterSize',defaultFilter);
    
    parse(p, input_img, varargin{:});
    
    Nround = numel(input_img);
    output_img = cell(Nround, 1);
    
    %thres_stack = uint8(zeros(size(input_stack)));
    %bw_stack = zeros(size(input_stack));
    
    for r = 1:Nround
        fprintf(sprintf("Thresholding Round %d...\n", r));
        curr_round = input_img{r};
        Nchannel = size(curr_round, 4);
        thres_stack = uint8(zeros(size(curr_round)));
        for c = 1:Nchannel
            curr_channel = curr_round(:,:,:,c);
            curr_max = max(curr_channel, [], 'all');
            curr_boundary = curr_max * p.Results.factor;
            curr_bw = curr_channel > curr_boundary;
            curr_thres = curr_channel;
            curr_thres(~curr_bw) = 0;
            
            %bw_stack(:,:,:,c,r) = curr_bw;
            
            if p.Results.DoSmoothing
                fprintf(sprintf("Smoothing Channel %d...\n", c));
                thres_stack(:,:,:,c) = imgaussfilt3(uint8(curr_thres), p.Results.Sigma, 'FilterSize', p.Results.FilterSize);
            else
                thres_stack(:,:,:,c) = uint8(curr_thres);
            end
            
        end
        
        output_img{r} = thres_stack;
    end
    
    %bw_stack = logical(bw_stack);
    
end


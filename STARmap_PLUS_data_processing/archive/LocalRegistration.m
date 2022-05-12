function output_img = LocalRegistration( input_img, varargin )
%LocalRegistration is used to do local (non-rigid) registration for images
%in all rounds
%   -----IO-----
%   input_img: input cell with image stacks for all rounds 

    % Input parser
    p = inputParser;
    % Defaults
    defaultRef = 1;
    defaultMethod = "max";
    defaultIter = 60;
    defaultAFS = 1;

    addRequired(p,'input_img');
    addOptional(p,'ref_round',defaultRef);
    addOptional(p,'method',defaultMethod);
    addParameter(p,'Iterations',defaultIter);
    addParameter(p,'AccumulatedFieldSmoothing',defaultAFS);

    parse(p, input_img, varargin{:});
    
    p.Results
    % Function 
    output_img = {};
    output_img{p.Results.ref_round} = input_img{p.Results.ref_round};

    ref_og = input_img{p.Results.ref_round};
    
    if p.Results.method == "max"
        ref_og = max(ref_og, [], 4);
    else
        ref_og = sum(ref_og / 4, 4);
    end

    Nround = numel(input_img);
    for r=1:Nround
        if r ~= p.Results.ref_round
            
            fprintf(sprintf("Round %d vs. Round %d...\n", r, p.Results.ref_round));
            
            curr_og = input_img{r};
            pyd_level = floor(log2(size(curr_og, 3))); 
            
            if p.Results.method == "max"
                max_og = max(curr_og, [], 4);
            else
                max_og = sum(curr_og / 4, 4);
            end
            
            % Non-rigid registration
            [D, moving_reg] = imregdemons(max_og, ref_og, p.Results.Iterations, 'PyramidLevels', pyd_level, 'AccumulatedFieldSmoothing', p.Results.AccumulatedFieldSmoothing, 'DisplayWaitbar', false);

            % Apply displacement field on each channel
            curr_reg = zeros(size(curr_og));
            for c=1:4
                curr_ch_reg = imwarp(curr_og(:,:,:,c), D);
                curr_reg(:,:,:,c) = curr_ch_reg;
            end

            output_img{r} = uint8(curr_reg);
        end
    end 


end


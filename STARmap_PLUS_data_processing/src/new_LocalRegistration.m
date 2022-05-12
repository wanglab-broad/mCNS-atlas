function obj = new_LocalRegistration( obj, varargin )
%LocalRegistration is used to do local (non-rigid) registration for images
%in all rounds
%   -----IO-----
%   registeredImages: input mat with image stacks for all rounds 

    % Input parser
    p = inputParser;
    % Defaults
    defaultRef = 1;
    defaultMethod = "max";
    defaultIter = 60;
    defaultAFS = 1;

    addOptional(p,'ref_round',defaultRef);
    addOptional(p,'Method',defaultMethod);
    addParameter(p,'Iterations',defaultIter);
    addParameter(p,'AccumulatedFieldSmoothing',defaultAFS);

    parse(p, varargin{:});
    
    %p.Results

    % Function 
    if p.Results.Method == "max"
        ref_og = max(obj.gpuImages(:,:,:,:,p.Results.ref_round), [], 4);
    else
        ref_og = sum(obj.gpuImages(:,:,:,:,p.Results.ref_round) / 4, 4);
    end

    for r=1:obj.Nround
        tic
        if r ~= p.Results.ref_round
            
            fprintf(sprintf("Round %d vs. Round %d...", r, p.Results.ref_round));
            
            if p.Results.Method == "max"
                curr_og = max(obj.gpuImages(:,:,:,:,r), [], 4);
            else
                curr_og = sum(obj.gpuImages(:,:,:,:,r) / 4, 4);
            end
            
            pyd_level = floor(log2(obj.dimZ)); 
            
            % Non-rigid registration
            [D, ~] = imregdemons(curr_og, ref_og, p.Results.Iterations, ...
                'PyramidLevels', pyd_level, ...
                'AccumulatedFieldSmoothing', p.Results.AccumulatedFieldSmoothing, ...
                'DisplayWaitbar', false);

            if obj.useGPU
                d = gather(D);
            else
                d = D;
            end
            
            % Apply displacement field on each channel
            for c=1:obj.Nchannel
                obj.registeredImages(:,:,:,c,r) = imwarp(obj.registeredImages(:,:,:,c,r), d);
            end
            
        end
        fprintf(sprintf('[time = %.2f s]\n', toc));
    end 

    if obj.useGPU
        obj.gpuImages = gpuArray(obj.registeredImages);
    else
        obj.gpuImages = obj.registeredImages;
    end
    
end


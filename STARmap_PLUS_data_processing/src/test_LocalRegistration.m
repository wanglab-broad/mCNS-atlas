function obj = test_LocalRegistration( obj, varargin )
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

    % Get ref round img
%     if p.Results.Method == "max"
%         ref_og = max(obj.registeredImages(:,:,:,:,p.Results.ref_round), [], 4);
%     else
%         ref_og = sum(obj.registeredImages(:,:,:,:,p.Results.ref_round) / 4, 4);
%     end
%     
%     if obj.useGPU
%         ref_og = gpuArray(ref_og);
%     end

    for r=1:obj.Nround
        tic
        if r ~= p.Results.ref_round
            
            % Get reg round img
            if p.Results.Method == "max"
                ref_og = max(obj.registeredImages(:,:,:,:,p.Results.ref_round), [], 4);
                curr_og = max(obj.registeredImages(:,:,:,:,r), [], 4);
            else
                ref_og = sum(obj.registeredImages(:,:,:,:,p.Results.ref_round) / 4, 4);
                curr_og = sum(obj.registeredImages(:,:,:,:,r) / 4, 4);
            end
            
            if obj.useGPU
                ref_og = gpuArray(ref_og);
                curr_og = gpuArray(curr_og);
            end
    
            fprintf(sprintf("Round %d vs. Round %d...", r, p.Results.ref_round));
            
            
            pyd_level = floor(log2(obj.dimZ)); 
            if pyd_level == 0
                pyd_level = 1;
            end
%             pyd_level = 5; 
            
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
        reset(gpuDevice);
        fprintf(sprintf('[time = %.2f s]\n', toc));
    end 

    
end


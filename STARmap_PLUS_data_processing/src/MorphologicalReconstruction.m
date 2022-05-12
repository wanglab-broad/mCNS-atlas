function [ output_img ] = MorphologicalReconstruction( input_img, varargin )
%MorphologicalReconstruction

    % Input parser
    p = inputParser;
    % Defaults
    defaultRadius = 6;
    defaultHeight = 3;


    addRequired(p,'input_img');
    addOptional(p,'radius',defaultRadius);
    addOptional(p,'height',defaultHeight);


    parse(p, input_img, varargin{:});

    p.Results
    
    Nround = numel(input_img);
    output_img = cell(Nround, 1);
    ms = offsetstrel('ball', p.Results.radius, p.Results.height);
    se = strel('sphere',2);
    
    for r=1:Nround 
        fprintf(sprintf("Morphological Reconstruction in Round %d...\n", r));
        curr_round = input_img{r};
        output_stack = zeros(size(curr_round));
        Nchannel = size(curr_round, 4);
        for c=1:Nchannel
            curr_channel = curr_round(:,:,:,c);
            marker = imerode(curr_channel, ms);
            obr = imreconstruct(marker, curr_channel);
            curr_out = curr_channel - obr;
            mask = imbinarize(curr_out,0.06);
            %curr_out(~mask) = 0;
            
            bw = imopen(mask, se);
            curr_out(~bw) = 0;
            
            curr_out = imsubtract(imadd(curr_out, imtophat(curr_out, ms)),imbothat(curr_out, ms));
            
            output_stack(:,:,:,c) = curr_out;
        end
        output_img{r} = uint8(output_stack);
    end
    
    
end


function input_img = new_EqualizeHist3D( input_img, method )
% 

Nround = size(input_img, 5);

switch method
    case "intra_round"
        for t=1:Nround 
            fprintf('Equalizing round %d\n', t);
            currStack = input_img(:,:,:,:,t);   

            for i=1:4                
                currStack(:,:,:,i) = uint8(imhistmatchn(uint8(currStack(:,:,:,i)), uint8(currStack(:,:,:,1)), 256));
                %currStack(:,:,:,i) = gather(histeq(gpuArray(currStack(:,:,:,i))));
            end    

            input_img(:,:,:,:,t) = currStack;
        end
        
    case "inter_round"
        for t=1:4 
            fprintf('Equalizing channel %d\n', t);
            currStack = input_img(:,:,:,t,:);   

            for i=1:Nround               
                % currStack = input_img{i}(:,:,:,t);   
                % input_img{i}(:,:,:,t) = uint8(imhistmatchn(uint8(currStack), uint8(input_img{1}(:,:,:,t)), 256));
                currStack(:,:,:,:,i) = uint8(imhistmatchn(uint8(currStack(:,:,:,:,i)), uint8(currStack(:,:,:,:,1)), 256));
            end    

            input_img(:,:,:,t,:) = currStack;
        end
        
end

end


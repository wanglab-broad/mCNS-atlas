function FinalImage = new_LoadMultipageTiff( fname, formatIn, formatOut, useGPU )

    % Suppress all warnings 
    warning('off','all');
    
    if nargin < 2
        formatIn = 'uint8';
        formatOut = 'uint8';
    end

    InfoImage=imfinfo(fname);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);

 
    if useGPU
        %FinalImage=zeros(nImage,mImage,NumberImages,format, 'gpuArray');
        FinalImage=zeros(nImage ,mImage ,NumberImages, formatIn, 'gpuArray');
    else
        %FinalImage=zeros(nImage,mImage,NumberImages,format);
        FinalImage=zeros(nImage, mImage, NumberImages, formatIn);
    end

    TifLink = Tiff(fname, 'r');
    for i=1:NumberImages
       TifLink.setDirectory(i);
       FinalImage(:,:,i)=TifLink.read();
    end
    
    if ~strcmp(formatIn, formatOut)
        % Convert to uint8
        if strcmp(formatOut, 'uint8')
            FinalImage = im2uint8(FinalImage);
        end
    end
    
    TifLink.close();
    
end


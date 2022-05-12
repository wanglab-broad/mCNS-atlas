function FinalImage = LoadMultipageTiff( fname, format, useGPU )

    % Suppress all warnings 
    warning('off','all');
    
    if nargin < 2
        format = 'uint8';
    end
    
    if nargin < 3
        useGPU = false;
    end

    InfoImage=imfinfo(fname);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);

    if useGPU
        FinalImage=zeros(nImage,mImage,NumberImages,format, 'gpuArray');
    else
        FinalImage=zeros(nImage,mImage,NumberImages,format);
    end

    TifLink = Tiff(fname, 'r');
    for i=1:NumberImages
       TifLink.setDirectory(i);
       FinalImage(:,:,i)=TifLink.read();
    end
    TifLink.close();
    
end


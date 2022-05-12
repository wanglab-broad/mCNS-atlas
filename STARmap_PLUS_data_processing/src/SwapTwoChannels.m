function input_img = SwapTwoChannels( input_img, channel_1, channel_2 )
%SwapTwoChannels

%     Nround = size(input_img, 5);
%     
%     for r=1:Nround
%         % temp = input_img{r};
%         temp = input_img(:,:,:,:,r);
%         swap_1 = temp(:,:,:,channel_1);
%         swap_2 = temp(:,:,:,channel_2);
%         temp(:,:,:,channel_2) = swap_1;
%         temp(:,:,:,channel_1) = swap_2;
%         input_img(:,:,:,:,r) = temp;
%         % input_img{r} = temp;
%         
%     end

    temp = input_img(:,:,:,channel_1,:);
    input_img(:,:,:,channel_1,:) = input_img(:,:,:,channel_2,:);
    input_img(:,:,:,channel_2,:) = temp;

end

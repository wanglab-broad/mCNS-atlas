% to use Miji you need to put the scripts directory of the Fiji app in
% Matlab search path and you also need to increase the memory of java
% instance in matlab(go to preference/general/java heap memory)
tic
clear all
file_path = "Data/STARmap_data/160genes/mPFC/seg_test.tif";
I = LoadMultipageTiff(file_path, 'uint8', false); % load the image here 
%% preprocess raw image
parfor i = 1:27
 B(:,:,i) = seg(I(:,:,i),3,4,400);
end
%% write binary image
B = double(B)*255;
for j=1:27        
 imwrite(B(:,:,j), "seg_test_preproessed.tif", 'writemode', 'append');        
end
%% use Miji for 3D distance watershed
Miji;
MIJ.run('Open...', 'path = [seg_test_preproessed.tif]'); % in [] put the path of the binary image of B just saved
MIJ.run("Distance Transform Watershed 3D", "distances=[Borgefors (3,4,5)] output=[16 bits] normalize dynamic=5.5 connectivity=6")
I2 = MIJ.getCurrentImage; % I2 Return the segmented label of the image
MIJ.run("Set Label Map", "colormap=[Golden angle] background=Black shuffle");
MIJ.run("Clear Results")
MIJ.exit
toc
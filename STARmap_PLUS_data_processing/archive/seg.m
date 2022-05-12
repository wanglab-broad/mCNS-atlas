function B = seg(A,p,q,r)
A = imbinarize(A,0.1);
se = strel('disk',p);
binary2 = imfill(A,q,'holes');
binary3 = imclose(binary2,se);
binary4 = imfill(binary3,q,'holes');
B = bwareaopen(binary4,r);
end
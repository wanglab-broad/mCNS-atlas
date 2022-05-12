function spot_img = get_spot_img(input_img, points, esize)
%get_spot_stak 
% esize is [x,y,z] vector

    if nargin < 3
        esize = [1 1];
    end

    points = int16(points);
    [X ,Y] = size(input_img);
    Npoints = size(points,1);
    spot_img = zeros(X, Y);
    

    for i=1:Npoints
        
        p = points(i,:);
        extentsX = GetExtents(p(2), esize(1), X);
        extentsY = GetExtents(p(1), esize(2), Y);                    

        % 3x3x3
        spot_img(extentsX, extentsY) = 255;

    end
    
    spot_img = uint8(spot_img);


    function e = GetExtents(pos, esize, lim)

    if pos-esize < 1 
        e1 = 1;
    else
        e1 = pos-esize;
    end

    if pos+esize > lim
        e2 = lim;
    else
        e2 = pos+esize;
    end

    e = e1:e2;


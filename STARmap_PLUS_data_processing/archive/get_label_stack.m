function spot_stack = get_label_stack(input_size, points, esize)
%get_spot_stak 
% esize is [x,y,z] vector

    if nargin < 3
        esize = [1 1 1];
    end

    points = int32(points);
    X = input_size(1);
    Y = input_size(2);
    Z = input_size(3);
    Npoints = size(points,1);
    spot_stack = zeros(X, Y, Z);
    

    for i=1:Npoints
        
        p = points(i,:);

        extentsX = GetExtents(p(2), esize(1), X);
        extentsY = GetExtents(p(1), esize(2), Y);                    
        % extentsZ = GetExtents(p(3), esize(3), Z);

        % 3x3x3
        spot_stack(extentsX, extentsY) = points(i,3);

    end
    
    spot_stack = uint16(spot_stack);


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


function e = GetExtents(pos, ksize, lim)

    if pos-ksize < 1 
        e1 = 1;
    else
        e1 = pos-ksize;
    end

    if pos+ksize > lim
        e2 = lim;
    else
        e2 = pos+ksize;
    end

    e = e1:e2;
end


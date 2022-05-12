function readsLocation = GetReadsLocation( goodSpots, labelImage )
%GetReadsLocation

    Nreads = size(goodSpots, 1);
    readsLocation = zeros(Nreads, 1);
    
    for i = 1:Nreads
        curr_loc = labelImage(goodSpots(i, 2), goodSpots(i, 1), goodSpots(i, 3)); 
        readsLocation(i) = curr_loc;
    end
    
end


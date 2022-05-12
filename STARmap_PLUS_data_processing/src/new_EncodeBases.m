function colorSeq = new_EncodeBases( seq )
% new_EncodeSOLID

    % construct hash table for encoding
    k = {'AT','CT','GT','TT',...
        'AG', 'CG', 'GG', 'TG',...
        'AC', 'CC', 'GC', 'TC',...
        'AA', 'CA', 'GA', 'TA'};
    v = {4,3,2,1,3,4,1,2,2,1,4,3,1,2,3,4};
    coding = containers.Map(k,v);
    
    start = 1;
    back = start + 1;
    colorSeq = "";
    while back <= strlength(seq)
        curr_str = extractBetween(seq, start, back);
        colorSeq = colorSeq + coding(curr_str);
        start = start + 1;
        back = start + 1;
    end
    
end


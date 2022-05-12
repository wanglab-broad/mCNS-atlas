function [color_seq] = encode_solid(seq, to_str)
%ENCODE_SOLID is used to convert barcode sequence of each gene in codebook
%to color strings 
%   -----IO-----
%   seq = barcode sequence for each gene, ie: CATCATC
%   color_seq = 243243 or '243243'
    
    if nargin < 2
        to_str = true;
    end
    
    % Base table
    k = {'AA','CA','GA','TA',...
        'AC', 'CC', 'GC', 'TC',...
        'AG', 'CG', 'GG', 'TG',...
        'AT', 'CT', 'GT', 'TT'};
    
    % Color table
    v = {1,2,3,4,...
        2,1,4,3,...
        3,4,1,2,...
        4,3,2,1};
    
    coding = containers.Map(k,v);
    color_seq = [];
    
    for i=1:(numel(seq) - 1)
        s = seq(i : i+1);
        color_seq = [color_seq coding(s)];
    end
    
    if to_str
        s = '';
        for i=1:numel(color_seq)
            s = [s num2str(color_seq(i))];
        end
        color_seq = s;
    end

end


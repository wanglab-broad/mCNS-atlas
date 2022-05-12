function bases = new_DecodeCS( basecsMat, startBase )
%new_DecodeCS

    % get dims
    [Npoint, Nround] = size(basecsMat);
    
    % preallocating
    bases = cell(Npoint,1);
    
    % construct reverse hash table for decoding
    %prb = textprogressbar(Npoint);
    k = {1,2,3,4};
    v = {{'AA','CC','GG','TT'},...
        {'AC','CA','GT','TG'},....
        {'AG','CT','GA','TC'},....
        {'AT','CG','GC','TA'}};

    reverseMap = containers.Map(k,v);     

%     parfor i=1:Npoint
%         currSeq = basecsMat(i,:);
%         currStr = '';
%         possible = reverseMap(currSeq(1));
%         for j=1:4
%             p = possible{j};
%             if strcmp(p(1), startBase)
%                 currStr = p;
%             end
%         end
%         for kk=2:numel(currSeq)
%             possible = reverseMap(currSeq(kk));
%             for j=1:4
%                 p = possible{j};
%                 if strcmp(p(1), currStr(kk))
%                     currStr = [currStr p(2)];
%                 end
%             end
%         end
%         bases{i} = currStr;
% 
%     end
    
    
    for i=1:Npoint
        currSeq = basecsMat(i,:);
        currStr = '';
        ref_base = startBase;
        for j=1:numel(currSeq)
            % currSeq(j)
            possible = reverseMap(currSeq(j));
            for n=1:4
                p = possible{n};
                if strcmp(p(1), ref_base) && isempty(currStr)
                    currStr = p;
                elseif strcmp(p(1), ref_base)
                    currStr = [currStr p(2)];
                end
            end 
            % currStr
            ref_base = currStr(end);
            
        end
        
        bases{i} = currStr;
        
    end
    
    %{
    function p = PossibleBases(c)
    switch c
        case 1
            p = {'AA', 'CC', 'GG', 'TT'};
        case 2
            p = {'AC', 'CA', 'GT', 'TG'};
        case 3
            p = {'AG', 'CT', 'GA', 'TC'};
        case 4
            p = {'AT', 'CG', 'GC', 'TA'};
    end
    %}
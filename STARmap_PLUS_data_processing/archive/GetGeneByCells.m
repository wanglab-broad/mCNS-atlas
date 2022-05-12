function count_matrix = GetGeneByCells( obj )
%GetGeneByCells

    count_matrix = table(convertCharsToStrings(obj.barcodeNames'), 'VariableNames', "Gene");
    
    for n = 1:obj.Ncells
        curr_idx = obj.goodReadsLoc == n;
        curr_reads = obj.goodReads(curr_idx);
        curr_reads = cellfun(@(x) obj.seqToGene(x), curr_reads, 'UniformOutput', false);
        curr_reads = convertCharsToStrings(curr_reads);
        
        [C,~,ic] = unique(curr_reads);
        temp_counts = accumarray(ic,1);
        %curr_counts = [C, temp_counts];
        curr_name = sprintf('Cell_%d', n);
        curr_table = table(C, temp_counts, 'VariableNames', {'Gene', curr_name});
        %curr_table(isnan(curr_table)) = 0;
        %curr_table = table(curr_counts(:,1), curr_counts(:,2), 'VariableNames', {'Gene', curr_name});
        count_matrix = outerjoin(count_matrix, curr_table, 'Type', 'left', 'MergeKeys', true);
        
    end
    
    count_matrix = fillmissing(count_matrix,'constant',0,'DataVariables',@isnumeric);
end


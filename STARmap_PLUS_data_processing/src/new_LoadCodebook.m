function [ geneToSeq, seqToGene ] = new_LoadCodebook( inputPath, remove_index, doReverse )
% new_LoadCodebook

    % load file
    fname = fullfile(inputPath, 'genes.csv');
    f = readmatrix(fname, 'OutputType', 'string', "Delimiter", ',');

    % load gene name and sequence 
    % f(:,1) - gene name, f(:,2) - gene barcode 
    if doReverse
        f(:,2) = reverse(f(:,2));
    end

    for i=1:size(f, 1)
        f(i,2) = new_EncodeBases(f(i,2));
    end
    
    if ~isempty(remove_index)
        f(:,2) = eraseBetween(f(:,2), remove_index, remove_index);

        % flip
        front = extractAfter(f(:,2), remove_index-1);
        back = extractBefore(f(:,2), remove_index);
        f(:,2) = front + back;

    end
    

    
    seqToGene = containers.Map(f(:,2), f(:,1));
    geneToSeq = containers.Map(f(:,1), f(:,2));


function [ geneToSeq, seqToGene ] = new_LoadCodebookDuo( inputPath, doReverse )
% new_LoadCodebook


    fname = fullfile(inputPath, 'genes_duo.csv');
    f = readmatrix(fname, 'OutputType', 'string', "Delimiter", ',');

    % load gene name and sequence 

    if doReverse
        f(:,2) = reverse(f(:,2));
        f(:,3) = reverse(f(:,3));
    end

    for i=1:size(f, 1)
        f(i,2) = new_EncodeBases(f(i,2));
        f(i,3) = new_EncodeBases(f(i,3));
    end


    seqToGene = {containers.Map(f(:,2), f(:,1)), containers.Map(f(:,3), f(:,1))};
    geneToSeq = {containers.Map(f(:,1), f(:,2)), containers.Map(f(:,1), f(:,3))};


function counts = SaveGoodReads( obj, input_id )
%

    counts = containers.Map(obj.barcodeNames, zeros(numel(obj.barcodeNames),1));
    for i=1:numel(obj.goodReads)               
        counts(obj.seqToGene(obj.goodReads{i})) = counts(obj.seqToGene(obj.goodReads{i})) + 1;
    end
    
    fname = strcat('gene_counts', input_id, '.csv');
    fid = fopen(fullfile(obj.outputPath, fname),'w');
    for j=counts.keys
        fwrite(fid, sprintf('%s,%s,%d\n', j{1}, obj.geneToSeq(j{1}), counts(j{1})));
    end
    fclose(fid);


end


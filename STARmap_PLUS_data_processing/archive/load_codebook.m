function [gene_to_seq_dict, seq_to_gene_dict] = load_codebook(input_dir, flip)
%LOAD_CODEBOOK is used to load Codebook (genes.csv) from the
%input directory and convert it into two dictionaries
%   -----IO-----
%   input_dir = path/to/genes.csv
%   flip = true/false Flip array left to rightcollapse, default: true
%   gene_to_seq_dict = dictionary where genes are keys and color_seq are
%   values 
%   seq_to_gene_dict = dictionary where color_seq are keys and genes are
%   values 

    if nargin < 2
        flip = true;
    end
    
    genes = {};
    seqs = {};
    fname = fullfile(input_dir, 'genes.csv');
    f = fopen(fname,'r');
    f_line = fgets(f);
    fields = strsplit(f_line,',');

    % Get gene
    genes{end+1} = fields{1}; 
    
    % Get sequence
    if flip
        seqs{end+1} = fliplr(strtrim(fields{2}));
    else
        seqs{end+1} = strtrim(fields{2});
    end

    
    while ischar(f_line)
        f_line = fgets(f);    
        if f_line ~= -1
            fields = strsplit(f_line,',');
            if isempty(fields{1}) || isempty(fields{2})
                break;
            else
                genes{end+1} = fields{1};
                if flip
                    seqs{end+1} = fliplr(strtrim(fields{2}));
                else
                    seqs{end+1} = strtrim(fields{2});
                end
            end
        end    
    end

    fclose(f);

    % Get color string for each barcode 
    for i=1:numel(seqs)
        seqs{i} = encode_solid(seqs{i});
    end
    
    gene_to_seq_dict = containers.Map(genes, seqs);
    seq_to_gene_dict = containers.Map(seqs, genes);

end


function [bases, maxColors, baseScores] = new_GetBaseSeq( colorSeq )
%GetBaseSeq

    [Npoint, Nround, ~] = size(colorSeq);
    maxColors = zeros(Npoint, Nround);
    baseScores = zeros(Npoint, Nround);

    % determine maximal color for each round
    fprintf('Geting max color...\n');

    count = 0;
    for i=1:Npoint
        multi_max = false;
        for j=1:Nround
            currMax = max(colorSeq(i,j,:), [], 3);
            if ~isnan(currMax)
                m = find(colorSeq(i,j,:) == currMax);
%                 maxColors(i,j) = m(1);
%                 baseScores(i,j) = -log(currMax);
                if numel(m) ~= 1
                    multi_max = true;
                    %count = count + 1;
                    maxColors(i,j) = -1;
                    baseScores(i,j) = Inf;
                else
                    maxColors(i,j) = m(1);
                    baseScores(i,j) = -log(currMax);
                end
            else
                maxColors(i,j) = -1;
                baseScores(i,j) = Inf;
            end
        end
        
        if multi_max
            count = count + 1;
        end
    end
    fprintf('%d \n', count);
    
    
    % alexa decoding
    bases = {Npoint, 1};
    fprintf('Decoding...\n');
    parfor i=1:Npoint
        if ~any(isinf(baseScores(i, :)))
            bases{i, 1} = Colorseq2Str(maxColors(i, :));
        end
    end
    fprintf('\n');


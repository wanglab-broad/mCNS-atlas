function obj = new_FilterReads_Duo( obj, endBases, split_loc, q_score_thers, showPlots )
% FilterReads_Duo

%     if nargin < 4
%         correctErrors = false;
%     end

%     % remove reads with any infinite values
%     finiteScores = ~any(isinf(obj.allScores), 2);            

    % remove reads with bad quality scores
    if showPlots
        figure(1);
        histogram(mean(obj.allScores, 2), 100)
        xlabel('Average scores'); ylabel('Count');
    end
    
    % export_fig is an add-on
    % export_fig(fullfile(obj.outPath, 'average_scores.png'));


    belowScoreThresh = mean(obj.allScores, 2) < q_score_thers; % 0.5; 
    s = sprintf('%f [%d / %d] percent of reads are below score thresh %d\n',...
        sum(belowScoreThresh)/numel(belowScoreThresh),...
        sum(belowScoreThresh), ...
        numel(belowScoreThresh), ...
        q_score_thers);
    fprintf(s);
    
    if ~isempty(obj.log)
        fprintf(obj.log, s);
    end
    
    obj.allReads = obj.allReads(belowScoreThresh);
    obj.allScores = obj.allScores(belowScoreThresh,:);
    obj.allSpots = obj.allSpots(belowScoreThresh,:);
    obj.basecsMat = obj.basecsMat(belowScoreThresh,:);
    
    if showPlots
        figure(1);
        errorbar(mean(obj.allScores), std(obj.allScores),'ko-'); 
        xlim([0 obj.Nround+1]); 
        xlabel('Round'); ylabel('Average qual score');
    end
    
% =========
    % Filter reads by whether they are in the codebook
    % This is just as a sanity check; reads are actually filtered
    % by whether they are in the codebook
    filtBases = obj.basecsMat;
    correctSeqs = zeros(size(filtBases,1), 2); % count sequences that are of correct form

    cs_front = filtBases(:, 1:split_loc);
    cs_back = filtBases(:, split_loc+1:obj.Nround);

    bases_front = new_DecodeCS(cs_front, endBases(1));
    bases_back = new_DecodeCS(cs_back, endBases(2));

    for i=1:numel(bases_front)
        currSeq = bases_front{i};
        if currSeq(1) == endBases(1) && currSeq(end) == endBases(1)
            correctSeqs(i, 1) = 1;
        end
    end

    for i=1:numel(bases_back)
        currSeq = bases_back{i};
        if currSeq(1) == endBases(2) && currSeq(end) == endBases(2)
            correctSeqs(i, 2) = 1;
        end
    end

    correctSeqs = correctSeqs(:,1) .* correctSeqs(:,2);

    fprintf('Filtration Statistics:\n');
    % fprintf(obj.log, 'Filtration Statistics:\n');
    score_1 = sum(correctSeqs)/size(filtBases,1);
    s1 = sprintf('%f [%d / %d] percent of good reads match barcode pattern %sN..N%s - %sN..N%s\n',...
        sum(correctSeqs)/size(filtBases,1),...
        sum(correctSeqs),...
        size(filtBases,1),...
        endBases(1),...
        endBases(1),...
        endBases(2),...
        endBases(2));
    fprintf(s1);

    % filter reads based on codebook
    Nreads = numel(obj.allReads);
    inCodebook = zeros(Nreads, 1);
    codebookSeqs = obj.barcodeSeqs;

    for s=1:Nreads
        str = obj.allReads{s};
        if ismember(str, codebookSeqs)         
            inCodebook(s) = 1;
        else
            inCodebook(s) = 0;
        end
    end

    % sum(inCodebook)

    readsToKeep = inCodebook==1;
    obj.goodSpots = obj.allSpots(readsToKeep,:);
    obj.goodReads = obj.allReads(readsToKeep);
    obj.goodScores = obj.allScores(readsToKeep,:);

    if showPlots
        figure(2);
        errorbar(mean(obj.goodScores), std(obj.goodScores),'ko-'); 
        xlim([0 obj.Nround+1]); 
        xlabel('Round'); ylabel('Average qual score');
    end

    score_2 = sum(readsToKeep)/Nreads;
    s2 = sprintf('%f [%d / %d] percent of good reads are in codebook\n',...
        sum(readsToKeep)/Nreads,...
        sum(readsToKeep),...
        Nreads);
    fprintf(s2);

    score_3 = sum(readsToKeep)/sum(correctSeqs);
    s3 = sprintf('%f [%d / %d] percent of form matched reads are in codebook\n',...
        sum(readsToKeep)/sum(correctSeqs),...
        sum(readsToKeep), ...
        sum(correctSeqs));            
    fprintf(s3);

    obj.FilterScores = [score_1 score_2 score_3]; 
    
    if ~isempty(obj.log)
        fprintf(obj.log, s1);
        fprintf(obj.log, s2);
        fprintf(obj.log, s3);
    end
        
end
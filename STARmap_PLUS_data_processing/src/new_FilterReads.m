function obj = new_FilterReads( obj, endBases, q_score_thers, showPlots )
% FilterReads
  
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

% ===========

    % Filter reads by whether they are in the codebook
    % This is just as a sanity check; reads are actually filtered
    % by whether they are in the codebook
    filtBases = obj.basecsMat;
    correctSeqs = zeros(size(filtBases, 1), 1); % count sequences that are of correct form

    
    % older code for filtering reads in sequence space
    bases = new_DecodeCS(filtBases, endBases(1));
    for i=1:numel(bases)
        currSeq = bases{i};
        if currSeq(1) == endBases(1) && currSeq(end) == endBases(2)
            correctSeqs(i) = 1;
        end
    end
    
    fprintf('Filtration Statistics:\n');
    fprintf(obj.log, 'Filtration Statistics:\n');
    score_1 = sum(correctSeqs)/size(filtBases,1);
    s = sprintf('%f [%d / %d] percent of good reads are %sNNNNN%s\n',...
        sum(correctSeqs)/size(filtBases,1),...
        sum(correctSeqs),...
        size(filtBases,1),...
        endBases(1),...
        endBases(2));
    fprintf(s);
    fprintf(obj.log, s);


    % filter reads based on codebook
    Nreads = numel(obj.allReads);
    inCodebook = zeros(Nreads,1);
    codebookSeqs = obj.barcodeSeqs;
    for s=1:Nreads
        str = obj.allReads{s};
        if ismember(str, codebookSeqs)         
            inCodebook(s) = 1;
        else
            inCodebook(s) = 0;
        end
    end


    % correct things in colospace
%     if correctErrors
%         Ncorrected = 0;            
%         upd = textprogressbar(sum(inCodebook==0));
%         correctionCount = 0;
% 
%         qualScores = obj.allScores;
%         basecsMat = obj.basecsMat;
%         barcodeMat = obj.barcodeMat;
%         Nround = obj.Nround;
% 
%         correctedIdx = zeros(Nreads,1);
%         for s=1:Nreads                   
%             if inCodebook(s) == 0
%                 correctionCount = correctionCount + 1;
%                 if mod(s,100) == 0
%                     upd(correctionCount);
%                 end
% 
%                 currBase = basecsMat(s,:);
% 
%                 dists = pdist2(currBase, barcodeMat, 'hamming')*Nround; % calculate hamming distance between currBase and each barcode
%                 minDist = min(dists);
%                 if minDist == 1
%                     potentialMatches = find(dists == minDist);
%                     if numel(potentialMatches) > 1 % if multiple matches, choose the one with the lowest quality score for the mismatched position
%                         %{
%                         currQuals = qualScores(s,:);
%                         mismatchedPositionScores = zeros(numel(potentialMatches),1);       
%                         for matchIdx=1:numel(potentialMatches)
%                             currBarcode = barcodeCS(potentialMatches(matchIdx),:);
%                             for baseIdx=1:numel(currBarcode)
%                                 if currBarcode(baseIdx) ~= currBase(baseIdx)
%                                     mismatchedPositionScores(matchIdx) = currQuals(baseIdx);
%                                 end
%                             end
%                             % compare 
%                         end
%                         correctedIdx(s) = potentialMatches(mismatchedPositionScores == min(mismatchedPositionScores));
%                         inCodebook(s) = 1;
%                         %}
%                     else
%                         correctedIdx(s) = potentialMatches; % look up correct seq                                                                        
%                         inCodebook(s) = 1;
%                     end                
%                 end
%             end
%         end
%         
%         Ncorrected = sum(correctedIdx > 0);
%         disp(Ncorrected);
% 
%         % actually reassign bases
%         reassigned = find(inCodebook(s) == 0 & correctedIdx(s) > 0);
%         for s=reassigned                    
%             obj.bases{s} = obj.barcodeSeqs(correctedIdx(s));                    
%         end
%     end
%     
    % only save reads that are below the score thresh and are in
    % the codebook
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
    s = sprintf('%f [%d / %d] percent of good reads are in codebook\n',...
        sum(readsToKeep)/Nreads,...
        sum(readsToKeep),...
        Nreads);
    fprintf(s);
    fprintf(obj.log, s);

%     if correctErrors
%         s = sprintf('%f [%d / %d] reads in codebook were corrected\n', ...
%             Ncorrected/sum(inCodebook), ...
%             Ncorrected, ...
%             sum(inCodebook));
%         fprintf(s); fprintf(fid, s);
%     end

    score_3 = sum(readsToKeep)/sum(correctSeqs);
    s = sprintf('%f [%d / %d] percent of %sNNNNN%s reads are in codebook\n',...
        sum(readsToKeep)/sum(correctSeqs),...
        sum(readsToKeep), ...
        sum(correctSeqs),...
        endBases(1),...
        endBases(2));            
    fprintf(s);
    fprintf(obj.log, s);

    obj.FilterScores = [score_1 score_2 score_3]; 
    
end